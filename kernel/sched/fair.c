// SPDX-License-Identifier: GPL-2.0
/*
 * Completely Fair Scheduling (CFS) Class (SCHED_NORMAL/SCHED_BATCH)
 *
 *  Copyright (C) 2007 Red Hat, Inc., Ingo Molnar <mingo@redhat.com>
 *
 *  Interactivity improvements by Mike Galbraith
 *  (C) 2007 Mike Galbraith <efault@gmx.de>
 *
 *  Various enhancements by Dmitry Adamushko.
 *  (C) 2007 Dmitry Adamushko <dmitry.adamushko@gmail.com>
 *
 *  Group scheduling enhancements by Srivatsa Vaddagiri
 *  Copyright IBM Corporation, 2007
 *  Author: Srivatsa Vaddagiri <vatsa@linux.vnet.ibm.com>
 *
 *  Scaled math optimizations by Thomas Gleixner
 *  Copyright (C) 2007, Thomas Gleixner <tglx@linutronix.de>
 *
 *  Adaptive scheduling granularity, math enhancements by Peter Zijlstra
 *  Copyright (C) 2007 Red Hat, Inc., Peter Zijlstra
 */
#include <linux/rbtree_augmented.h>
#include "sched.h"

#include <trace/events/sched.h>

#include "walt.h"

#ifdef CONFIG_SMP
static inline bool task_fits_max(struct task_struct *p, int cpu);
static inline unsigned long boosted_task_util(struct task_struct *task);
#endif /* CONFIG_SMP */

#ifdef CONFIG_SCHED_WALT
static void walt_fixup_sched_stats_fair(struct rq *rq, struct task_struct *p,
					u16 updated_demand_scaled,
					u16 updated_pred_demand_scaled);
static void walt_fixup_nr_big_tasks(struct rq *rq, struct task_struct *p,
					int delta, bool inc);
#endif /* CONFIG_SCHED_WALT */

#if defined(CONFIG_SCHED_WALT) && defined(CONFIG_CFS_BANDWIDTH)

static void walt_init_cfs_rq_stats(struct cfs_rq *cfs_rq);
static void walt_inc_cfs_rq_stats(struct cfs_rq *cfs_rq,
				  struct task_struct *p);
static void walt_dec_cfs_rq_stats(struct cfs_rq *cfs_rq,
				  struct task_struct *p);
static void walt_inc_throttled_cfs_rq_stats(struct walt_sched_stats *stats,
					    struct cfs_rq *cfs_rq);
static void walt_dec_throttled_cfs_rq_stats(struct walt_sched_stats *stats,
					    struct cfs_rq *cfs_rq);
#else
static inline void walt_init_cfs_rq_stats(struct cfs_rq *cfs_rq) {}
static inline void
walt_inc_cfs_rq_stats(struct cfs_rq *cfs_rq, struct task_struct *p) {}
static inline void
walt_dec_cfs_rq_stats(struct cfs_rq *cfs_rq, struct task_struct *p) {}

#define walt_inc_throttled_cfs_rq_stats(...)
#define walt_dec_throttled_cfs_rq_stats(...)

#endif

/*
 * Enable/disable honoring sync flag in energy-aware wakeups.
 */
unsigned int sysctl_sched_sync_hint_enable = 1;

/*
 * Enable/disable using cstate knowledge in idle sibling selection
 */
unsigned int sysctl_sched_cstate_aware = 1;

/*
 * The initial- and re-scaling of tunables is configurable
 *
 * Options are:
 *
 *   SCHED_TUNABLESCALING_NONE - unscaled, always *1
 *   SCHED_TUNABLESCALING_LOG - scaled logarithmical, *1+ilog(ncpus)
 *   SCHED_TUNABLESCALING_LINEAR - scaled linear, *ncpus
 *
 * (default SCHED_TUNABLESCALING_LOG = *(1+ilog(ncpus))
 */
enum sched_tunable_scaling sysctl_sched_tunable_scaling = SCHED_TUNABLESCALING_LOG;

/*
 * Minimal preemption granularity for CPU-bound tasks:
 *
 * (default: 3 msec * (1 + ilog(ncpus)), units: nanoseconds)
 */
unsigned int sysctl_sched_base_slice			= 3000000ULL;
static unsigned int normalized_sysctl_sched_base_slice	= 3000000ULL;

const_debug unsigned int sysctl_sched_migration_cost	= 500000UL;
DEFINE_PER_CPU_READ_MOSTLY(int, sched_load_boost);

#ifdef CONFIG_SMP
/*
 * For asym packing, by default the lower numbered CPU has higher priority.
 */
int __weak arch_asym_cpu_priority(int cpu)
{
	return -cpu;
}

/*
 * The margin used when comparing CPU capacities.
 * is 'cap1' noticeably greater than 'cap2'
 *
 * (default: ~5%)
 */
#define capacity_greater(cap1, cap2) ((cap1) * 1024 > (cap2) * 1078)

/*
 * The margin used when comparing utilization with CPU capacity.
 *
 * (default: ~20%)
 */
#define fits_capacity(cap, max)	((cap) * 1280 < (max) * 1024)

#endif

#ifdef CONFIG_CFS_BANDWIDTH
/*
 * Amount of runtime to allocate from global (tg) to local (per-cfs_rq) pool
 * each time a cfs_rq requests quota.
 *
 * Note: in the case that the slice exceeds the runtime remaining (either due
 * to consumption or the quota being specified to be smaller than the slice)
 * we will always only issue the remaining available time.
 *
 * (default: 5 msec, units: microseconds)
 */
unsigned int sysctl_sched_cfs_bandwidth_slice		= 5000UL;
#endif

#ifdef CONFIG_SCHED_WALT
unsigned int sched_capacity_margin_up[CPU_NR] = {
			[0 ... CPU_NR-1] = 1078}; /* ~5% margin */
unsigned int sched_capacity_margin_down[CPU_NR] = {
			[0 ... CPU_NR-1] = 1205}; /* ~15% margin */
/* 1ms default for 20ms window size scaled to 1024 */
unsigned int sysctl_sched_min_task_util_for_boost = 51;
/* 0.68ms default for 20ms window size scaled to 1024 */
unsigned int sysctl_sched_min_task_util_for_colocation = 35;
__read_mostly unsigned int sysctl_sched_prefer_spread;
unsigned int sysctl_walt_rtg_cfs_boost_prio = 99; /* disabled by default */
unsigned int sysctl_walt_low_latency_task_threshold; /* disabled by default */
__read_mostly unsigned int sysctl_sched_force_lb_enable = 1;
#endif

static inline void update_load_add(struct load_weight *lw, unsigned long inc)
{
	lw->weight += inc;
	lw->inv_weight = 0;
}

static inline void update_load_sub(struct load_weight *lw, unsigned long dec)
{
	lw->weight -= dec;
	lw->inv_weight = 0;
}

static inline void update_load_set(struct load_weight *lw, unsigned long w)
{
	lw->weight = w;
	lw->inv_weight = 0;
}

/*
 * Increase the granularity value when there are more CPUs,
 * because with more CPUs the 'effective latency' as visible
 * to users decreases. But the relationship is not linear,
 * so pick a second-best guess by going with the log2 of the
 * number of CPUs.
 *
 * This idea comes from the SD scheduler of Con Kolivas:
 */
static unsigned int get_update_sysctl_factor(void)
{
	unsigned int cpus = min_t(unsigned int, num_online_cpus(), 8);
	unsigned int factor;

	switch (sysctl_sched_tunable_scaling) {
	case SCHED_TUNABLESCALING_NONE:
		factor = 1;
		break;
	case SCHED_TUNABLESCALING_LINEAR:
		factor = cpus;
		break;
	case SCHED_TUNABLESCALING_LOG:
	default:
		factor = 1 + ilog2(cpus);
		break;
	}

	return factor;
}

static void update_sysctl(void)
{
	unsigned int factor = get_update_sysctl_factor();

#define SET_SYSCTL(name) \
	(sysctl_##name = (factor) * normalized_sysctl_##name)
	SET_SYSCTL(sched_base_slice);
#undef SET_SYSCTL
}

void sched_init_granularity(void)
{
	update_sysctl();
}

#define WMULT_CONST	(~0U)
#define WMULT_SHIFT	32

static void __update_inv_weight(struct load_weight *lw)
{
	unsigned long w;

	if (likely(lw->inv_weight))
		return;

	w = scale_load_down(lw->weight);

	if (BITS_PER_LONG > 32 && unlikely(w >= WMULT_CONST))
		lw->inv_weight = 1;
	else if (unlikely(!w))
		lw->inv_weight = WMULT_CONST;
	else
		lw->inv_weight = WMULT_CONST / w;
}

/*
 * delta_exec * weight / lw.weight
 *   OR
 * (delta_exec * (weight * lw->inv_weight)) >> WMULT_SHIFT
 *
 * Either weight := NICE_0_LOAD and lw \e sched_prio_to_wmult[], in which case
 * we're guaranteed shift stays positive because inv_weight is guaranteed to
 * fit 32 bits, and NICE_0_LOAD gives another 10 bits; therefore shift >= 22.
 *
 * Or, weight =< lw.weight (because lw.weight is the runqueue weight), thus
 * weight/lw.weight <= 1, and therefore our shift will also be positive.
 */
static u64 __calc_delta(u64 delta_exec, unsigned long weight, struct load_weight *lw)
{
	u64 fact = scale_load_down(weight);
	int shift = WMULT_SHIFT;

	__update_inv_weight(lw);

	if (unlikely(fact >> 32)) {
		while (fact >> 32) {
			fact >>= 1;
			shift--;
		}
	}

	/* hint to use a 32x32->64 mul */
	fact = (u64)(u32)fact * lw->inv_weight;

	while (fact >> 32) {
		fact >>= 1;
		shift--;
	}

	return mul_u64_u32_shr(delta_exec, fact, shift);
}

/*
 * delta /= w
 */
static inline u64 calc_delta_fair(u64 delta, struct sched_entity *se)
{
	if (unlikely(se->load.weight != NICE_0_LOAD))
		delta = __calc_delta(delta, NICE_0_LOAD, &se->load);

	return delta;
}

const struct sched_class fair_sched_class;

/**************************************************************
 * CFS operations on generic schedulable entities:
 */

#ifdef CONFIG_FAIR_GROUP_SCHED
static inline struct task_struct *task_of(struct sched_entity *se)
{
	SCHED_WARN_ON(!entity_is_task(se));
	return container_of(se, struct task_struct, se);
}

/* Walk up scheduling entities hierarchy */
#define for_each_sched_entity(se) \
		for (; se; se = se->parent)

static inline struct cfs_rq *task_cfs_rq(struct task_struct *p)
{
	return p->se.cfs_rq;
}

/* runqueue on which this entity is (to be) queued */
static inline struct cfs_rq *cfs_rq_of(struct sched_entity *se)
{
	return se->cfs_rq;
}

/* runqueue "owned" by this group */
static inline struct cfs_rq *group_cfs_rq(struct sched_entity *grp)
{
	return grp->my_q;
}

static inline bool list_add_leaf_cfs_rq(struct cfs_rq *cfs_rq)
{
	struct rq *rq = rq_of(cfs_rq);
	int cpu = cpu_of(rq);

	if (cfs_rq->on_list)
		return rq->tmp_alone_branch == &rq->leaf_cfs_rq_list;

	cfs_rq->on_list = 1;

	/*
	 * Ensure we either appear before our parent (if already
	 * enqueued) or force our parent to appear after us when it is
	 * enqueued. The fact that we always enqueue bottom-up
	 * reduces this to two cases and a special case for the root
	 * cfs_rq. Furthermore, it also means that we will always reset
	 * tmp_alone_branch either when the branch is connected
	 * to a tree or when we reach the top of the tree
	 */
	if (cfs_rq->tg->parent &&
	    cfs_rq->tg->parent->cfs_rq[cpu]->on_list) {
		/*
		 * If parent is already on the list, we add the child
		 * just before. Thanks to circular linked property of
		 * the list, this means to put the child at the tail
		 * of the list that starts by parent.
		 */
		list_add_tail_rcu(&cfs_rq->leaf_cfs_rq_list,
			&(cfs_rq->tg->parent->cfs_rq[cpu]->leaf_cfs_rq_list));
		/*
		 * The branch is now connected to its tree so we can
		 * reset tmp_alone_branch to the beginning of the
		 * list.
		 */
		rq->tmp_alone_branch = &rq->leaf_cfs_rq_list;
		return true;
	}

	if (!cfs_rq->tg->parent) {
		/*
		 * cfs rq without parent should be put
		 * at the tail of the list.
		 */
		list_add_tail_rcu(&cfs_rq->leaf_cfs_rq_list,
			&rq->leaf_cfs_rq_list);
		/*
		 * We have reach the top of a tree so we can reset
		 * tmp_alone_branch to the beginning of the list.
		 */
		rq->tmp_alone_branch = &rq->leaf_cfs_rq_list;
		return true;
	}

	/*
	 * The parent has not already been added so we want to
	 * make sure that it will be put after us.
	 * tmp_alone_branch points to the begin of the branch
	 * where we will add parent.
	 */
	list_add_rcu(&cfs_rq->leaf_cfs_rq_list, rq->tmp_alone_branch);
	/*
	 * update tmp_alone_branch to points to the new begin
	 * of the branch
	 */
	rq->tmp_alone_branch = &cfs_rq->leaf_cfs_rq_list;
	return false;
}

static inline void list_del_leaf_cfs_rq(struct cfs_rq *cfs_rq)
{
	if (cfs_rq->on_list) {
		struct rq *rq = rq_of(cfs_rq);

		/*
		 * With cfs_rq being unthrottled/throttled during an enqueue,
		 * it can happen the tmp_alone_branch points the a leaf that
		 * we finally want to del. In this case, tmp_alone_branch moves
		 * to the prev element but it will point to rq->leaf_cfs_rq_list
		 * at the end of the enqueue.
		 */
		if (rq->tmp_alone_branch == &cfs_rq->leaf_cfs_rq_list)
			rq->tmp_alone_branch = cfs_rq->leaf_cfs_rq_list.prev;

		list_del_rcu(&cfs_rq->leaf_cfs_rq_list);
		cfs_rq->on_list = 0;
	}
}

static inline void assert_list_leaf_cfs_rq(struct rq *rq)
{
	SCHED_WARN_ON(rq->tmp_alone_branch != &rq->leaf_cfs_rq_list);
}

/* Iterate thr' all leaf cfs_rq's on a runqueue */
#define for_each_leaf_cfs_rq_safe(rq, cfs_rq, pos)			\
	list_for_each_entry_safe(cfs_rq, pos, &rq->leaf_cfs_rq_list,	\
				 leaf_cfs_rq_list)

/* Do the two (enqueued) entities belong to the same group ? */
static inline struct cfs_rq *
is_same_group(struct sched_entity *se, struct sched_entity *pse)
{
	if (se->cfs_rq == pse->cfs_rq)
		return se->cfs_rq;

	return NULL;
}

static inline struct sched_entity *parent_entity(struct sched_entity *se)
{
	return se->parent;
}

static void
find_matching_se(struct sched_entity **se, struct sched_entity **pse)
{
	int se_depth, pse_depth;

	/*
	 * preemption test can be made between sibling entities who are in the
	 * same cfs_rq i.e who have a common parent. Walk up the hierarchy of
	 * both tasks until we find their ancestors who are siblings of common
	 * parent.
	 */

	/* First walk up until both entities are at same depth */
	se_depth = (*se)->depth;
	pse_depth = (*pse)->depth;

	while (se_depth > pse_depth) {
		se_depth--;
		*se = parent_entity(*se);
	}

	while (pse_depth > se_depth) {
		pse_depth--;
		*pse = parent_entity(*pse);
	}

	while (!is_same_group(*se, *pse)) {
		*se = parent_entity(*se);
		*pse = parent_entity(*pse);
	}
}

#else	/* !CONFIG_FAIR_GROUP_SCHED */

static inline struct task_struct *task_of(struct sched_entity *se)
{
	return container_of(se, struct task_struct, se);
}

#define for_each_sched_entity(se) \
		for (; se; se = NULL)

static inline struct cfs_rq *task_cfs_rq(struct task_struct *p)
{
	return &task_rq(p)->cfs;
}

static inline struct cfs_rq *cfs_rq_of(struct sched_entity *se)
{
	struct task_struct *p = task_of(se);
	struct rq *rq = task_rq(p);

	return &rq->cfs;
}

/* runqueue "owned" by this group */
static inline struct cfs_rq *group_cfs_rq(struct sched_entity *grp)
{
	return NULL;
}

static inline bool list_add_leaf_cfs_rq(struct cfs_rq *cfs_rq)
{
	return true;
}

static inline void list_del_leaf_cfs_rq(struct cfs_rq *cfs_rq)
{
}

static inline void assert_list_leaf_cfs_rq(struct rq *rq)
{
}

#define for_each_leaf_cfs_rq_safe(rq, cfs_rq, pos)	\
		for (cfs_rq = &rq->cfs, pos = NULL; cfs_rq; cfs_rq = pos)

static inline struct sched_entity *parent_entity(struct sched_entity *se)
{
	return NULL;
}

static inline void
find_matching_se(struct sched_entity **se, struct sched_entity **pse)
{
}

#endif	/* CONFIG_FAIR_GROUP_SCHED */

static __always_inline
void account_cfs_rq_runtime(struct cfs_rq *cfs_rq, u64 delta_exec);

/**************************************************************
 * Scheduling class tree data structure manipulation methods:
 */

static inline u64 max_vruntime(u64 max_vruntime, u64 vruntime)
{
	s64 delta = (s64)(vruntime - max_vruntime);
	if (delta > 0)
		max_vruntime = vruntime;

	return max_vruntime;
}

static inline u64 min_vruntime(u64 min_vruntime, u64 vruntime)
{
	s64 delta = (s64)(vruntime - min_vruntime);
	if (delta < 0)
		min_vruntime = vruntime;

	return min_vruntime;
}

static inline bool entity_before(struct sched_entity *a,
				struct sched_entity *b)
{
	/*
	 * Tiebreak on vruntime seems unnecessary since it can
	 * hardly happen.
	 */
	return (s64)(a->deadline - b->deadline) < 0;
}

static inline s64 entity_key(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	return (s64)(se->vruntime - cfs_rq->min_vruntime);
}

#define __node_2_se(node) \
	rb_entry((node), struct sched_entity, run_node)

/*
 * Compute virtual time from the per-task service numbers:
 *
 * Fair schedulers conserve lag:
 *
 *   \Sum lag_i = 0
 *
 * Where lag_i is given by:
 *
 *   lag_i = S - s_i = w_i * (V - v_i)
 *
 * Where S is the ideal service time and V is it's virtual time counterpart.
 * Therefore:
 *
 *   \Sum lag_i = 0
 *   \Sum w_i * (V - v_i) = 0
 *   \Sum w_i * V - w_i * v_i = 0
 *
 * From which we can solve an expression for V in v_i (which we have in
 * se->vruntime):
 *
 *       \Sum v_i * w_i   \Sum v_i * w_i
 *   V = -------------- = --------------
 *          \Sum w_i            W
 *
 * Specifically, this is the weighted average of all entity virtual runtimes.
 *
 * [[ NOTE: this is only equal to the ideal scheduler under the condition
 *          that join/leave operations happen at lag_i = 0, otherwise the
 *          virtual time has non-continguous motion equivalent to:
 *
 *	      V +-= lag_i / W
 *
 *	    Also see the comment in place_entity() that deals with this. ]]
 *
 * However, since v_i is u64, and the multiplcation could easily overflow
 * transform it into a relative form that uses smaller quantities:
 *
 * Substitute: v_i == (v_i - v0) + v0
 *
 *     \Sum ((v_i - v0) + v0) * w_i   \Sum (v_i - v0) * w_i
 * V = ---------------------------- = --------------------- + v0
 *                  W                            W
 *
 * Which we track using:
 *
 *                    v0 := cfs_rq->min_vruntime
 * \Sum (v_i - v0) * w_i := cfs_rq->avg_vruntime
 *              \Sum w_i := cfs_rq->avg_load
 *
 * Since min_vruntime is a monotonic increasing variable that closely tracks
 * the per-task service, these deltas: (v_i - v), will be in the order of the
 * maximal (virtual) lag induced in the system due to quantisation.
 *
 * Also, we use scale_load_down() to reduce the size.
 *
 * As measured, the max (key * weight) value was ~44 bits for a kernel build.
 */
static void
avg_vruntime_add(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	unsigned long weight = scale_load_down(se->load.weight);
	s64 key = entity_key(cfs_rq, se);

	cfs_rq->avg_vruntime += key * weight;
	cfs_rq->avg_load += weight;
}

static void
avg_vruntime_sub(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	unsigned long weight = scale_load_down(se->load.weight);
	s64 key = entity_key(cfs_rq, se);

	cfs_rq->avg_vruntime -= key * weight;
	cfs_rq->avg_load -= weight;
}

static inline
void avg_vruntime_update(struct cfs_rq *cfs_rq, s64 delta)
{
	/*
	 * v' = v + d ==> avg_vruntime' = avg_runtime - d*avg_load
	 */
	cfs_rq->avg_vruntime -= cfs_rq->avg_load * delta;
}

/*
 * Specifically: avg_runtime() + 0 must result in entity_eligible() := true
 * For this to be so, the result of this function must have a left bias.
 */
u64 avg_vruntime(struct cfs_rq *cfs_rq)
{
	struct sched_entity *curr = cfs_rq->curr;
	s64 avg = cfs_rq->avg_vruntime;
	long load = cfs_rq->avg_load;

	if (curr && curr->on_rq) {
		unsigned long weight = scale_load_down(curr->load.weight);

		avg += entity_key(cfs_rq, curr) * weight;
		load += weight;
	}

	if (load) {
		/* sign flips effective floor / ceil */
		if (avg < 0)
			avg -= (load - 1);
		avg = div_s64(avg, load);
	}

	return cfs_rq->min_vruntime + avg;
}

/*
 * lag_i = S - s_i = w_i * (V - v_i)
 *
 * However, since V is approximated by the weighted average of all entities it
 * is possible -- by addition/removal/reweight to the tree -- to move V around
 * and end up with a larger lag than we started with.
 *
 * Limit this to either double the slice length with a minimum of TICK_NSEC
 * since that is the timing granularity.
 *
 * EEVDF gives the following limit for a steady state system:
 *
 *   -r_max < lag < max(r_max, q)
 *
 * XXX could add max_slice to the augmented data to track this.
 */
static s64 entity_lag(u64 avruntime, struct sched_entity *se)
{
	s64 vlag, limit;

	vlag = avruntime - se->vruntime;
	limit = calc_delta_fair(max_t(u64, 2*se->slice, TICK_NSEC), se);

	return clamp(vlag, -limit, limit);
}

static void update_entity_lag(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	SCHED_WARN_ON(!se->on_rq);

	se->vlag = entity_lag(avg_vruntime(cfs_rq), se);
}

/*
 * Entity is eligible once it received less service than it ought to have,
 * eg. lag >= 0.
 *
 * lag_i = S - s_i = w_i*(V - v_i)
 *
 * lag_i >= 0 -> V >= v_i
 *
 *     \Sum (v_i - v)*w_i
 * V = ------------------ + v
 *          \Sum w_i
 *
 * lag_i >= 0 -> \Sum (v_i - v)*w_i >= (v_i - v)*(\Sum w_i)
 *
 * Note: using 'avg_vruntime() > se->vruntime' is inacurate due
 *       to the loss in precision caused by the division.
 */
static int vruntime_eligible(struct cfs_rq *cfs_rq, u64 vruntime)
{
	struct sched_entity *curr = cfs_rq->curr;
	s64 avg = cfs_rq->avg_vruntime;
	long load = cfs_rq->avg_load;

	if (curr && curr->on_rq) {
		unsigned long weight = scale_load_down(curr->load.weight);

		avg += entity_key(cfs_rq, curr) * weight;
		load += weight;
	}

	return avg >= (s64)(vruntime - cfs_rq->min_vruntime) * load;
}

int entity_eligible(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	if (!sched_feat(ENFORCE_ELIGIBILITY))
		return 1;

	return vruntime_eligible(cfs_rq, se->vruntime);
}

static u64 __update_min_vruntime(struct cfs_rq *cfs_rq, u64 vruntime)
{
	u64 min_vruntime = cfs_rq->min_vruntime;
	/*
	 * open coded max_vruntime() to allow updating avg_vruntime
	 */
	s64 delta = (s64)(vruntime - min_vruntime);
	if (delta > 0) {
		avg_vruntime_update(cfs_rq, delta);
		min_vruntime = vruntime;
	}
	return min_vruntime;
}

static void update_min_vruntime(struct cfs_rq *cfs_rq)
{
	struct sched_entity *se = __pick_root_entity(cfs_rq);
	struct sched_entity *curr = cfs_rq->curr;
	u64 vruntime = cfs_rq->min_vruntime;

	if (curr) {
		if (curr->on_rq)
			vruntime = curr->vruntime;
		else
			curr = NULL;
	}

	if (se) {
		if (!curr)
			vruntime = se->min_vruntime;
		else
			vruntime = min_vruntime(vruntime, se->min_vruntime);
	}

	/* ensure we never gain time by being placed backwards. */
	u64_u32_store(cfs_rq->min_vruntime,
		      __update_min_vruntime(cfs_rq, vruntime));
}

static inline bool __entity_less(struct rb_node *a, const struct rb_node *b)
{
	return entity_before(__node_2_se(a), __node_2_se(b));
}

#define vruntime_gt(field, lse, rse) ({ (s64)((lse)->field - (rse)->field) > 0; })

static inline void __min_vruntime_update(struct sched_entity *se, struct rb_node *node)
{
	if (node) {
		struct sched_entity *rse = __node_2_se(node);
		if (vruntime_gt(min_vruntime, se, rse))
			se->min_vruntime = rse->min_vruntime;
	}
}

/*
 * se->min_vruntime = min(se->vruntime, {left,right}->min_vruntime)
 */
static inline bool min_vruntime_update(struct sched_entity *se, bool exit)
{
	u64 old_min_vruntime = se->min_vruntime;
	struct rb_node *node = &se->run_node;

	se->min_vruntime = se->vruntime;
	__min_vruntime_update(se, node->rb_right);
	__min_vruntime_update(se, node->rb_left);

	return se->min_vruntime == old_min_vruntime;
}

RB_DECLARE_CALLBACKS(static, min_vruntime_cb, struct sched_entity,
		     run_node, min_vruntime, min_vruntime_update);

/*
 * Enqueue an entity into the rb-tree:
 */
static void __enqueue_entity(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	avg_vruntime_add(cfs_rq, se);
	se->min_vruntime = se->vruntime;
	rb_add_augmented_cached(&se->run_node, &cfs_rq->tasks_timeline,
				__entity_less, &min_vruntime_cb);
}

static void __dequeue_entity(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	rb_erase_augmented_cached(&se->run_node, &cfs_rq->tasks_timeline,
				  &min_vruntime_cb);
	avg_vruntime_sub(cfs_rq, se);
}

struct sched_entity *__pick_root_entity(struct cfs_rq *cfs_rq)
{
	struct rb_node *root = cfs_rq->tasks_timeline.rb_root.rb_node;

	if (!root)
		return NULL;

	return __node_2_se(root);
}

struct sched_entity *__pick_first_entity(struct cfs_rq *cfs_rq)
{
	struct rb_node *left = rb_first_cached(&cfs_rq->tasks_timeline);

	if (!left)
		return NULL;

	return __node_2_se(left);
}

/*
 * Earliest Eligible Virtual Deadline First
 *
 * In order to provide latency guarantees for different request sizes
 * EEVDF selects the best runnable task from two criteria:
 *
 *  1) the task must be eligible (must be owed service)
 *
 *  2) from those tasks that meet 1), we select the one
 *     with the earliest virtual deadline.
 *
 * We can do this in O(log n) time due to an augmented RB-tree. The
 * tree keeps the entries sorted on deadline, but also functions as a
 * heap based on the vruntime by keeping:
 *
 *  se->min_vruntime = min(se->vruntime, se->{left,right}->min_vruntime)
 *
 * Which allows tree pruning through eligibility.
 */
static struct sched_entity *pick_eevdf(struct cfs_rq *cfs_rq)
{
	struct rb_node *node = cfs_rq->tasks_timeline.rb_root.rb_node;
	struct sched_entity *se = __pick_first_entity(cfs_rq);
	struct sched_entity *curr = cfs_rq->curr;
	struct sched_entity *best = NULL;

	/*
	 * We can safely skip eligibility check if there is only one entity
	 * in this cfs_rq, saving some cycles.
	 */
	if (cfs_rq->nr_running == 1)
		return curr && curr->on_rq ? curr : se;

	if (curr && (!curr->on_rq || !entity_eligible(cfs_rq, curr)))
		curr = NULL;

	/*
	 * Once selected, run a task until it either becomes non-eligible or
	 * until it gets a new slice. See the HACK in set_next_entity().
	 */
	if (sched_feat(RUN_TO_PARITY) && curr && curr->vlag == curr->deadline)
		return curr;

	/* Pick the leftmost entity if it's eligible */
	if (se && entity_eligible(cfs_rq, se)) {
		best = se;
		goto found;
	}

	/* Heap search for the EEVD entity */
	while (node) {
		struct rb_node *left = node->rb_left;

		/*
		 * Eligible entities in left subtree are always better
		 * choices, since they have earlier deadlines.
		 */
		if (left && vruntime_eligible(cfs_rq,
					__node_2_se(left)->min_vruntime)) {
			node = left;
			continue;
		}

		se = __node_2_se(node);

		/*
		 * The left subtree either is empty or has no eligible
		 * entity, so check the current node since it is the one
		 * with earliest deadline that might be eligible.
		 */
		if (entity_eligible(cfs_rq, se)) {
			best = se;
			break;
		}

		node = node->rb_right;
	}
found:
	if (!best || (curr && entity_before(curr, best)))
		best = curr;

	return best;
}

#ifdef CONFIG_SCHED_DEBUG
struct sched_entity *__pick_last_entity(struct cfs_rq *cfs_rq)
{
	struct rb_node *last = rb_last(&cfs_rq->tasks_timeline.rb_root);

	if (!last)
		return NULL;

	return __node_2_se(last);
}

/**************************************************************
 * Scheduling class statistics methods:
 */

int sched_proc_update_handler(struct ctl_table *table, int write,
		void __user *buffer, size_t *lenp,
		loff_t *ppos)
{
	int ret = proc_dointvec_minmax(table, write, buffer, lenp, ppos);
	unsigned int factor = get_update_sysctl_factor();

	if (ret || !write)
		return ret;

#define WRT_SYSCTL(name) \
	(normalized_sysctl_##name = sysctl_##name / (factor))
	WRT_SYSCTL(sched_base_slice);
#undef WRT_SYSCTL

	return 0;
}
#endif

static void clear_buddies(struct cfs_rq *cfs_rq, struct sched_entity *se);

/*
 * XXX: strictly: vd_i += N*r_i/w_i such that: vd_i > ve_i
 * this is probably good enough.
 */
static void update_deadline(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	if ((s64)(se->vruntime - se->deadline) < 0)
		return;

	/*
	 * For EEVDF the virtual time slope is determined by w_i (iow.
	 * nice) while the request time r_i is determined by
	 * sysctl_sched_base_slice.
	 */
	se->slice = sysctl_sched_base_slice;

	/*
	 * EEVDF: vd_i = ve_i + r_i / w_i
	 */
	se->deadline = se->vruntime + calc_delta_fair(se->slice, se);

	/*
	 * The task has consumed its request, reschedule.
	 */
	if (cfs_rq->nr_running > 1) {
		resched_curr(rq_of(cfs_rq));
		clear_buddies(cfs_rq, se);
	}
}

#include "pelt.h"
#ifdef CONFIG_SMP

static int select_idle_sibling(struct task_struct *p, int prev_cpu, int cpu);
static unsigned long task_h_load(struct task_struct *p);
static unsigned long capacity_of(int cpu);

/* Give new sched_entity start runnable values to heavy its load in infant time */
void init_entity_runnable_average(struct sched_entity *se)
{
	struct sched_avg *sa = &se->avg;

	memset(sa, 0, sizeof(*sa));

	/*
	 * Tasks are intialized with full load to be seen as heavy tasks until
	 * they get a chance to stabilize to their real load level.
	 * Group entities are intialized with zero load to reflect the fact that
	 * nothing has been attached to the task group yet.
	 */
	if (entity_is_task(se))
		sa->runnable_load_avg = sa->load_avg = scale_load_down(se->load.weight);

	se->runnable_weight = se->load.weight;

	/* when this task enqueue'ed, it will contribute to its cfs_rq's load_avg */
}

static inline u64 cfs_rq_clock_task(struct cfs_rq *cfs_rq);
static void attach_entity_cfs_rq(struct sched_entity *se);

/*
 * With new tasks being created, their initial util_avgs are extrapolated
 * based on the cfs_rq's current util_avg:
 *
 *   util_avg = cfs_rq->avg.util_avg / (cfs_rq->avg.load_avg + 1)
 *		* se_weight(se)
 *
 * However, in many cases, the above util_avg does not give a desired
 * value. Moreover, the sum of the util_avgs may be divergent, such
 * as when the series is a harmonic series.
 *
 * To solve this problem, we also cap the util_avg of successive tasks to
 * only 1/2 of the left utilization budget:
 *
 *   util_avg_cap = (cpu_scale - cfs_rq->avg.util_avg) / 2^n
 *
 * where n denotes the nth task and cpu_scale the CPU capacity.
 *
 * For example, for a CPU with 1024 of capacity, a simplest series from
 * the beginning would be like:
 *
 *  task  util_avg: 512, 256, 128,  64,  32,   16,    8, ...
 * cfs_rq util_avg: 512, 768, 896, 960, 992, 1008, 1016, ...
 *
 * Finally, that extrapolated util_avg is clamped to the cap (util_avg_cap)
 * if util_avg > util_avg_cap.
 */
void post_init_entity_util_avg(struct sched_entity *se)
{
	struct cfs_rq *cfs_rq = cfs_rq_of(se);
	struct sched_avg *sa = &se->avg;
	long cpu_scale = arch_scale_cpu_capacity(NULL, cpu_of(rq_of(cfs_rq)));
	long cap = (long)(cpu_scale - cfs_rq->avg.util_avg) / 2;

	if (cap > 0) {
		if (cfs_rq->avg.util_avg != 0) {
			sa->util_avg  = cfs_rq->avg.util_avg * se_weight(se);
			sa->util_avg /= (cfs_rq->avg.load_avg + 1);

			if (sa->util_avg > cap)
				sa->util_avg = cap;
		} else {
			sa->util_avg = cap;
		}
	}

	if (entity_is_task(se)) {
		struct task_struct *p = task_of(se);
		if (p->sched_class != &fair_sched_class) {
			/*
			 * For !fair tasks do:
			 *
			update_cfs_rq_load_avg(now, cfs_rq);
			attach_entity_load_avg(cfs_rq, se, 0);
			switched_from_fair(rq, p);
			 *
			 * such that the next switched_to_fair() has the
			 * expected state.
			 */
			se->avg.last_update_time = cfs_rq_clock_pelt(cfs_rq);
			return;
		}
	}

	attach_entity_cfs_rq(se);
}

#else /* !CONFIG_SMP */
void init_entity_runnable_average(struct sched_entity *se)
{
}
void post_init_entity_util_avg(struct sched_entity *se)
{
}
static void update_tg_load_avg(struct cfs_rq *cfs_rq, int force)
{
}
#endif /* CONFIG_SMP */

/*
 * Update the current task's runtime statistics.
 */
static void update_curr(struct cfs_rq *cfs_rq)
{
	struct sched_entity *curr = cfs_rq->curr;
	u64 now = rq_clock_task(rq_of(cfs_rq));
	u64 delta_exec;

	if (unlikely(!curr))
		return;

	delta_exec = now - curr->exec_start;
	if (unlikely((s64)delta_exec <= 0))
		return;

	curr->exec_start = now;

	schedstat_set(curr->statistics.exec_max,
		      max(delta_exec, curr->statistics.exec_max));

	curr->sum_exec_runtime += delta_exec;
	schedstat_add(cfs_rq->exec_clock, delta_exec);

	curr->vruntime += calc_delta_fair(delta_exec, curr);
	update_deadline(cfs_rq, curr);
	update_min_vruntime(cfs_rq);

	if (entity_is_task(curr)) {
		struct task_struct *curtask = task_of(curr);

		trace_sched_stat_runtime(curtask, delta_exec, curr->vruntime);
		cgroup_account_cputime(curtask, delta_exec);
		account_group_exec_runtime(curtask, delta_exec);
	}

	account_cfs_rq_runtime(cfs_rq, delta_exec);
}

static void update_curr_fair(struct rq *rq)
{
	update_curr(cfs_rq_of(&rq->curr->se));
}

static inline void
update_stats_wait_start(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	u64 wait_start, prev_wait_start;

	if (!schedstat_enabled())
		return;

	wait_start = rq_clock(rq_of(cfs_rq));
	prev_wait_start = schedstat_val(se->statistics.wait_start);

	if (entity_is_task(se) && task_on_rq_migrating(task_of(se)) &&
	    likely(wait_start > prev_wait_start))
		wait_start -= prev_wait_start;

	__schedstat_set(se->statistics.wait_start, wait_start);
}

static inline void
update_stats_wait_end(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	struct task_struct *p;
	u64 delta;

	if (!schedstat_enabled())
		return;

	delta = rq_clock(rq_of(cfs_rq)) - schedstat_val(se->statistics.wait_start);

	if (entity_is_task(se)) {
		p = task_of(se);
		if (task_on_rq_migrating(p)) {
			/*
			 * Preserve migrating task's wait time so wait_start
			 * time stamp can be adjusted to accumulate wait time
			 * prior to migration.
			 */
			__schedstat_set(se->statistics.wait_start, delta);
			return;
		}
		trace_sched_stat_wait(p, delta);
	}

	__schedstat_set(se->statistics.wait_max,
		      max(schedstat_val(se->statistics.wait_max), delta));
	__schedstat_inc(se->statistics.wait_count);
	__schedstat_add(se->statistics.wait_sum, delta);
	__schedstat_set(se->statistics.wait_start, 0);
}

static inline void
update_stats_enqueue_sleeper(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	struct task_struct *tsk = NULL;
	u64 sleep_start, block_start;

	if (!schedstat_enabled())
		return;

	sleep_start = schedstat_val(se->statistics.sleep_start);
	block_start = schedstat_val(se->statistics.block_start);

	if (entity_is_task(se))
		tsk = task_of(se);

	if (sleep_start) {
		u64 delta = rq_clock(rq_of(cfs_rq)) - sleep_start;

		if ((s64)delta < 0)
			delta = 0;

		if (unlikely(delta > schedstat_val(se->statistics.sleep_max)))
			__schedstat_set(se->statistics.sleep_max, delta);

		__schedstat_set(se->statistics.sleep_start, 0);
		__schedstat_add(se->statistics.sum_sleep_runtime, delta);

		if (tsk) {
			account_scheduler_latency(tsk, delta >> 10, 1);
			trace_sched_stat_sleep(tsk, delta);
		}
	}
	if (block_start) {
		u64 delta = rq_clock(rq_of(cfs_rq)) - block_start;

		if ((s64)delta < 0)
			delta = 0;

		if (unlikely(delta > schedstat_val(se->statistics.block_max)))
			__schedstat_set(se->statistics.block_max, delta);

		__schedstat_set(se->statistics.block_start, 0);
		__schedstat_add(se->statistics.sum_sleep_runtime, delta);

		if (tsk) {
			if (tsk->in_iowait) {
				__schedstat_add(se->statistics.iowait_sum, delta);
				__schedstat_inc(se->statistics.iowait_count);
				trace_sched_stat_iowait(tsk, delta);
			}

			trace_sched_stat_blocked(tsk, delta);

			/*
			 * Blocking time is in units of nanosecs, so shift by
			 * 20 to get a milliseconds-range estimation of the
			 * amount of time that the task spent sleeping:
			 */
			if (unlikely(prof_on == SLEEP_PROFILING)) {
				profile_hits(SLEEP_PROFILING,
						(void *)get_wchan(tsk),
						delta >> 20);
			}
			account_scheduler_latency(tsk, delta >> 10, 0);
		}
	}
}

/*
 * Task is being enqueued - update stats:
 */
static inline void
update_stats_enqueue(struct cfs_rq *cfs_rq, struct sched_entity *se, int flags)
{
	if (!schedstat_enabled())
		return;

	/*
	 * Are we enqueueing a waiting task? (for current tasks
	 * a dequeue/enqueue event is a NOP)
	 */
	if (se != cfs_rq->curr)
		update_stats_wait_start(cfs_rq, se);

	if (flags & ENQUEUE_WAKEUP)
		update_stats_enqueue_sleeper(cfs_rq, se);
}

static inline void
update_stats_dequeue(struct cfs_rq *cfs_rq, struct sched_entity *se, int flags)
{

	if (!schedstat_enabled())
		return;

	/*
	 * Mark the end of the wait period if dequeueing a
	 * waiting task:
	 */
	if (se != cfs_rq->curr)
		update_stats_wait_end(cfs_rq, se);

	if ((flags & DEQUEUE_SLEEP) && entity_is_task(se)) {
		struct task_struct *tsk = task_of(se);

		if (tsk->state & TASK_INTERRUPTIBLE)
			__schedstat_set(se->statistics.sleep_start,
				      rq_clock(rq_of(cfs_rq)));
		if (tsk->state & TASK_UNINTERRUPTIBLE)
			__schedstat_set(se->statistics.block_start,
				      rq_clock(rq_of(cfs_rq)));
	}
}

/*
 * We are picking a new current task - update its stats:
 */
static inline void
update_stats_curr_start(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	/*
	 * We are starting a new run period:
	 */
	se->exec_start = rq_clock_task(rq_of(cfs_rq));
}

/**************************************************
 * Scheduling class queueing methods:
 */

#ifdef CONFIG_NUMA_BALANCING
/*
 * Approximate time to scan a full NUMA task in ms. The task scan period is
 * calculated based on the tasks virtual memory size and
 * numa_balancing_scan_size.
 */
unsigned int sysctl_numa_balancing_scan_period_min = 1000;
unsigned int sysctl_numa_balancing_scan_period_max = 60000;

/* Portion of address space to scan in MB */
unsigned int sysctl_numa_balancing_scan_size = 256;

/* Scan @scan_size MB every @scan_period after an initial @scan_delay in ms */
unsigned int sysctl_numa_balancing_scan_delay = 1000;

struct numa_group {
	atomic_t refcount;

	spinlock_t lock; /* nr_tasks, tasks */
	int nr_tasks;
	pid_t gid;
	int active_nodes;

	struct rcu_head rcu;
	unsigned long total_faults;
	unsigned long max_faults_cpu;
	/*
	 * Faults_cpu is used to decide whether memory should move
	 * towards the CPU. As a consequence, these stats are weighted
	 * more by CPU use than by memory faults.
	 */
	unsigned long *faults_cpu;
	unsigned long faults[0];
};

/*
 * For functions that can be called in multiple contexts that permit reading
 * ->numa_group (see struct task_struct for locking rules).
 */
static struct numa_group *deref_task_numa_group(struct task_struct *p)
{
	return rcu_dereference_check(p->numa_group, p == current ||
		(lockdep_is_held(&task_rq(p)->lock) && !READ_ONCE(p->on_cpu)));
}

static struct numa_group *deref_curr_numa_group(struct task_struct *p)
{
	return rcu_dereference_protected(p->numa_group, p == current);
}

static inline unsigned long group_faults_priv(struct numa_group *ng);
static inline unsigned long group_faults_shared(struct numa_group *ng);

static unsigned int task_nr_scan_windows(struct task_struct *p)
{
	unsigned long rss = 0;
	unsigned long nr_scan_pages;

	/*
	 * Calculations based on RSS as non-present and empty pages are skipped
	 * by the PTE scanner and NUMA hinting faults should be trapped based
	 * on resident pages
	 */
	nr_scan_pages = sysctl_numa_balancing_scan_size << (20 - PAGE_SHIFT);
	rss = get_mm_rss(p->mm);
	if (!rss)
		rss = nr_scan_pages;

	rss = round_up(rss, nr_scan_pages);
	return rss / nr_scan_pages;
}

/* For sanitys sake, never scan more PTEs than MAX_SCAN_WINDOW MB/sec. */
#define MAX_SCAN_WINDOW 2560

static unsigned int task_scan_min(struct task_struct *p)
{
	unsigned int scan_size = READ_ONCE(sysctl_numa_balancing_scan_size);
	unsigned int scan, floor;
	unsigned int windows = 1;

	if (scan_size < MAX_SCAN_WINDOW)
		windows = MAX_SCAN_WINDOW / scan_size;
	floor = 1000 / windows;

	scan = sysctl_numa_balancing_scan_period_min / task_nr_scan_windows(p);
	return max_t(unsigned int, floor, scan);
}

static unsigned int task_scan_start(struct task_struct *p)
{
	unsigned long smin = task_scan_min(p);
	unsigned long period = smin;
	struct numa_group *ng;

	/* Scale the maximum scan period with the amount of shared memory. */
	rcu_read_lock();
	ng = rcu_dereference(p->numa_group);
	if (ng) {
		unsigned long shared = group_faults_shared(ng);
		unsigned long private = group_faults_priv(ng);

		period *= atomic_read(&ng->refcount);
		period *= shared + 1;
		period /= private + shared + 1;
	}
	rcu_read_unlock();

	return max(smin, period);
}

static unsigned int task_scan_max(struct task_struct *p)
{
	unsigned long smin = task_scan_min(p);
	unsigned long smax;
	struct numa_group *ng;

	/* Watch for min being lower than max due to floor calculations */
	smax = sysctl_numa_balancing_scan_period_max / task_nr_scan_windows(p);

	/* Scale the maximum scan period with the amount of shared memory. */
	ng = deref_curr_numa_group(p);
	if (ng) {
		unsigned long shared = group_faults_shared(ng);
		unsigned long private = group_faults_priv(ng);
		unsigned long period = smax;

		period *= atomic_read(&ng->refcount);
		period *= shared + 1;
		period /= private + shared + 1;

		smax = max(smax, period);
	}

	return max(smin, smax);
}

void init_numa_balancing(unsigned long clone_flags, struct task_struct *p)
{
	int mm_users = 0;
	struct mm_struct *mm = p->mm;

	if (mm) {
		mm_users = atomic_read(&mm->mm_users);
		if (mm_users == 1) {
			mm->numa_next_scan = jiffies + msecs_to_jiffies(sysctl_numa_balancing_scan_delay);
			mm->numa_scan_seq = 0;
		}
	}
	p->node_stamp			= 0;
	p->numa_scan_seq		= mm ? mm->numa_scan_seq : 0;
	p->numa_scan_period		= sysctl_numa_balancing_scan_delay;
	p->numa_work.next		= &p->numa_work;
	p->numa_faults			= NULL;
	RCU_INIT_POINTER(p->numa_group, NULL);
	p->last_task_numa_placement	= 0;
	p->last_sum_exec_runtime	= 0;

	/* New address space, reset the preferred nid */
	if (!(clone_flags & CLONE_VM)) {
		p->numa_preferred_nid = -1;
		return;
	}

	/*
	 * New thread, keep existing numa_preferred_nid which should be copied
	 * already by arch_dup_task_struct but stagger when scans start.
	 */
	if (mm) {
		unsigned int delay;

		delay = min_t(unsigned int, task_scan_max(current),
			current->numa_scan_period * mm_users * NSEC_PER_MSEC);
		delay += 2 * TICK_NSEC;
		p->node_stamp = delay;
	}
}

static void account_numa_enqueue(struct rq *rq, struct task_struct *p)
{
	rq->nr_numa_running += (p->numa_preferred_nid != -1);
	rq->nr_preferred_running += (p->numa_preferred_nid == task_node(p));
}

static void account_numa_dequeue(struct rq *rq, struct task_struct *p)
{
	rq->nr_numa_running -= (p->numa_preferred_nid != -1);
	rq->nr_preferred_running -= (p->numa_preferred_nid == task_node(p));
}

/* Shared or private faults. */
#define NR_NUMA_HINT_FAULT_TYPES 2

/* Memory and CPU locality */
#define NR_NUMA_HINT_FAULT_STATS (NR_NUMA_HINT_FAULT_TYPES * 2)

/* Averaged statistics, and temporary buffers. */
#define NR_NUMA_HINT_FAULT_BUCKETS (NR_NUMA_HINT_FAULT_STATS * 2)

pid_t task_numa_group_id(struct task_struct *p)
{
	struct numa_group *ng;
	pid_t gid = 0;

	rcu_read_lock();
	ng = rcu_dereference(p->numa_group);
	if (ng)
		gid = ng->gid;
	rcu_read_unlock();

	return gid;
}

/*
 * The averaged statistics, shared & private, memory & CPU,
 * occupy the first half of the array. The second half of the
 * array is for current counters, which are averaged into the
 * first set by task_numa_placement.
 */
static inline int task_faults_idx(enum numa_faults_stats s, int nid, int priv)
{
	return NR_NUMA_HINT_FAULT_TYPES * (s * nr_node_ids + nid) + priv;
}

static inline unsigned long task_faults(struct task_struct *p, int nid)
{
	if (!p->numa_faults)
		return 0;

	return p->numa_faults[task_faults_idx(NUMA_MEM, nid, 0)] +
		p->numa_faults[task_faults_idx(NUMA_MEM, nid, 1)];
}

static inline unsigned long group_faults(struct task_struct *p, int nid)
{
	struct numa_group *ng = deref_task_numa_group(p);

	if (!ng)
		return 0;

	return ng->faults[task_faults_idx(NUMA_MEM, nid, 0)] +
		ng->faults[task_faults_idx(NUMA_MEM, nid, 1)];
}

static inline unsigned long group_faults_cpu(struct numa_group *group, int nid)
{
	return group->faults_cpu[task_faults_idx(NUMA_MEM, nid, 0)] +
		group->faults_cpu[task_faults_idx(NUMA_MEM, nid, 1)];
}

static inline unsigned long group_faults_priv(struct numa_group *ng)
{
	unsigned long faults = 0;
	int node;

	for_each_online_node(node) {
		faults += ng->faults[task_faults_idx(NUMA_MEM, node, 1)];
	}

	return faults;
}

static inline unsigned long group_faults_shared(struct numa_group *ng)
{
	unsigned long faults = 0;
	int node;

	for_each_online_node(node) {
		faults += ng->faults[task_faults_idx(NUMA_MEM, node, 0)];
	}

	return faults;
}

/*
 * A node triggering more than 1/3 as many NUMA faults as the maximum is
 * considered part of a numa group's pseudo-interleaving set. Migrations
 * between these nodes are slowed down, to allow things to settle down.
 */
#define ACTIVE_NODE_FRACTION 3

static bool numa_is_active_node(int nid, struct numa_group *ng)
{
	return group_faults_cpu(ng, nid) * ACTIVE_NODE_FRACTION > ng->max_faults_cpu;
}

/* Handle placement on systems where not all nodes are directly connected. */
static unsigned long score_nearby_nodes(struct task_struct *p, int nid,
					int maxdist, bool task)
{
	unsigned long score = 0;
	int node;

	/*
	 * All nodes are directly connected, and the same distance
	 * from each other. No need for fancy placement algorithms.
	 */
	if (sched_numa_topology_type == NUMA_DIRECT)
		return 0;

	/*
	 * This code is called for each node, introducing N^2 complexity,
	 * which should be ok given the number of nodes rarely exceeds 8.
	 */
	for_each_online_node(node) {
		unsigned long faults;
		int dist = node_distance(nid, node);

		/*
		 * The furthest away nodes in the system are not interesting
		 * for placement; nid was already counted.
		 */
		if (dist == sched_max_numa_distance || node == nid)
			continue;

		/*
		 * On systems with a backplane NUMA topology, compare groups
		 * of nodes, and move tasks towards the group with the most
		 * memory accesses. When comparing two nodes at distance
		 * "hoplimit", only nodes closer by than "hoplimit" are part
		 * of each group. Skip other nodes.
		 */
		if (sched_numa_topology_type == NUMA_BACKPLANE &&
					dist >= maxdist)
			continue;

		/* Add up the faults from nearby nodes. */
		if (task)
			faults = task_faults(p, node);
		else
			faults = group_faults(p, node);

		/*
		 * On systems with a glueless mesh NUMA topology, there are
		 * no fixed "groups of nodes". Instead, nodes that are not
		 * directly connected bounce traffic through intermediate
		 * nodes; a numa_group can occupy any set of nodes.
		 * The further away a node is, the less the faults count.
		 * This seems to result in good task placement.
		 */
		if (sched_numa_topology_type == NUMA_GLUELESS_MESH) {
			faults *= (sched_max_numa_distance - dist);
			faults /= (sched_max_numa_distance - LOCAL_DISTANCE);
		}

		score += faults;
	}

	return score;
}

/*
 * These return the fraction of accesses done by a particular task, or
 * task group, on a particular numa node.  The group weight is given a
 * larger multiplier, in order to group tasks together that are almost
 * evenly spread out between numa nodes.
 */
static inline unsigned long task_weight(struct task_struct *p, int nid,
					int dist)
{
	unsigned long faults, total_faults;

	if (!p->numa_faults)
		return 0;

	total_faults = p->total_numa_faults;

	if (!total_faults)
		return 0;

	faults = task_faults(p, nid);
	faults += score_nearby_nodes(p, nid, dist, true);

	return 1000 * faults / total_faults;
}

static inline unsigned long group_weight(struct task_struct *p, int nid,
					 int dist)
{
	struct numa_group *ng = deref_task_numa_group(p);
	unsigned long faults, total_faults;

	if (!ng)
		return 0;

	total_faults = ng->total_faults;

	if (!total_faults)
		return 0;

	faults = group_faults(p, nid);
	faults += score_nearby_nodes(p, nid, dist, false);

	return 1000 * faults / total_faults;
}

bool should_numa_migrate_memory(struct task_struct *p, struct page * page,
				int src_nid, int dst_cpu)
{
	struct numa_group *ng = deref_curr_numa_group(p);
	int dst_nid = cpu_to_node(dst_cpu);
	int last_cpupid, this_cpupid;

	this_cpupid = cpu_pid_to_cpupid(dst_cpu, current->pid);
	last_cpupid = page_cpupid_xchg_last(page, this_cpupid);

	/*
	 * Allow first faults or private faults to migrate immediately early in
	 * the lifetime of a task. The magic number 4 is based on waiting for
	 * two full passes of the "multi-stage node selection" test that is
	 * executed below.
	 */
	if ((p->numa_preferred_nid == -1 || p->numa_scan_seq <= 4) &&
	    (cpupid_pid_unset(last_cpupid) || cpupid_match_pid(p, last_cpupid)))
		return true;

	/*
	 * Multi-stage node selection is used in conjunction with a periodic
	 * migration fault to build a temporal task<->page relation. By using
	 * a two-stage filter we remove short/unlikely relations.
	 *
	 * Using P(p) ~ n_p / n_t as per frequentist probability, we can equate
	 * a task's usage of a particular page (n_p) per total usage of this
	 * page (n_t) (in a given time-span) to a probability.
	 *
	 * Our periodic faults will sample this probability and getting the
	 * same result twice in a row, given these samples are fully
	 * independent, is then given by P(n)^2, provided our sample period
	 * is sufficiently short compared to the usage pattern.
	 *
	 * This quadric squishes small probabilities, making it less likely we
	 * act on an unlikely task<->page relation.
	 */
	if (!cpupid_pid_unset(last_cpupid) &&
				cpupid_to_nid(last_cpupid) != dst_nid)
		return false;

	/* Always allow migrate on private faults */
	if (cpupid_match_pid(p, last_cpupid))
		return true;

	/* A shared fault, but p->numa_group has not been set up yet. */
	if (!ng)
		return true;

	/*
	 * Destination node is much more heavily used than the source
	 * node? Allow migration.
	 */
	if (group_faults_cpu(ng, dst_nid) > group_faults_cpu(ng, src_nid) *
					ACTIVE_NODE_FRACTION)
		return true;

	/*
	 * Distribute memory according to CPU & memory use on each node,
	 * with 3/4 hysteresis to avoid unnecessary memory migrations:
	 *
	 * faults_cpu(dst)   3   faults_cpu(src)
	 * --------------- * - > ---------------
	 * faults_mem(dst)   4   faults_mem(src)
	 */
	return group_faults_cpu(ng, dst_nid) * group_faults(p, src_nid) * 3 >
	       group_faults_cpu(ng, src_nid) * group_faults(p, dst_nid) * 4;
}

static unsigned long weighted_cpuload(struct rq *rq);
static unsigned long source_load(int cpu, int type);
static unsigned long target_load(int cpu, int type);

/* Cached statistics for all CPUs within a node */
struct numa_stats {
	unsigned long load;

	/* Total compute capacity of CPUs on a node */
	unsigned long compute_capacity;

	unsigned int nr_running;
};

/*
 * XXX borrowed from update_sg_lb_stats
 */
static void update_numa_stats(struct numa_stats *ns, int nid)
{
	int smt, cpu, cpus = 0;
	unsigned long capacity;

	memset(ns, 0, sizeof(*ns));
	for_each_cpu(cpu, cpumask_of_node(nid)) {
		struct rq *rq = cpu_rq(cpu);

		ns->nr_running += rq->nr_running;
		ns->load += weighted_cpuload(rq);
		ns->compute_capacity += capacity_of(cpu);

		cpus++;
	}

	/*
	 * If we raced with hotplug and there are no CPUs left in our mask
	 * the @ns structure is NULL'ed and task_numa_compare() will
	 * not find this node attractive.
	 *
	 * We'll detect a huge imbalance and bail there.
	 */
	if (!cpus)
		return;

	/* smt := ceil(cpus / capacity), assumes: 1 < smt_power < 2 */
	smt = DIV_ROUND_UP(SCHED_CAPACITY_SCALE * cpus, ns->compute_capacity);
	capacity = cpus / smt; /* cores */

	capacity = min_t(unsigned, capacity,
		DIV_ROUND_CLOSEST(ns->compute_capacity, SCHED_CAPACITY_SCALE));
}

struct task_numa_env {
	struct task_struct *p;

	int src_cpu, src_nid;
	int dst_cpu, dst_nid;

	struct numa_stats src_stats, dst_stats;

	int imbalance_pct;
	int dist;

	struct task_struct *best_task;
	long best_imp;
	int best_cpu;
};

static void task_numa_assign(struct task_numa_env *env,
			     struct task_struct *p, long imp)
{
	struct rq *rq = cpu_rq(env->dst_cpu);

	/* Bail out if run-queue part of active NUMA balance. */
	if (xchg(&rq->numa_migrate_on, 1))
		return;

	/*
	 * Clear previous best_cpu/rq numa-migrate flag, since task now
	 * found a better CPU to move/swap.
	 */
	if (env->best_cpu != -1) {
		rq = cpu_rq(env->best_cpu);
		WRITE_ONCE(rq->numa_migrate_on, 0);
	}

	if (env->best_task)
		put_task_struct(env->best_task);
	if (p)
		get_task_struct(p);

	env->best_task = p;
	env->best_imp = imp;
	env->best_cpu = env->dst_cpu;
}

static bool load_too_imbalanced(long src_load, long dst_load,
				struct task_numa_env *env)
{
	long imb, old_imb;
	long orig_src_load, orig_dst_load;
	long src_capacity, dst_capacity;

	/*
	 * The load is corrected for the CPU capacity available on each node.
	 *
	 * src_load        dst_load
	 * ------------ vs ---------
	 * src_capacity    dst_capacity
	 */
	src_capacity = env->src_stats.compute_capacity;
	dst_capacity = env->dst_stats.compute_capacity;

	imb = abs(dst_load * src_capacity - src_load * dst_capacity);

	orig_src_load = env->src_stats.load;
	orig_dst_load = env->dst_stats.load;

	old_imb = abs(orig_dst_load * src_capacity - orig_src_load * dst_capacity);

	/* Would this change make things worse? */
	return (imb > old_imb);
}

/*
 * Maximum NUMA importance can be 1998 (2*999);
 * SMALLIMP @ 30 would be close to 1998/64.
 * Used to deter task migration.
 */
#define SMALLIMP	30

/*
 * This checks if the overall compute and NUMA accesses of the system would
 * be improved if the source tasks was migrated to the target dst_cpu taking
 * into account that it might be best if task running on the dst_cpu should
 * be exchanged with the source task
 */
static void task_numa_compare(struct task_numa_env *env,
			      long taskimp, long groupimp, bool maymove)
{
	struct numa_group *cur_ng, *p_ng = deref_curr_numa_group(env->p);
	struct rq *dst_rq = cpu_rq(env->dst_cpu);
	long imp = p_ng ? groupimp : taskimp;
	struct task_struct *cur;
	long src_load, dst_load;
	int dist = env->dist;
	long moveimp = imp;
	long load;

	if (READ_ONCE(dst_rq->numa_migrate_on))
		return;

	rcu_read_lock();
	cur = task_rcu_dereference(&dst_rq->curr);
	if (cur && ((cur->flags & PF_EXITING) || is_idle_task(cur)))
		cur = NULL;

	/*
	 * Because we have preemption enabled we can get migrated around and
	 * end try selecting ourselves (current == env->p) as a swap candidate.
	 */
	if (cur == env->p)
		goto unlock;

	if (!cur) {
		if (maymove && moveimp >= env->best_imp)
			goto assign;
		else
			goto unlock;
	}

	/*
	 * "imp" is the fault differential for the source task between the
	 * source and destination node. Calculate the total differential for
	 * the source task and potential destination task. The more negative
	 * the value is, the more remote accesses that would be expected to
	 * be incurred if the tasks were swapped.
	 */
	/* Skip this swap candidate if cannot move to the source cpu */
	if (!cpumask_test_cpu(env->src_cpu, &cur->cpus_allowed))
		goto unlock;

	/*
	 * If dst and source tasks are in the same NUMA group, or not
	 * in any group then look only at task weights.
	 */
	cur_ng = rcu_dereference(cur->numa_group);
	if (cur_ng == p_ng) {
		imp = taskimp + task_weight(cur, env->src_nid, dist) -
		      task_weight(cur, env->dst_nid, dist);
		/*
		 * Add some hysteresis to prevent swapping the
		 * tasks within a group over tiny differences.
		 */
		if (cur_ng)
			imp -= imp / 16;
	} else {
		/*
		 * Compare the group weights. If a task is all by itself
		 * (not part of a group), use the task weight instead.
		 */
		if (cur_ng && p_ng)
			imp += group_weight(cur, env->src_nid, dist) -
			       group_weight(cur, env->dst_nid, dist);
		else
			imp += task_weight(cur, env->src_nid, dist) -
			       task_weight(cur, env->dst_nid, dist);
	}

	if (maymove && moveimp > imp && moveimp > env->best_imp) {
		imp = moveimp;
		cur = NULL;
		goto assign;
	}

	/*
	 * If the NUMA importance is less than SMALLIMP,
	 * task migration might only result in ping pong
	 * of tasks and also hurt performance due to cache
	 * misses.
	 */
	if (imp < SMALLIMP || imp <= env->best_imp + SMALLIMP / 2)
		goto unlock;

	/*
	 * In the overloaded case, try and keep the load balanced.
	 */
	load = task_h_load(env->p) - task_h_load(cur);
	if (!load)
		goto assign;

	dst_load = env->dst_stats.load + load;
	src_load = env->src_stats.load - load;

	if (load_too_imbalanced(src_load, dst_load, env))
		goto unlock;

assign:
	/*
	 * One idle CPU per node is evaluated for a task numa move.
	 * Call select_idle_sibling to maybe find a better one.
	 */
	if (!cur) {
		/*
		 * select_idle_siblings() uses an per-CPU cpumask that
		 * can be used from IRQ context.
		 */
		local_irq_disable();
		env->dst_cpu = select_idle_sibling(env->p, env->src_cpu,
						   env->dst_cpu);
		local_irq_enable();
	}

	task_numa_assign(env, cur, imp);
unlock:
	rcu_read_unlock();
}

static void task_numa_find_cpu(struct task_numa_env *env,
				long taskimp, long groupimp)
{
	long src_load, dst_load, load;
	bool maymove = false;
	int cpu;

	load = task_h_load(env->p);
	dst_load = env->dst_stats.load + load;
	src_load = env->src_stats.load - load;

	/*
	 * If the improvement from just moving env->p direction is better
	 * than swapping tasks around, check if a move is possible.
	 */
	maymove = !load_too_imbalanced(src_load, dst_load, env);

	for_each_cpu(cpu, cpumask_of_node(env->dst_nid)) {
		/* Skip this CPU if the source task cannot migrate */
		if (!cpumask_test_cpu(cpu, &env->p->cpus_allowed))
			continue;

		env->dst_cpu = cpu;
		task_numa_compare(env, taskimp, groupimp, maymove);
	}
}

static int task_numa_migrate(struct task_struct *p)
{
	struct task_numa_env env = {
		.p = p,

		.src_cpu = task_cpu(p),
		.src_nid = task_node(p),

		.imbalance_pct = 112,

		.best_task = NULL,
		.best_imp = 0,
		.best_cpu = -1,
	};
	unsigned long taskweight, groupweight;
	struct sched_domain *sd;
	long taskimp, groupimp;
	struct numa_group *ng;
	struct rq *best_rq;
	int nid, ret, dist;

	/*
	 * Pick the lowest SD_NUMA domain, as that would have the smallest
	 * imbalance and would be the first to start moving tasks about.
	 *
	 * And we want to avoid any moving of tasks about, as that would create
	 * random movement of tasks -- counter the numa conditions we're trying
	 * to satisfy here.
	 */
	rcu_read_lock();
	sd = rcu_dereference(per_cpu(sd_numa, env.src_cpu));
	if (sd)
		env.imbalance_pct = 100 + (sd->imbalance_pct - 100) / 2;
	rcu_read_unlock();

	/*
	 * Cpusets can break the scheduler domain tree into smaller
	 * balance domains, some of which do not cross NUMA boundaries.
	 * Tasks that are "trapped" in such domains cannot be migrated
	 * elsewhere, so there is no point in (re)trying.
	 */
	if (unlikely(!sd)) {
		sched_setnuma(p, task_node(p));
		return -EINVAL;
	}

	env.dst_nid = p->numa_preferred_nid;
	dist = env.dist = node_distance(env.src_nid, env.dst_nid);
	taskweight = task_weight(p, env.src_nid, dist);
	groupweight = group_weight(p, env.src_nid, dist);
	update_numa_stats(&env.src_stats, env.src_nid);
	taskimp = task_weight(p, env.dst_nid, dist) - taskweight;
	groupimp = group_weight(p, env.dst_nid, dist) - groupweight;
	update_numa_stats(&env.dst_stats, env.dst_nid);

	/* Try to find a spot on the preferred nid. */
	task_numa_find_cpu(&env, taskimp, groupimp);

	/*
	 * Look at other nodes in these cases:
	 * - there is no space available on the preferred_nid
	 * - the task is part of a numa_group that is interleaved across
	 *   multiple NUMA nodes; in order to better consolidate the group,
	 *   we need to check other locations.
	 */
	ng = deref_curr_numa_group(p);
	if (env.best_cpu == -1 || (ng && ng->active_nodes > 1)) {
		for_each_online_node(nid) {
			if (nid == env.src_nid || nid == p->numa_preferred_nid)
				continue;

			dist = node_distance(env.src_nid, env.dst_nid);
			if (sched_numa_topology_type == NUMA_BACKPLANE &&
						dist != env.dist) {
				taskweight = task_weight(p, env.src_nid, dist);
				groupweight = group_weight(p, env.src_nid, dist);
			}

			/* Only consider nodes where both task and groups benefit */
			taskimp = task_weight(p, nid, dist) - taskweight;
			groupimp = group_weight(p, nid, dist) - groupweight;
			if (taskimp < 0 && groupimp < 0)
				continue;

			env.dist = dist;
			env.dst_nid = nid;
			update_numa_stats(&env.dst_stats, env.dst_nid);
			task_numa_find_cpu(&env, taskimp, groupimp);
		}
	}

	/*
	 * If the task is part of a workload that spans multiple NUMA nodes,
	 * and is migrating into one of the workload's active nodes, remember
	 * this node as the task's preferred numa node, so the workload can
	 * settle down.
	 * A task that migrated to a second choice node will be better off
	 * trying for a better one later. Do not set the preferred node here.
	 */
	if (ng) {
		if (env.best_cpu == -1)
			nid = env.src_nid;
		else
			nid = cpu_to_node(env.best_cpu);

		if (nid != p->numa_preferred_nid)
			sched_setnuma(p, nid);
	}

	/* No better CPU than the current one was found. */
	if (env.best_cpu == -1)
		return -EAGAIN;

	best_rq = cpu_rq(env.best_cpu);
	if (env.best_task == NULL) {
		ret = migrate_task_to(p, env.best_cpu);
		WRITE_ONCE(best_rq->numa_migrate_on, 0);
		if (ret != 0)
			trace_sched_stick_numa(p, env.src_cpu, env.best_cpu);
		return ret;
	}

	ret = migrate_swap(p, env.best_task, env.best_cpu, env.src_cpu);
	WRITE_ONCE(best_rq->numa_migrate_on, 0);

	if (ret != 0)
		trace_sched_stick_numa(p, env.src_cpu, task_cpu(env.best_task));
	put_task_struct(env.best_task);
	return ret;
}

/* Attempt to migrate a task to a CPU on the preferred node. */
static void numa_migrate_preferred(struct task_struct *p)
{
	unsigned long interval = HZ;

	/* This task has no NUMA fault statistics yet */
	if (unlikely(p->numa_preferred_nid == -1 || !p->numa_faults))
		return;

	/* Periodically retry migrating the task to the preferred node */
	interval = min(interval, msecs_to_jiffies(p->numa_scan_period) / 16);
	p->numa_migrate_retry = jiffies + interval;

	/* Success if task is already running on preferred CPU */
	if (task_node(p) == p->numa_preferred_nid)
		return;

	/* Otherwise, try migrate to a CPU on the preferred node */
	task_numa_migrate(p);
}

/*
 * Find out how many nodes on the workload is actively running on. Do this by
 * tracking the nodes from which NUMA hinting faults are triggered. This can
 * be different from the set of nodes where the workload's memory is currently
 * located.
 */
static void numa_group_count_active_nodes(struct numa_group *numa_group)
{
	unsigned long faults, max_faults = 0;
	int nid, active_nodes = 0;

	for_each_online_node(nid) {
		faults = group_faults_cpu(numa_group, nid);
		if (faults > max_faults)
			max_faults = faults;
	}

	for_each_online_node(nid) {
		faults = group_faults_cpu(numa_group, nid);
		if (faults * ACTIVE_NODE_FRACTION > max_faults)
			active_nodes++;
	}

	numa_group->max_faults_cpu = max_faults;
	numa_group->active_nodes = active_nodes;
}

/*
 * When adapting the scan rate, the period is divided into NUMA_PERIOD_SLOTS
 * increments. The more local the fault statistics are, the higher the scan
 * period will be for the next scan window. If local/(local+remote) ratio is
 * below NUMA_PERIOD_THRESHOLD (where range of ratio is 1..NUMA_PERIOD_SLOTS)
 * the scan period will decrease. Aim for 70% local accesses.
 */
#define NUMA_PERIOD_SLOTS 10
#define NUMA_PERIOD_THRESHOLD 7

/*
 * Increase the scan period (slow down scanning) if the majority of
 * our memory is already on our local node, or if the majority of
 * the page accesses are shared with other processes.
 * Otherwise, decrease the scan period.
 */
static void update_task_scan_period(struct task_struct *p,
			unsigned long shared, unsigned long private)
{
	unsigned int period_slot;
	int lr_ratio, ps_ratio;
	int diff;

	unsigned long remote = p->numa_faults_locality[0];
	unsigned long local = p->numa_faults_locality[1];

	/*
	 * If there were no record hinting faults then either the task is
	 * completely idle or all activity is areas that are not of interest
	 * to automatic numa balancing. Related to that, if there were failed
	 * migration then it implies we are migrating too quickly or the local
	 * node is overloaded. In either case, scan slower
	 */
	if (local + shared == 0 || p->numa_faults_locality[2]) {
		p->numa_scan_period = min(p->numa_scan_period_max,
			p->numa_scan_period << 1);

		p->mm->numa_next_scan = jiffies +
			msecs_to_jiffies(p->numa_scan_period);

		return;
	}

	/*
	 * Prepare to scale scan period relative to the current period.
	 *	 == NUMA_PERIOD_THRESHOLD scan period stays the same
	 *       <  NUMA_PERIOD_THRESHOLD scan period decreases (scan faster)
	 *	 >= NUMA_PERIOD_THRESHOLD scan period increases (scan slower)
	 */
	period_slot = DIV_ROUND_UP(p->numa_scan_period, NUMA_PERIOD_SLOTS);
	lr_ratio = (local * NUMA_PERIOD_SLOTS) / (local + remote);
	ps_ratio = (private * NUMA_PERIOD_SLOTS) / (private + shared);

	if (ps_ratio >= NUMA_PERIOD_THRESHOLD) {
		/*
		 * Most memory accesses are local. There is no need to
		 * do fast NUMA scanning, since memory is already local.
		 */
		int slot = ps_ratio - NUMA_PERIOD_THRESHOLD;
		if (!slot)
			slot = 1;
		diff = slot * period_slot;
	} else if (lr_ratio >= NUMA_PERIOD_THRESHOLD) {
		/*
		 * Most memory accesses are shared with other tasks.
		 * There is no point in continuing fast NUMA scanning,
		 * since other tasks may just move the memory elsewhere.
		 */
		int slot = lr_ratio - NUMA_PERIOD_THRESHOLD;
		if (!slot)
			slot = 1;
		diff = slot * period_slot;
	} else {
		/*
		 * Private memory faults exceed (SLOTS-THRESHOLD)/SLOTS,
		 * yet they are not on the local NUMA node. Speed up
		 * NUMA scanning to get the memory moved over.
		 */
		int ratio = max(lr_ratio, ps_ratio);
		diff = -(NUMA_PERIOD_THRESHOLD - ratio) * period_slot;
	}

	p->numa_scan_period = clamp(p->numa_scan_period + diff,
			task_scan_min(p), task_scan_max(p));
	memset(p->numa_faults_locality, 0, sizeof(p->numa_faults_locality));
}

/*
 * Get the fraction of time the task has been running since the last
 * NUMA placement cycle. The scheduler keeps similar statistics, but
 * decays those on a 32ms period, which is orders of magnitude off
 * from the dozens-of-seconds NUMA balancing period. Use the scheduler
 * stats only if the task is so new there are no NUMA statistics yet.
 */
static u64 numa_get_avg_runtime(struct task_struct *p, u64 *period)
{
	u64 runtime, delta, now;
	/* Use the start of this time slice to avoid calculations. */
	now = p->se.exec_start;
	runtime = p->se.sum_exec_runtime;

	if (p->last_task_numa_placement) {
		delta = runtime - p->last_sum_exec_runtime;
		*period = now - p->last_task_numa_placement;

		/* Avoid time going backwards, prevent potential divide error: */
		if (unlikely((s64)*period < 0))
			*period = 0;
	} else {
		delta = p->se.avg.load_sum;
		*period = LOAD_AVG_MAX;
	}

	p->last_sum_exec_runtime = runtime;
	p->last_task_numa_placement = now;

	return delta;
}

/*
 * Determine the preferred nid for a task in a numa_group. This needs to
 * be done in a way that produces consistent results with group_weight,
 * otherwise workloads might not converge.
 */
static int preferred_group_nid(struct task_struct *p, int nid)
{
	nodemask_t nodes;
	int dist;

	/* Direct connections between all NUMA nodes. */
	if (sched_numa_topology_type == NUMA_DIRECT)
		return nid;

	/*
	 * On a system with glueless mesh NUMA topology, group_weight
	 * scores nodes according to the number of NUMA hinting faults on
	 * both the node itself, and on nearby nodes.
	 */
	if (sched_numa_topology_type == NUMA_GLUELESS_MESH) {
		unsigned long score, max_score = 0;
		int node, max_node = nid;

		dist = sched_max_numa_distance;

		for_each_online_node(node) {
			score = group_weight(p, node, dist);
			if (score > max_score) {
				max_score = score;
				max_node = node;
			}
		}
		return max_node;
	}

	/*
	 * Finding the preferred nid in a system with NUMA backplane
	 * interconnect topology is more involved. The goal is to locate
	 * tasks from numa_groups near each other in the system, and
	 * untangle workloads from different sides of the system. This requires
	 * searching down the hierarchy of node groups, recursively searching
	 * inside the highest scoring group of nodes. The nodemask tricks
	 * keep the complexity of the search down.
	 */
	nodes = node_online_map;
	for (dist = sched_max_numa_distance; dist > LOCAL_DISTANCE; dist--) {
		unsigned long max_faults = 0;
		nodemask_t max_group = NODE_MASK_NONE;
		int a, b;

		/* Are there nodes at this distance from each other? */
		if (!find_numa_distance(dist))
			continue;

		for_each_node_mask(a, nodes) {
			unsigned long faults = 0;
			nodemask_t this_group;
			nodes_clear(this_group);

			/* Sum group's NUMA faults; includes a==b case. */
			for_each_node_mask(b, nodes) {
				if (node_distance(a, b) < dist) {
					faults += group_faults(p, b);
					node_set(b, this_group);
					node_clear(b, nodes);
				}
			}

			/* Remember the top group. */
			if (faults > max_faults) {
				max_faults = faults;
				max_group = this_group;
				/*
				 * subtle: at the smallest distance there is
				 * just one node left in each "group", the
				 * winner is the preferred nid.
				 */
				nid = a;
			}
		}
		/* Next round, evaluate the nodes within max_group. */
		if (!max_faults)
			break;
		nodes = max_group;
	}
	return nid;
}

static void task_numa_placement(struct task_struct *p)
{
	int seq, nid, max_nid = -1;
	unsigned long max_faults = 0;
	unsigned long fault_types[2] = { 0, 0 };
	unsigned long total_faults;
	u64 runtime, period;
	spinlock_t *group_lock = NULL;
	struct numa_group *ng;

	/*
	 * The p->mm->numa_scan_seq field gets updated without
	 * exclusive access. Use READ_ONCE() here to ensure
	 * that the field is read in a single access:
	 */
	seq = READ_ONCE(p->mm->numa_scan_seq);
	if (p->numa_scan_seq == seq)
		return;
	p->numa_scan_seq = seq;
	p->numa_scan_period_max = task_scan_max(p);

	total_faults = p->numa_faults_locality[0] +
		       p->numa_faults_locality[1];
	runtime = numa_get_avg_runtime(p, &period);

	/* If the task is part of a group prevent parallel updates to group stats */
	ng = deref_curr_numa_group(p);
	if (ng) {
		group_lock = &ng->lock;
		spin_lock_irq(group_lock);
	}

	/* Find the node with the highest number of faults */
	for_each_online_node(nid) {
		/* Keep track of the offsets in numa_faults array */
		int mem_idx, membuf_idx, cpu_idx, cpubuf_idx;
		unsigned long faults = 0, group_faults = 0;
		int priv;

		for (priv = 0; priv < NR_NUMA_HINT_FAULT_TYPES; priv++) {
			long diff, f_diff, f_weight;

			mem_idx = task_faults_idx(NUMA_MEM, nid, priv);
			membuf_idx = task_faults_idx(NUMA_MEMBUF, nid, priv);
			cpu_idx = task_faults_idx(NUMA_CPU, nid, priv);
			cpubuf_idx = task_faults_idx(NUMA_CPUBUF, nid, priv);

			/* Decay existing window, copy faults since last scan */
			diff = p->numa_faults[membuf_idx] - p->numa_faults[mem_idx] / 2;
			fault_types[priv] += p->numa_faults[membuf_idx];
			p->numa_faults[membuf_idx] = 0;

			/*
			 * Normalize the faults_from, so all tasks in a group
			 * count according to CPU use, instead of by the raw
			 * number of faults. Tasks with little runtime have
			 * little over-all impact on throughput, and thus their
			 * faults are less important.
			 */
			f_weight = div64_u64(runtime << 16, period + 1);
			f_weight = (f_weight * p->numa_faults[cpubuf_idx]) /
				   (total_faults + 1);
			f_diff = f_weight - p->numa_faults[cpu_idx] / 2;
			p->numa_faults[cpubuf_idx] = 0;

			p->numa_faults[mem_idx] += diff;
			p->numa_faults[cpu_idx] += f_diff;
			faults += p->numa_faults[mem_idx];
			p->total_numa_faults += diff;
			if (ng) {
				/*
				 * safe because we can only change our own group
				 *
				 * mem_idx represents the offset for a given
				 * nid and priv in a specific region because it
				 * is at the beginning of the numa_faults array.
				 */
				ng->faults[mem_idx] += diff;
				ng->faults_cpu[mem_idx] += f_diff;
				ng->total_faults += diff;
				group_faults += ng->faults[mem_idx];
			}
		}

		if (!ng) {
			if (faults > max_faults) {
				max_faults = faults;
				max_nid = nid;
			}
		} else if (group_faults > max_faults) {
			max_faults = group_faults;
			max_nid = nid;
		}
	}

	if (ng) {
		numa_group_count_active_nodes(ng);
		spin_unlock_irq(group_lock);
		max_nid = preferred_group_nid(p, max_nid);
	}

	if (max_faults) {
		/* Set the new preferred node */
		if (max_nid != p->numa_preferred_nid)
			sched_setnuma(p, max_nid);
	}

	update_task_scan_period(p, fault_types[0], fault_types[1]);
}

static inline int get_numa_group(struct numa_group *grp)
{
	return atomic_inc_not_zero(&grp->refcount);
}

static inline void put_numa_group(struct numa_group *grp)
{
	if (atomic_dec_and_test(&grp->refcount))
		kfree_rcu(grp, rcu);
}

static void task_numa_group(struct task_struct *p, int cpupid, int flags,
			int *priv)
{
	struct numa_group *grp, *my_grp;
	struct task_struct *tsk;
	bool join = false;
	int cpu = cpupid_to_cpu(cpupid);
	int i;

	if (unlikely(!deref_curr_numa_group(p))) {
		unsigned int size = sizeof(struct numa_group) +
				    4*nr_node_ids*sizeof(unsigned long);

		grp = kzalloc(size, GFP_KERNEL | __GFP_NOWARN);
		if (!grp)
			return;

		atomic_set(&grp->refcount, 1);
		grp->active_nodes = 1;
		grp->max_faults_cpu = 0;
		spin_lock_init(&grp->lock);
		grp->gid = p->pid;
		/* Second half of the array tracks nids where faults happen */
		grp->faults_cpu = grp->faults + NR_NUMA_HINT_FAULT_TYPES *
						nr_node_ids;

		for (i = 0; i < NR_NUMA_HINT_FAULT_STATS * nr_node_ids; i++)
			grp->faults[i] = p->numa_faults[i];

		grp->total_faults = p->total_numa_faults;

		grp->nr_tasks++;
		rcu_assign_pointer(p->numa_group, grp);
	}

	rcu_read_lock();
	tsk = READ_ONCE(cpu_rq(cpu)->curr);

	if (!cpupid_match_pid(tsk, cpupid))
		goto no_join;

	grp = rcu_dereference(tsk->numa_group);
	if (!grp)
		goto no_join;

	my_grp = deref_curr_numa_group(p);
	if (grp == my_grp)
		goto no_join;

	/*
	 * Only join the other group if its bigger; if we're the bigger group,
	 * the other task will join us.
	 */
	if (my_grp->nr_tasks > grp->nr_tasks)
		goto no_join;

	/*
	 * Tie-break on the grp address.
	 */
	if (my_grp->nr_tasks == grp->nr_tasks && my_grp > grp)
		goto no_join;

	/* Always join threads in the same process. */
	if (tsk->mm == current->mm)
		join = true;

	/* Simple filter to avoid false positives due to PID collisions */
	if (flags & TNF_SHARED)
		join = true;

	/* Update priv based on whether false sharing was detected */
	*priv = !join;

	if (join && !get_numa_group(grp))
		goto no_join;

	rcu_read_unlock();

	if (!join)
		return;

	BUG_ON(irqs_disabled());
	double_lock_irq(&my_grp->lock, &grp->lock);

	for (i = 0; i < NR_NUMA_HINT_FAULT_STATS * nr_node_ids; i++) {
		my_grp->faults[i] -= p->numa_faults[i];
		grp->faults[i] += p->numa_faults[i];
	}
	my_grp->total_faults -= p->total_numa_faults;
	grp->total_faults += p->total_numa_faults;

	my_grp->nr_tasks--;
	grp->nr_tasks++;

	spin_unlock(&my_grp->lock);
	spin_unlock_irq(&grp->lock);

	rcu_assign_pointer(p->numa_group, grp);

	put_numa_group(my_grp);
	return;

no_join:
	rcu_read_unlock();
	return;
}

/*
 * Get rid of NUMA staticstics associated with a task (either current or dead).
 * If @final is set, the task is dead and has reached refcount zero, so we can
 * safely free all relevant data structures. Otherwise, there might be
 * concurrent reads from places like load balancing and procfs, and we should
 * reset the data back to default state without freeing ->numa_faults.
 */
void task_numa_free(struct task_struct *p, bool final)
{
	/* safe: p either is current or is being freed by current */
	struct numa_group *grp = rcu_dereference_raw(p->numa_group);
	unsigned long *numa_faults = p->numa_faults;
	unsigned long flags;
	int i;

	if (!numa_faults)
		return;

	if (grp) {
		spin_lock_irqsave(&grp->lock, flags);
		for (i = 0; i < NR_NUMA_HINT_FAULT_STATS * nr_node_ids; i++)
			grp->faults[i] -= p->numa_faults[i];
		grp->total_faults -= p->total_numa_faults;

		grp->nr_tasks--;
		spin_unlock_irqrestore(&grp->lock, flags);
		RCU_INIT_POINTER(p->numa_group, NULL);
		put_numa_group(grp);
	}

	if (final) {
		p->numa_faults = NULL;
		kfree(numa_faults);
	} else {
		p->total_numa_faults = 0;
		for (i = 0; i < NR_NUMA_HINT_FAULT_STATS * nr_node_ids; i++)
			numa_faults[i] = 0;
	}
}

/*
 * Got a PROT_NONE fault for a page on @node.
 */
void task_numa_fault(int last_cpupid, int mem_node, int pages, int flags)
{
	struct task_struct *p = current;
	bool migrated = flags & TNF_MIGRATED;
	int cpu_node = task_node(current);
	int local = !!(flags & TNF_FAULT_LOCAL);
	struct numa_group *ng;
	int priv;

	if (!static_branch_likely(&sched_numa_balancing))
		return;

	/* for example, ksmd faulting in a user's mm */
	if (!p->mm)
		return;

	/* Allocate buffer to track faults on a per-node basis */
	if (unlikely(!p->numa_faults)) {
		int size = sizeof(*p->numa_faults) *
			   NR_NUMA_HINT_FAULT_BUCKETS * nr_node_ids;

		p->numa_faults = kzalloc(size, GFP_KERNEL|__GFP_NOWARN);
		if (!p->numa_faults)
			return;

		p->total_numa_faults = 0;
		memset(p->numa_faults_locality, 0, sizeof(p->numa_faults_locality));
	}

	/*
	 * First accesses are treated as private, otherwise consider accesses
	 * to be private if the accessing pid has not changed
	 */
	if (unlikely(last_cpupid == (-1 & LAST_CPUPID_MASK))) {
		priv = 1;
	} else {
		priv = cpupid_match_pid(p, last_cpupid);
		if (!priv && !(flags & TNF_NO_GROUP))
			task_numa_group(p, last_cpupid, flags, &priv);
	}

	/*
	 * If a workload spans multiple NUMA nodes, a shared fault that
	 * occurs wholly within the set of nodes that the workload is
	 * actively using should be counted as local. This allows the
	 * scan rate to slow down when a workload has settled down.
	 */
	ng = deref_curr_numa_group(p);
	if (!priv && !local && ng && ng->active_nodes > 1 &&
				numa_is_active_node(cpu_node, ng) &&
				numa_is_active_node(mem_node, ng))
		local = 1;

	/*
	 * Retry task to preferred node migration periodically, in case it
	 * case it previously failed, or the scheduler moved us.
	 */
	if (time_after(jiffies, p->numa_migrate_retry)) {
		task_numa_placement(p);
		numa_migrate_preferred(p);
	}

	if (migrated)
		p->numa_pages_migrated += pages;
	if (flags & TNF_MIGRATE_FAIL)
		p->numa_faults_locality[2] += pages;

	p->numa_faults[task_faults_idx(NUMA_MEMBUF, mem_node, priv)] += pages;
	p->numa_faults[task_faults_idx(NUMA_CPUBUF, cpu_node, priv)] += pages;
	p->numa_faults_locality[local] += pages;
}

static void reset_ptenuma_scan(struct task_struct *p)
{
	/*
	 * We only did a read acquisition of the mmap sem, so
	 * p->mm->numa_scan_seq is written to without exclusive access
	 * and the update is not guaranteed to be atomic. That's not
	 * much of an issue though, since this is just used for
	 * statistical sampling. Use READ_ONCE/WRITE_ONCE, which are not
	 * expensive, to avoid any form of compiler optimizations:
	 */
	WRITE_ONCE(p->mm->numa_scan_seq, READ_ONCE(p->mm->numa_scan_seq) + 1);
	p->mm->numa_scan_offset = 0;
}

/*
 * The expensive part of numa migration is done from task_work context.
 * Triggered from task_tick_numa().
 */
void task_numa_work(struct callback_head *work)
{
	unsigned long migrate, next_scan, now = jiffies;
	struct task_struct *p = current;
	struct mm_struct *mm = p->mm;
	u64 runtime = p->se.sum_exec_runtime;
	struct vm_area_struct *vma;
	unsigned long start, end;
	unsigned long nr_pte_updates = 0;
	long pages, virtpages;

	SCHED_WARN_ON(p != container_of(work, struct task_struct, numa_work));

	work->next = work; /* protect against double add */
	/*
	 * Who cares about NUMA placement when they're dying.
	 *
	 * NOTE: make sure not to dereference p->mm before this check,
	 * exit_task_work() happens _after_ exit_mm() so we could be called
	 * without p->mm even though we still had it when we enqueued this
	 * work.
	 */
	if (p->flags & PF_EXITING)
		return;

	if (!mm->numa_next_scan) {
		mm->numa_next_scan = now +
			msecs_to_jiffies(sysctl_numa_balancing_scan_delay);
	}

	/*
	 * Enforce maximal scan/migration frequency..
	 */
	migrate = mm->numa_next_scan;
	if (time_before(now, migrate))
		return;

	if (p->numa_scan_period == 0) {
		p->numa_scan_period_max = task_scan_max(p);
		p->numa_scan_period = task_scan_start(p);
	}

	next_scan = now + msecs_to_jiffies(p->numa_scan_period);
	if (cmpxchg(&mm->numa_next_scan, migrate, next_scan) != migrate)
		return;

	/*
	 * Delay this task enough that another task of this mm will likely win
	 * the next time around.
	 */
	p->node_stamp += 2 * TICK_NSEC;

	start = mm->numa_scan_offset;
	pages = sysctl_numa_balancing_scan_size;
	pages <<= 20 - PAGE_SHIFT; /* MB in pages */
	virtpages = pages * 8;	   /* Scan up to this much virtual space */
	if (!pages)
		return;


	if (!down_read_trylock(&mm->mmap_sem))
		return;
	vma = find_vma(mm, start);
	if (!vma) {
		reset_ptenuma_scan(p);
		start = 0;
		vma = mm->mmap;
	}
	for (; vma; vma = vma->vm_next) {
		if (!vma_migratable(vma) || !vma_policy_mof(vma) ||
			is_vm_hugetlb_page(vma) || (vma->vm_flags & VM_MIXEDMAP)) {
			continue;
		}

		/*
		 * Shared library pages mapped by multiple processes are not
		 * migrated as it is expected they are cache replicated. Avoid
		 * hinting faults in read-only file-backed mappings or the vdso
		 * as migrating the pages will be of marginal benefit.
		 */
		if (!vma->vm_mm ||
		    (vma->vm_file && (vma->vm_flags & (VM_READ|VM_WRITE)) == (VM_READ)))
			continue;

		/*
		 * Skip inaccessible VMAs to avoid any confusion between
		 * PROT_NONE and NUMA hinting ptes
		 */
		if (!vma_is_accessible(vma))
			continue;

		do {
			start = max(start, vma->vm_start);
			end = ALIGN(start + (pages << PAGE_SHIFT), HPAGE_SIZE);
			end = min(end, vma->vm_end);
			nr_pte_updates = change_prot_numa(vma, start, end);

			/*
			 * Try to scan sysctl_numa_balancing_size worth of
			 * hpages that have at least one present PTE that
			 * is not already pte-numa. If the VMA contains
			 * areas that are unused or already full of prot_numa
			 * PTEs, scan up to virtpages, to skip through those
			 * areas faster.
			 */
			if (nr_pte_updates)
				pages -= (end - start) >> PAGE_SHIFT;
			virtpages -= (end - start) >> PAGE_SHIFT;

			start = end;
			if (pages <= 0 || virtpages <= 0)
				goto out;

			cond_resched();
		} while (end != vma->vm_end);
	}

out:
	/*
	 * It is possible to reach the end of the VMA list but the last few
	 * VMAs are not guaranteed to the vma_migratable. If they are not, we
	 * would find the !migratable VMA on the next scan but not reset the
	 * scanner to the start so check it now.
	 */
	if (vma)
		mm->numa_scan_offset = start;
	else
		reset_ptenuma_scan(p);
	up_read(&mm->mmap_sem);

	/*
	 * Make sure tasks use at least 32x as much time to run other code
	 * than they used here, to limit NUMA PTE scanning overhead to 3% max.
	 * Usually update_task_scan_period slows down scanning enough; on an
	 * overloaded system we need to limit overhead on a per task basis.
	 */
	if (unlikely(p->se.sum_exec_runtime != runtime)) {
		u64 diff = p->se.sum_exec_runtime - runtime;
		p->node_stamp += 32 * diff;
	}
}

/*
 * Drive the periodic memory faults..
 */
void task_tick_numa(struct rq *rq, struct task_struct *curr)
{
	struct callback_head *work = &curr->numa_work;
	u64 period, now;

	/*
	 * We don't care about NUMA placement if we don't have memory.
	 */
	if ((curr->flags & (PF_EXITING | PF_KTHREAD)) || work->next != work)
		return;

	/*
	 * Using runtime rather than walltime has the dual advantage that
	 * we (mostly) drive the selection from busy threads and that the
	 * task needs to have done some actual work before we bother with
	 * NUMA placement.
	 */
	now = curr->se.sum_exec_runtime;
	period = (u64)curr->numa_scan_period * NSEC_PER_MSEC;

	if (now > curr->node_stamp + period) {
		if (!curr->node_stamp)
			curr->numa_scan_period = task_scan_start(curr);
		curr->node_stamp += period;

		if (!time_before(jiffies, curr->mm->numa_next_scan)) {
			init_task_work(work, task_numa_work); /* TODO: move this into sched_fork() */
			task_work_add(curr, work, true);
		}
	}
}

static void update_scan_period(struct task_struct *p, int new_cpu)
{
	int src_nid = cpu_to_node(task_cpu(p));
	int dst_nid = cpu_to_node(new_cpu);

	if (!static_branch_likely(&sched_numa_balancing))
		return;

	if (!p->mm || !p->numa_faults || (p->flags & PF_EXITING))
		return;

	if (src_nid == dst_nid)
		return;

	/*
	 * Allow resets if faults have been trapped before one scan
	 * has completed. This is most likely due to a new task that
	 * is pulled cross-node due to wakeups or load balancing.
	 */
	if (p->numa_scan_seq) {
		/*
		 * Avoid scan adjustments if moving to the preferred
		 * node or if the task was not previously running on
		 * the preferred node.
		 */
		if (dst_nid == p->numa_preferred_nid ||
		    (p->numa_preferred_nid != -1 && src_nid != p->numa_preferred_nid))
			return;
	}

	p->numa_scan_period = task_scan_start(p);
}

#else
static void task_tick_numa(struct rq *rq, struct task_struct *curr)
{
}

static inline void account_numa_enqueue(struct rq *rq, struct task_struct *p)
{
}

static inline void account_numa_dequeue(struct rq *rq, struct task_struct *p)
{
}

static inline void update_scan_period(struct task_struct *p, int new_cpu)
{
}

#endif /* CONFIG_NUMA_BALANCING */

static void
account_entity_enqueue(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	update_load_add(&cfs_rq->load, se->load.weight);
	if (!parent_entity(se))
		update_load_add(&rq_of(cfs_rq)->load, se->load.weight);
#ifdef CONFIG_SMP
	if (entity_is_task(se)) {
		struct rq *rq = rq_of(cfs_rq);

		account_numa_enqueue(rq, task_of(se));
		list_add(&se->group_node, &rq->cfs_tasks);
	}
#endif
	cfs_rq->nr_running++;
}

static void
account_entity_dequeue(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	update_load_sub(&cfs_rq->load, se->load.weight);
	if (!parent_entity(se))
		update_load_sub(&rq_of(cfs_rq)->load, se->load.weight);
#ifdef CONFIG_SMP
	if (entity_is_task(se)) {
		account_numa_dequeue(rq_of(cfs_rq), task_of(se));
		list_del_init(&se->group_node);
	}
#endif
	cfs_rq->nr_running--;
}

/*
 * Signed add and clamp on underflow.
 *
 * Explicitly do a load-store to ensure the intermediate value never hits
 * memory. This allows lockless observations without ever seeing the negative
 * values.
 */
#define add_positive(_ptr, _val) do {                           \
	typeof(_ptr) ptr = (_ptr);                              \
	typeof(_val) val = (_val);                              \
	typeof(*ptr) res, var = READ_ONCE(*ptr);                \
								\
	res = var + val;                                        \
								\
	if (val < 0 && res > var)                               \
		res = 0;                                        \
								\
	WRITE_ONCE(*ptr, res);                                  \
} while (0)

/*
 * Unsigned subtract and clamp on underflow.
 *
 * Explicitly do a load-store to ensure the intermediate value never hits
 * memory. This allows lockless observations without ever seeing the negative
 * values.
 */
#define sub_positive(_ptr, _val) do {				\
	typeof(_ptr) ptr = (_ptr);				\
	typeof(*ptr) val = (_val);				\
	typeof(*ptr) res, var = READ_ONCE(*ptr);		\
	res = var - val;					\
	if (res > var)						\
		res = 0;					\
	WRITE_ONCE(*ptr, res);					\
} while (0)

#ifdef CONFIG_SMP
static inline void
enqueue_runnable_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	cfs_rq->runnable_weight += se->runnable_weight;

	cfs_rq->avg.runnable_load_avg += se->avg.runnable_load_avg;
	cfs_rq->avg.runnable_load_sum += se_runnable(se) * se->avg.runnable_load_sum;
}

static inline void
dequeue_runnable_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	cfs_rq->runnable_weight -= se->runnable_weight;

	sub_positive(&cfs_rq->avg.runnable_load_avg, se->avg.runnable_load_avg);
	sub_positive(&cfs_rq->avg.runnable_load_sum,
		     se_runnable(se) * se->avg.runnable_load_sum);
}

static inline void
enqueue_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	cfs_rq->avg.load_avg += se->avg.load_avg;
	cfs_rq->avg.load_sum += se_weight(se) * se->avg.load_sum;
}

static inline void
dequeue_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	sub_positive(&cfs_rq->avg.load_avg, se->avg.load_avg);
	sub_positive(&cfs_rq->avg.load_sum, se_weight(se) * se->avg.load_sum);
}
#else
static inline void
enqueue_runnable_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se) { }
static inline void
dequeue_runnable_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se) { }
static inline void
enqueue_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se) { }
static inline void
dequeue_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se) { }
#endif

static void reweight_eevdf(struct sched_entity *se, u64 avruntime,
			   unsigned long weight)
{
	unsigned long old_weight = se->load.weight;
	s64 vlag, vslice;

	/*
	 * VRUNTIME
	 * ========
	 *
	 * COROLLARY #1: The virtual runtime of the entity needs to be
	 * adjusted if re-weight at !0-lag point.
	 *
	 * Proof: For contradiction assume this is not true, so we can
	 * re-weight without changing vruntime at !0-lag point.
	 *
	 *             Weight	VRuntime   Avg-VRuntime
	 *     before    w          v            V
	 *      after    w'         v'           V'
	 *
	 * Since lag needs to be preserved through re-weight:
	 *
	 *	lag = (V - v)*w = (V'- v')*w', where v = v'
	 *	==>	V' = (V - v)*w/w' + v		(1)
	 *
	 * Let W be the total weight of the entities before reweight,
	 * since V' is the new weighted average of entities:
	 *
	 *	V' = (WV + w'v - wv) / (W + w' - w)	(2)
	 *
	 * by using (1) & (2) we obtain:
	 *
	 *	(WV + w'v - wv) / (W + w' - w) = (V - v)*w/w' + v
	 *	==> (WV-Wv+Wv+w'v-wv)/(W+w'-w) = (V - v)*w/w' + v
	 *	==> (WV - Wv)/(W + w' - w) + v = (V - v)*w/w' + v
	 *	==>	(V - v)*W/(W + w' - w) = (V - v)*w/w' (3)
	 *
	 * Since we are doing at !0-lag point which means V != v, we
	 * can simplify (3):
	 *
	 *	==>	W / (W + w' - w) = w / w'
	 *	==>	Ww' = Ww + ww' - ww
	 *	==>	W * (w' - w) = w * (w' - w)
	 *	==>	W = w	(re-weight indicates w' != w)
	 *
	 * So the cfs_rq contains only one entity, hence vruntime of
	 * the entity @v should always equal to the cfs_rq's weighted
	 * average vruntime @V, which means we will always re-weight
	 * at 0-lag point, thus breach assumption. Proof completed.
	 *
	 *
	 * COROLLARY #2: Re-weight does NOT affect weighted average
	 * vruntime of all the entities.
	 *
	 * Proof: According to corollary #1, Eq. (1) should be:
	 *
	 *	(V - v)*w = (V' - v')*w'
	 *	==>    v' = V' - (V - v)*w/w'		(4)
	 *
	 * According to the weighted average formula, we have:
	 *
	 *	V' = (WV - wv + w'v') / (W - w + w')
	 *	   = (WV - wv + w'(V' - (V - v)w/w')) / (W - w + w')
	 *	   = (WV - wv + w'V' - Vw + wv) / (W - w + w')
	 *	   = (WV + w'V' - Vw) / (W - w + w')
	 *
	 *	==>  V'*(W - w + w') = WV + w'V' - Vw
	 *	==>	V' * (W - w) = (W - w) * V	(5)
	 *
	 * If the entity is the only one in the cfs_rq, then reweight
	 * always occurs at 0-lag point, so V won't change. Or else
	 * there are other entities, hence W != w, then Eq. (5) turns
	 * into V' = V. So V won't change in either case, proof done.
	 *
	 *
	 * So according to corollary #1 & #2, the effect of re-weight
	 * on vruntime should be:
	 *
	 *	v' = V' - (V - v) * w / w'		(4)
	 *	   = V  - (V - v) * w / w'
	 *	   = V  - vl * w / w'
	 *	   = V  - vl'
	 */
	if (avruntime != se->vruntime) {
		vlag = entity_lag(avruntime, se);
		vlag = div_s64(vlag * old_weight, weight);
		se->vruntime = avruntime - vlag;
	}

	/*
	 * DEADLINE
	 * ========
	 *
	 * When the weight changes, the virtual time slope changes and
	 * we should adjust the relative virtual deadline accordingly.
	 *
	 *	d' = v' + (d - v)*w/w'
	 *	   = V' - (V - v)*w/w' + (d - v)*w/w'
	 *	   = V  - (V - v)*w/w' + (d - v)*w/w'
	 *	   = V  + (d - V)*w/w'
	 */
	vslice = (s64)(se->deadline - avruntime);
	vslice = div_s64(vslice * old_weight, weight);
	se->deadline = avruntime + vslice;
}

static void reweight_entity(struct cfs_rq *cfs_rq, struct sched_entity *se,
			    unsigned long weight, unsigned long runnable)
{
	bool curr = cfs_rq->curr == se;
	u64 avruntime;

	if (se->on_rq) {
		/* commit outstanding execution time */
		update_curr(cfs_rq);
		avruntime = avg_vruntime(cfs_rq);
		if (!curr)
			__dequeue_entity(cfs_rq, se);
                update_load_sub(&cfs_rq->load, se->load.weight);
		dequeue_runnable_load_avg(cfs_rq, se);
	}
	dequeue_load_avg(cfs_rq, se);

	se->runnable_weight = runnable;

	if (se->on_rq) {
			reweight_eevdf(se, avruntime, weight);
				} else {
		/*
		 * Because we keep se->vlag = V - v_i, while: lag_i = w_i*(V - v_i),
		 * we need to scale se->vlag when w_i changes.
		 */
		se->vlag = div_s64(se->vlag * se->load.weight, weight);
	}

	update_load_set(&se->load, weight);

#ifdef CONFIG_SMP
	do {
		u32 divider = get_pelt_divider(&se->avg);

		se->avg.load_avg = div_u64(se_weight(se) * se->avg.load_sum, divider);
		se->avg.runnable_load_avg =
			div_u64(se_runnable(se) * se->avg.runnable_load_sum, divider);
	} while (0);
#endif

	enqueue_load_avg(cfs_rq, se);
	if (se->on_rq) {
                update_load_add(&cfs_rq->load, se->load.weight);
		enqueue_runnable_load_avg(cfs_rq, se);
		if (!curr)
			__enqueue_entity(cfs_rq, se);

		/*
		 * The entity's vruntime has been adjusted, so let's check
		 * whether the rq-wide min_vruntime needs updated too. Since
		 * the calculations above require stable min_vruntime rather
		 * than up-to-date one, we do the update at the end of the
		 * reweight process.
		 */
		update_min_vruntime(cfs_rq);
	}
}

void reweight_task(struct task_struct *p, int prio)
{
	struct sched_entity *se = &p->se;
	struct cfs_rq *cfs_rq = cfs_rq_of(se);
	struct load_weight *load = &se->load;
	unsigned long weight = scale_load(sched_prio_to_weight[prio]);

	reweight_entity(cfs_rq, se, weight, weight);
	load->inv_weight = sched_prio_to_wmult[prio];
}

#ifdef CONFIG_FAIR_GROUP_SCHED
#ifdef CONFIG_SMP
/*
 * All this does is approximate the hierarchical proportion which includes that
 * global sum we all love to hate.
 *
 * That is, the weight of a group entity, is the proportional share of the
 * group weight based on the group runqueue weights. That is:
 *
 *                     tg->weight * grq->load.weight
 *   ge->load.weight = -----------------------------               (1)
 *                       \Sum grq->load.weight
 *
 * Now, because computing that sum is prohibitively expensive to compute (been
 * there, done that) we approximate it with this average stuff. The average
 * moves slower and therefore the approximation is cheaper and more stable.
 *
 * So instead of the above, we substitute:
 *
 *   grq->load.weight -> grq->avg.load_avg                         (2)
 *
 * which yields the following:
 *
 *                     tg->weight * grq->avg.load_avg
 *   ge->load.weight = ------------------------------              (3)
 *                             tg->load_avg
 *
 * Where: tg->load_avg ~= \Sum grq->avg.load_avg
 *
 * That is shares_avg, and it is right (given the approximation (2)).
 *
 * The problem with it is that because the average is slow -- it was designed
 * to be exactly that of course -- this leads to transients in boundary
 * conditions. In specific, the case where the group was idle and we start the
 * one task. It takes time for our CPU's grq->avg.load_avg to build up,
 * yielding bad latency etc..
 *
 * Now, in that special case (1) reduces to:
 *
 *                     tg->weight * grq->load.weight
 *   ge->load.weight = ----------------------------- = tg->weight   (4)
 *                         grp->load.weight
 *
 * That is, the sum collapses because all other CPUs are idle; the UP scenario.
 *
 * So what we do is modify our approximation (3) to approach (4) in the (near)
 * UP case, like:
 *
 *   ge->load.weight =
 *
 *              tg->weight * grq->load.weight
 *     ---------------------------------------------------         (5)
 *     tg->load_avg - grq->avg.load_avg + grq->load.weight
 *
 * But because grq->load.weight can drop to 0, resulting in a divide by zero,
 * we need to use grq->avg.load_avg as its lower bound, which then gives:
 *
 *
 *                     tg->weight * grq->load.weight
 *   ge->load.weight = -----------------------------		   (6)
 *                             tg_load_avg'
 *
 * Where:
 *
 *   tg_load_avg' = tg->load_avg - grq->avg.load_avg +
 *                  max(grq->load.weight, grq->avg.load_avg)
 *
 * And that is shares_weight and is icky. In the (near) UP case it approaches
 * (4) while in the normal case it approaches (3). It consistently
 * overestimates the ge->load.weight and therefore:
 *
 *   \Sum ge->load.weight >= tg->weight
 *
 * hence icky!
 */
static long calc_group_shares(struct cfs_rq *cfs_rq)
{
	long tg_weight, tg_shares, load, shares;
	struct task_group *tg = cfs_rq->tg;

	tg_shares = READ_ONCE(tg->shares);

	load = max(scale_load_down(cfs_rq->load.weight), cfs_rq->avg.load_avg);

	tg_weight = atomic_long_read(&tg->load_avg);

	/* Ensure tg_weight >= load */
	tg_weight -= cfs_rq->tg_load_avg_contrib;
	tg_weight += load;

	shares = (tg_shares * load);
	if (tg_weight)
		shares /= tg_weight;

	/*
	 * MIN_SHARES has to be unscaled here to support per-CPU partitioning
	 * of a group with small tg->shares value. It is a floor value which is
	 * assigned as a minimum load.weight to the sched_entity representing
	 * the group on a CPU.
	 *
	 * E.g. on 64-bit for a group with tg->shares of scale_load(15)=15*1024
	 * on an 8-core system with 8 tasks each runnable on one CPU shares has
	 * to be 15*1024*1/8=1920 instead of scale_load(MIN_SHARES)=2*1024. In
	 * case no task is runnable on a CPU MIN_SHARES=2 should be returned
	 * instead of 0.
	 */
	return clamp_t(long, shares, MIN_SHARES, tg_shares);
}

/*
 * This calculates the effective runnable weight for a group entity based on
 * the group entity weight calculated above.
 *
 * Because of the above approximation (2), our group entity weight is
 * an load_avg based ratio (3). This means that it includes blocked load and
 * does not represent the runnable weight.
 *
 * Approximate the group entity's runnable weight per ratio from the group
 * runqueue:
 *
 *					     grq->avg.runnable_load_avg
 *   ge->runnable_weight = ge->load.weight * -------------------------- (7)
 *						 grq->avg.load_avg
 *
 * However, analogous to above, since the avg numbers are slow, this leads to
 * transients in the from-idle case. Instead we use:
 *
 *   ge->runnable_weight = ge->load.weight *
 *
 *		max(grq->avg.runnable_load_avg, grq->runnable_weight)
 *		-----------------------------------------------------	(8)
 *		      max(grq->avg.load_avg, grq->load.weight)
 *
 * Where these max() serve both to use the 'instant' values to fix the slow
 * from-idle and avoid the /0 on to-idle, similar to (6).
 */
static long calc_group_runnable(struct cfs_rq *cfs_rq, long shares)
{
	long runnable, load_avg;

	load_avg = max(cfs_rq->avg.load_avg,
		       scale_load_down(cfs_rq->load.weight));

	runnable = max(cfs_rq->avg.runnable_load_avg,
		       scale_load_down(cfs_rq->runnable_weight));

	runnable *= shares;
	if (load_avg)
		runnable /= load_avg;

	return clamp_t(long, runnable, MIN_SHARES, shares);
}
#endif /* CONFIG_SMP */

static inline int throttled_hierarchy(struct cfs_rq *cfs_rq);

/*
 * Recomputes the group entity based on the current state of its group
 * runqueue.
 */
static void update_cfs_group(struct sched_entity *se)
{
	struct cfs_rq *gcfs_rq = group_cfs_rq(se);
	long shares, runnable;

	if (!gcfs_rq)
		return;

	if (throttled_hierarchy(gcfs_rq))
		return;

#ifndef CONFIG_SMP
	runnable = shares = READ_ONCE(gcfs_rq->tg->shares);

#else
	shares = calc_group_shares(gcfs_rq);
	runnable = calc_group_runnable(gcfs_rq, shares);
#endif

	if (unlikely(se->load.weight != shares))
		reweight_entity(cfs_rq_of(se), se, shares, runnable);
}

#else /* CONFIG_FAIR_GROUP_SCHED */
static inline void update_cfs_group(struct sched_entity *se)
{
}
#endif /* CONFIG_FAIR_GROUP_SCHED */

static inline void cfs_rq_util_change(struct cfs_rq *cfs_rq, int flags)
{
	struct rq *rq = rq_of(cfs_rq);

	if (&rq->cfs == cfs_rq || (flags & SCHED_CPUFREQ_MIGRATION)) {
		/*
		 * There are a few boundary cases this might miss but it should
		 * get called often enough that that should (hopefully) not be
		 * a real problem.
		 *
		 * It will not get called when we go idle, because the idle
		 * thread is a different class (!fair), nor will the utilization
		 * number include things like RT tasks.
		 *
		 * As is, the util number is not freq-invariant (we'd have to
		 * implement arch_scale_freq_capacity() for that).
		 *
		 * See cpu_util().
		 */
		cpufreq_update_util(rq, flags);
	}
}

static inline int per_task_boost(struct task_struct *p)
{
	if (p->boost_period) {
		if (sched_clock() > p->boost_expires) {
			p->boost_period = 0;
			p->boost_expires = 0;
			p->boost = 0;
		}
	}
	return p->boost;
}

#ifdef CONFIG_SMP
static inline u64 cfs_rq_last_update_time(struct cfs_rq *cfs_rq)
{
	return u64_u32_load_copy(cfs_rq->avg.last_update_time,
				 cfs_rq->last_update_time_copy);
}
#ifdef CONFIG_FAIR_GROUP_SCHED
/**
 * update_tg_load_avg - update the tg's load avg
 * @cfs_rq: the cfs_rq whose avg changed
 * @force: update regardless of how small the difference
 *
 * This function 'ensures': tg->load_avg := \Sum tg->cfs_rq[]->avg.load.
 * However, because tg->load_avg is a global value there are performance
 * considerations.
 *
 * In order to avoid having to look at the other cfs_rq's, we use a
 * differential update where we store the last value we propagated. This in
 * turn allows skipping updates if the differential is 'small'.
 *
 * Updating tg's load_avg is necessary before update_cfs_share().
 */
static inline void update_tg_load_avg(struct cfs_rq *cfs_rq, int force)
{
	long delta = cfs_rq->avg.load_avg - cfs_rq->tg_load_avg_contrib;

	/*
	 * No need to update load_avg for root_task_group as it is not used.
	 */
	if (cfs_rq->tg == &root_task_group)
		return;

	if (force || abs(delta) > cfs_rq->tg_load_avg_contrib / 64) {
		atomic_long_add(delta, &cfs_rq->tg->load_avg);
		cfs_rq->tg_load_avg_contrib = cfs_rq->avg.load_avg;

		trace_sched_load_tg(cfs_rq);
	}
}

/*
 * Called within set_task_rq() right before setting a task's CPU. The
 * caller only guarantees p->pi_lock is held; no other assumptions,
 * including the state of rq->lock, should be made.
 */
void set_task_rq_fair(struct sched_entity *se,
		      struct cfs_rq *prev, struct cfs_rq *next)
{
	u64 p_last_update_time;
	u64 n_last_update_time;

	if (!sched_feat(ATTACH_AGE_LOAD))
		return;

	/*
	 * We are supposed to update the task to "current" time, then its up to
	 * date and ready to go to new CPU/cfs_rq. But we have difficulty in
	 * getting what current time is, so simply throw away the out-of-date
	 * time. This will result in the wakee task is less decayed, but giving
	 * the wakee more load sounds not bad.
	 */
	if (!(se->avg.last_update_time && prev))
		return;

	p_last_update_time = cfs_rq_last_update_time(prev);
	n_last_update_time = cfs_rq_last_update_time(next);

	__update_load_avg_blocked_se(p_last_update_time, se);
	se->avg.last_update_time = n_last_update_time;
}


/*
 * When on migration a sched_entity joins/leaves the PELT hierarchy, we need to
 * propagate its contribution. The key to this propagation is the invariant
 * that for each group:
 *
 *   ge->avg == grq->avg						(1)
 *
 * _IFF_ we look at the pure running and runnable sums. Because they
 * represent the very same entity, just at different points in the hierarchy.
 *
 * Per the above update_tg_cfs_util() is trivial and simply copies the running
 * sum over (but still wrong, because the group entity and group rq do not have
 * their PELT windows aligned).
 *
 * However, update_tg_cfs_runnable() is more complex. So we have:
 *
 *   ge->avg.load_avg = ge->load.weight * ge->avg.runnable_avg		(2)
 *
 * And since, like util, the runnable part should be directly transferable,
 * the following would _appear_ to be the straight forward approach:
 *
 *   grq->avg.load_avg = grq->load.weight * grq->avg.runnable_avg	(3)
 *
 * And per (1) we have:
 *
 *   ge->avg.runnable_avg == grq->avg.runnable_avg
 *
 * Which gives:
 *
 *                      ge->load.weight * grq->avg.load_avg
 *   ge->avg.load_avg = -----------------------------------		(4)
 *                               grq->load.weight
 *
 * Except that is wrong!
 *
 * Because while for entities historical weight is not important and we
 * really only care about our future and therefore can consider a pure
 * runnable sum, runqueues can NOT do this.
 *
 * We specifically want runqueues to have a load_avg that includes
 * historical weights. Those represent the blocked load, the load we expect
 * to (shortly) return to us. This only works by keeping the weights as
 * integral part of the sum. We therefore cannot decompose as per (3).
 *
 * Another reason this doesn't work is that runnable isn't a 0-sum entity.
 * Imagine a rq with 2 tasks that each are runnable 2/3 of the time. Then the
 * rq itself is runnable anywhere between 2/3 and 1 depending on how the
 * runnable section of these tasks overlap (or not). If they were to perfectly
 * align the rq as a whole would be runnable 2/3 of the time. If however we
 * always have at least 1 runnable task, the rq as a whole is always runnable.
 *
 * So we'll have to approximate.. :/
 *
 * Given the constraint:
 *
 *   ge->avg.running_sum <= ge->avg.runnable_sum <= LOAD_AVG_MAX
 *
 * We can construct a rule that adds runnable to a rq by assuming minimal
 * overlap.
 *
 * On removal, we'll assume each task is equally runnable; which yields:
 *
 *   grq->avg.runnable_sum = grq->avg.load_sum / grq->load.weight
 *
 * XXX: only do this for the part of runnable > running ?
 *
 */

static inline void
update_tg_cfs_util(struct cfs_rq *cfs_rq, struct sched_entity *se, struct cfs_rq *gcfs_rq)
{
	long delta = gcfs_rq->avg.util_avg - se->avg.util_avg;

	/* Nothing to update */
	if (!delta)
		return;

	/*
	 * The relation between sum and avg is:
	 *
	 *   LOAD_AVG_MAX - 1024 + sa->period_contrib
	 *
	 * however, the PELT windows are not aligned between grq and gse.
	 */

	/* Set new sched_entity's utilization */
	se->avg.util_avg = gcfs_rq->avg.util_avg;
	se->avg.util_sum = se->avg.util_avg * LOAD_AVG_MAX;

	/* Update parent cfs_rq utilization */
	add_positive(&cfs_rq->avg.util_avg, delta);
	cfs_rq->avg.util_sum = cfs_rq->avg.util_avg * LOAD_AVG_MAX;
}

static inline void
update_tg_cfs_runnable(struct cfs_rq *cfs_rq, struct sched_entity *se, struct cfs_rq *gcfs_rq)
{
	long delta_avg, running_sum, runnable_sum = gcfs_rq->prop_runnable_sum;
	unsigned long runnable_load_avg, load_avg;
	u64 runnable_load_sum, load_sum = 0;
	s64 delta_sum;

	if (!runnable_sum)
		return;

	gcfs_rq->prop_runnable_sum = 0;

	if (runnable_sum >= 0) {
		/*
		 * Add runnable; clip at LOAD_AVG_MAX. Reflects that until
		 * the CPU is saturated running == runnable.
		 */
		runnable_sum += se->avg.load_sum;
		runnable_sum = min(runnable_sum, (long)LOAD_AVG_MAX);
	} else {
		/*
		 * Estimate the new unweighted runnable_sum of the gcfs_rq by
		 * assuming all tasks are equally runnable.
		 */
		if (scale_load_down(gcfs_rq->load.weight)) {
			load_sum = div_s64(gcfs_rq->avg.load_sum,
				scale_load_down(gcfs_rq->load.weight));
		}

		/* But make sure to not inflate se's runnable */
		runnable_sum = min(se->avg.load_sum, load_sum);
	}

	/*
	 * runnable_sum can't be lower than running_sum
	 * Rescale running sum to be in the same range as runnable sum
	 * running_sum is in [0 : LOAD_AVG_MAX <<  SCHED_CAPACITY_SHIFT]
	 * runnable_sum is in [0 : LOAD_AVG_MAX]
	 */
	running_sum = se->avg.util_sum >> SCHED_CAPACITY_SHIFT;
	runnable_sum = max(runnable_sum, running_sum);

	load_sum = (s64)se_weight(se) * runnable_sum;
	load_avg = div_s64(load_sum, LOAD_AVG_MAX);

	delta_sum = load_sum - (s64)se_weight(se) * se->avg.load_sum;
	delta_avg = load_avg - se->avg.load_avg;

	se->avg.load_sum = runnable_sum;
	se->avg.load_avg = load_avg;
	add_positive(&cfs_rq->avg.load_avg, delta_avg);
	add_positive(&cfs_rq->avg.load_sum, delta_sum);

	runnable_load_sum = (s64)se_runnable(se) * runnable_sum;
	runnable_load_avg = div_s64(runnable_load_sum, LOAD_AVG_MAX);

	if (se->on_rq) {
		delta_sum = runnable_load_sum -
				se_weight(se) * se->avg.runnable_load_sum;
		delta_avg = runnable_load_avg - se->avg.runnable_load_avg;
		add_positive(&cfs_rq->avg.runnable_load_avg, delta_avg);
		add_positive(&cfs_rq->avg.runnable_load_sum, delta_sum);
	}

	se->avg.runnable_load_sum = runnable_sum;
	se->avg.runnable_load_avg = runnable_load_avg;
}

static inline void add_tg_cfs_propagate(struct cfs_rq *cfs_rq, long runnable_sum)
{
	cfs_rq->propagate = 1;
	cfs_rq->prop_runnable_sum += runnable_sum;
}

/* Update task and its cfs_rq load average */
static inline int propagate_entity_load_avg(struct sched_entity *se)
{
	struct cfs_rq *cfs_rq, *gcfs_rq;

	if (entity_is_task(se))
		return 0;

	gcfs_rq = group_cfs_rq(se);
	if (!gcfs_rq->propagate)
		return 0;

	gcfs_rq->propagate = 0;

	cfs_rq = cfs_rq_of(se);

	add_tg_cfs_propagate(cfs_rq, gcfs_rq->prop_runnable_sum);

	update_tg_cfs_util(cfs_rq, se, gcfs_rq);
	update_tg_cfs_runnable(cfs_rq, se, gcfs_rq);

	trace_sched_load_cfs_rq(cfs_rq);
	trace_sched_load_se(se);

	return 1;
}

/*
 * Check if we need to update the load and the utilization of a blocked
 * group_entity:
 */
static inline bool skip_blocked_update(struct sched_entity *se)
{
	struct cfs_rq *gcfs_rq = group_cfs_rq(se);

	/*
	 * If sched_entity still have not zero load or utilization, we have to
	 * decay it:
	 */
	if (se->avg.load_avg || se->avg.util_avg)
		return false;

	/*
	 * If there is a pending propagation, we have to update the load and
	 * the utilization of the sched_entity:
	 */
	if (gcfs_rq->propagate)
		return false;

	/*
	 * Otherwise, the load and the utilization of the sched_entity is
	 * already zero and there is no pending propagation, so it will be a
	 * waste of time to try to decay it:
	 */
	return true;
}

#else /* CONFIG_FAIR_GROUP_SCHED */

static inline void update_tg_load_avg(struct cfs_rq *cfs_rq, int force) {}

static inline int propagate_entity_load_avg(struct sched_entity *se)
{
	return 0;
}

static inline void add_tg_cfs_propagate(struct cfs_rq *cfs_rq, long runnable_sum) {}

#endif /* CONFIG_FAIR_GROUP_SCHED */

/**
 * update_cfs_rq_load_avg - update the cfs_rq's load/util averages
 * @now: current time, as per cfs_rq_clock_pelt()
 * @cfs_rq: cfs_rq to update
 *
 * The cfs_rq avg is the direct sum of all its entities (blocked and runnable)
 * avg. The immediate corollary is that all (fair) tasks must be attached, see
 * post_init_entity_util_avg().
 *
 * cfs_rq->avg is used for task_h_load() and update_cfs_share() for example.
 *
 * Returns true if the load decayed or we removed load.
 *
 * Since both these conditions indicate a changed cfs_rq->avg.load we should
 * call update_tg_load_avg() when this function returns true.
 */
static inline int
update_cfs_rq_load_avg(u64 now, struct cfs_rq *cfs_rq)
{
	unsigned long removed_load = 0, removed_util = 0, removed_runnable_sum = 0;
	struct sched_avg *sa = &cfs_rq->avg;
	int decayed = 0;

	if (cfs_rq->removed.nr) {
		unsigned long r;
		u32 divider = get_pelt_divider(&cfs_rq->avg);

		raw_spin_lock(&cfs_rq->removed.lock);
		swap(cfs_rq->removed.util_avg, removed_util);
		swap(cfs_rq->removed.load_avg, removed_load);
		swap(cfs_rq->removed.runnable_sum, removed_runnable_sum);
		cfs_rq->removed.nr = 0;
		raw_spin_unlock(&cfs_rq->removed.lock);

		r = removed_load;
		sub_positive(&sa->load_avg, r);
		sub_positive(&sa->load_sum, r * divider);

		r = removed_util;
		sub_positive(&sa->util_avg, r);
		sub_positive(&sa->util_sum, r * divider);

		add_tg_cfs_propagate(cfs_rq, -(long)removed_runnable_sum);

		decayed = 1;
	}

	decayed |= __update_load_avg_cfs_rq(now, cfs_rq);
	u64_u32_store_copy(sa->last_update_time,
			   cfs_rq->last_update_time_copy,
			   sa->last_update_time);
	return decayed;
}

/**
 * attach_entity_load_avg - attach this entity to its cfs_rq load avg
 * @cfs_rq: cfs_rq to attach to
 * @se: sched_entity to attach
 * @flags: migration hints
 *
 * Must call update_cfs_rq_load_avg() before this, since we rely on
 * cfs_rq->avg.last_update_time being current.
 */
static void attach_entity_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se, int flags)
{
	/*
	 * cfs_rq->avg.period_contrib can be used for both cfs_rq and se.
	 * See ___update_load_avg() for details.
	 */
	u32 divider = get_pelt_divider(&cfs_rq->avg);

	/*
	 * When we attach the @se to the @cfs_rq, we must align the decay
	 * window because without that, really weird and wonderful things can
	 * happen.
	 *
	 * XXX illustrate
	 */
	se->avg.last_update_time = cfs_rq->avg.last_update_time;
	se->avg.period_contrib = cfs_rq->avg.period_contrib;

	/*
	 * Hell(o) Nasty stuff.. we need to recompute _sum based on the new
	 * period_contrib. This isn't strictly correct, but since we're
	 * entirely outside of the PELT hierarchy, nobody cares if we truncate
	 * _sum a little.
	 */
	se->avg.util_sum = se->avg.util_avg * divider;

	se->avg.load_sum = se->avg.load_avg * divider;
	if (se_weight(se) < se->avg.load_sum)
		se->avg.load_sum = div_u64(se->avg.load_sum, se_weight(se));
	else
		se->avg.load_sum = 1;

	se->avg.runnable_load_sum = se->avg.load_sum;

	enqueue_load_avg(cfs_rq, se);
	cfs_rq->avg.util_avg += se->avg.util_avg;
	cfs_rq->avg.util_sum += se->avg.util_sum;

	add_tg_cfs_propagate(cfs_rq, se->avg.load_sum);

	cfs_rq_util_change(cfs_rq, flags);

	trace_sched_load_cfs_rq(cfs_rq);
}

/**
 * detach_entity_load_avg - detach this entity from its cfs_rq load avg
 * @cfs_rq: cfs_rq to detach from
 * @se: sched_entity to detach
 *
 * Must call update_cfs_rq_load_avg() before this, since we rely on
 * cfs_rq->avg.last_update_time being current.
 */
static void detach_entity_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	dequeue_load_avg(cfs_rq, se);
	sub_positive(&cfs_rq->avg.util_avg, se->avg.util_avg);
	sub_positive(&cfs_rq->avg.util_sum, se->avg.util_sum);

	add_tg_cfs_propagate(cfs_rq, -se->avg.load_sum);

	cfs_rq_util_change(cfs_rq, 0);

	trace_sched_load_cfs_rq(cfs_rq);
}

/*
 * Optional action to be done while updating the load average
 */
#define UPDATE_TG	0x1
#define SKIP_AGE_LOAD	0x2
#define DO_ATTACH	0x4

/* Update task and its cfs_rq load average */
static inline void update_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se, int flags)
{
	u64 now = cfs_rq_clock_pelt(cfs_rq);
	int decayed;

	/*
	 * Track task load average for carrying it to new CPU after migrated, and
	 * track group sched_entity load average for task_h_load calc in migration
	 */
	if (se->avg.last_update_time && !(flags & SKIP_AGE_LOAD))
		__update_load_avg_se(now, cfs_rq, se);

	decayed  = update_cfs_rq_load_avg(now, cfs_rq);
	decayed |= propagate_entity_load_avg(se);

	if (!se->avg.last_update_time && (flags & DO_ATTACH)) {

		/*
		 * DO_ATTACH means we're here from enqueue_entity().
		 * !last_update_time means we've passed through
		 * migrate_task_rq_fair() indicating we migrated.
		 *
		 * IOW we're enqueueing a task on a new CPU.
		 */
		attach_entity_load_avg(cfs_rq, se, SCHED_CPUFREQ_MIGRATION);
		update_tg_load_avg(cfs_rq, 0);

	} else if (decayed) {
		cfs_rq_util_change(cfs_rq, 0);

		if (flags & UPDATE_TG)
			update_tg_load_avg(cfs_rq, 0);
	}
}

/*
 * Synchronize entity load avg of dequeued entity without locking
 * the previous rq.
 */
void sync_entity_load_avg(struct sched_entity *se)
{
	struct cfs_rq *cfs_rq = cfs_rq_of(se);
	u64 last_update_time;

	last_update_time = cfs_rq_last_update_time(cfs_rq);
	__update_load_avg_blocked_se(last_update_time, se);
}

/*
 * Task first catches up with cfs_rq, and then subtract
 * itself from the cfs_rq (task must be off the queue now).
 */
void remove_entity_load_avg(struct sched_entity *se)
{
	struct cfs_rq *cfs_rq = cfs_rq_of(se);
	unsigned long flags;

	/*
	 * tasks cannot exit without having gone through wake_up_new_task() ->
	 * post_init_entity_util_avg() which will have added things to the
	 * cfs_rq, so we can remove unconditionally.
	 *
	 * Similarly for groups, they will have passed through
	 * post_init_entity_util_avg() before unregister_sched_fair_group()
	 * calls this.
	 */

	sync_entity_load_avg(se);

	raw_spin_lock_irqsave(&cfs_rq->removed.lock, flags);
	++cfs_rq->removed.nr;
	cfs_rq->removed.util_avg	+= se->avg.util_avg;
	cfs_rq->removed.load_avg	+= se->avg.load_avg;
	cfs_rq->removed.runnable_sum	+= se->avg.load_sum; /* == runnable_sum */
	raw_spin_unlock_irqrestore(&cfs_rq->removed.lock, flags);
}

static inline unsigned long cfs_rq_runnable_load_avg(struct cfs_rq *cfs_rq)
{
	return cfs_rq->avg.runnable_load_avg;
}

static inline unsigned long cfs_rq_load_avg(struct cfs_rq *cfs_rq)
{
	return cfs_rq->avg.load_avg;
}

static inline unsigned long _task_util_est(struct task_struct *p)
{
	struct util_est ue = READ_ONCE(p->se.avg.util_est);

	return max(ue.ewma, ue.enqueued);
}

unsigned long task_util_est(struct task_struct *p)
{
#ifdef CONFIG_SCHED_WALT
	return p->ravg.demand_scaled;
#endif
	return max(task_util(p), _task_util_est(p));
}

#ifdef CONFIG_UCLAMP_TASK
static inline unsigned long uclamp_task_util(struct task_struct *p)
{
	return clamp(task_util_est(p),
		     uclamp_eff_value(p, UCLAMP_MIN),
		     uclamp_eff_value(p, UCLAMP_MAX));
}
#else
static inline unsigned long uclamp_task_util(struct task_struct *p)
{
	return task_util_est(p);
}
#endif

static inline void util_est_enqueue(struct cfs_rq *cfs_rq,
				    struct task_struct *p)
{
	unsigned int enqueued;

	if (!sched_feat(UTIL_EST))
		return;

	/* Update root cfs_rq's estimated utilization */
	enqueued  = cfs_rq->avg.util_est.enqueued;
	enqueued += (_task_util_est(p) | UTIL_AVG_UNCHANGED);
	WRITE_ONCE(cfs_rq->avg.util_est.enqueued, enqueued);

	/* Update plots for Task and CPU estimated utilization */
	trace_sched_util_est_task(p, &p->se.avg);
	trace_sched_util_est_cpu(cpu_of(rq_of(cfs_rq)), cfs_rq);
}

/*
 * Check if a (signed) value is within a specified (unsigned) margin,
 * based on the observation that:
 *
 *     abs(x) < y := (unsigned)(x + y - 1) < (2 * y - 1)
 *
 * NOTE: this only works when value + maring < INT_MAX.
 */
static inline bool within_margin(int value, int margin)
{
	return ((unsigned int)(value + margin - 1) < (2 * margin - 1));
}

static void
util_est_dequeue(struct cfs_rq *cfs_rq, struct task_struct *p, bool task_sleep)
{
	long last_ewma_diff;
	struct util_est ue;
	int cpu;

	if (!sched_feat(UTIL_EST))
		return;

	/* Update root cfs_rq's estimated utilization */
	ue.enqueued  = cfs_rq->avg.util_est.enqueued;
	ue.enqueued -= min_t(unsigned int, ue.enqueued,
			     (_task_util_est(p) | UTIL_AVG_UNCHANGED));
	WRITE_ONCE(cfs_rq->avg.util_est.enqueued, ue.enqueued);

	/* Update plots for CPU's estimated utilization */
	trace_sched_util_est_cpu(cpu_of(rq_of(cfs_rq)), cfs_rq);

	/*
	 * Skip update of task's estimated utilization when the task has not
	 * yet completed an activation, e.g. being migrated.
	 */
	if (!task_sleep)
		return;

	/*
	 * If the PELT values haven't changed since enqueue time,
	 * skip the util_est update.
	 */
	ue = p->se.avg.util_est;
	if (ue.enqueued & UTIL_AVG_UNCHANGED)
		return;

	/*
	 * Reset EWMA on utilization increases, the moving average is used only
	 * to smooth utilization decreases.
	 */
	ue.enqueued = (task_util(p) | UTIL_AVG_UNCHANGED);
	if (sched_feat(UTIL_EST_FASTUP)) {
		if (ue.ewma < ue.enqueued) {
			ue.ewma = ue.enqueued;
			goto done;
		}
	}

	/*
	 * Skip update of task's estimated utilization when its EWMA is
	 * already ~1% close to its last activation value.
	 */
	last_ewma_diff = ue.enqueued - ue.ewma;
	if (within_margin(last_ewma_diff, (SCHED_CAPACITY_SCALE / 100)))
		return;

	/*
	 * To avoid overestimation of actual task utilization, skip updates if
	 * we cannot grant there is idle time in this CPU.
	 */
	cpu = cpu_of(rq_of(cfs_rq));
	if (task_util(p) > capacity_orig_of(cpu))
		return;

	/*
	 * Update Task's estimated utilization
	 *
	 * When *p completes an activation we can consolidate another sample
	 * of the task size. This is done by storing the current PELT value
	 * as ue.enqueued and by using this value to update the Exponential
	 * Weighted Moving Average (EWMA):
	 *
	 *  ewma(t) = w *  task_util(p) + (1-w) * ewma(t-1)
	 *          = w *  task_util(p) +         ewma(t-1)  - w * ewma(t-1)
	 *          = w * (task_util(p) -         ewma(t-1)) +     ewma(t-1)
	 *          = w * (      last_ewma_diff            ) +     ewma(t-1)
	 *          = w * (last_ewma_diff  +  ewma(t-1) / w)
	 *
	 * Where 'w' is the weight of new samples, which is configured to be
	 * 0.25, thus making w=1/4 ( >>= UTIL_EST_WEIGHT_SHIFT)
	 */
	ue.ewma <<= UTIL_EST_WEIGHT_SHIFT;
	ue.ewma  += last_ewma_diff;
	ue.ewma >>= UTIL_EST_WEIGHT_SHIFT;
done:
	WRITE_ONCE(p->se.avg.util_est, ue);

	/* Update plots for Task's estimated utilization */
	trace_sched_util_est_task(p, &p->se.avg);
}
#ifdef CONFIG_SCHED_WALT
/*
 * Check whether cpu is in the fastest set of cpu's that p should run on.
 * If p is prefer_high_cap, prefer that p runs on a faster cpu; otherwise,
 * allow p to run on any cpu.
 */
static inline bool cpu_is_in_target_set(struct task_struct *p, int cpu)
{
	struct root_domain *rd = cpu_rq(cpu)->rd;
	int first_cpu, next_usable_cpu;
#ifdef CONFIG_SCHED_TUNE
	if (schedtune_prefer_high_cap(p)) {
#else
	if (uclamp_boosted(p)) {
#endif
		first_cpu = rd->mid_cap_orig_cpu != -1 ? rd->mid_cap_orig_cpu :
			    rd->max_cap_orig_cpu;

	} else {
		first_cpu = rd->min_cap_orig_cpu;
	}

	next_usable_cpu = cpumask_next(first_cpu - 1, &p->cpus_allowed);
	return cpu >= next_usable_cpu || next_usable_cpu >= nr_cpu_ids;
}
#endif

static inline bool
bias_to_this_cpu(struct task_struct *p, int cpu, int start_cpu)
{
	bool base_test = cpumask_test_cpu(cpu, &p->cpus_allowed) &&
#ifdef CONFIG_SCHED_WALT
			cpu_active(cpu) && cpu_is_in_target_set(p, cpu);
#else
						   cpu_active(cpu);
#endif
	bool start_cap_test = (capacity_orig_of(cpu) >=
					capacity_orig_of(start_cpu));

	return base_test && start_cap_test;
}

#ifdef CONFIG_SCHED_WALT
static inline bool task_fits_capacity(struct task_struct *p,
					long capacity,
					int cpu)
{
	unsigned int margin;

	/*
	 * Derive upmigration/downmigrate margin wrt the src/dest
	 * CPU.
	 */
	if (capacity_orig_of(task_cpu(p)) > capacity_orig_of(cpu))
		margin = schedtune_task_boost(p) > 0 &&
			  !schedtune_prefer_high_cap(p) &&
			   p->prio <= DEFAULT_PRIO ?
			sched_capacity_margin_down_boosted[task_cpu(p)] :
			sched_capacity_margin_down[task_cpu(p)];
	else
		margin = schedtune_task_boost(p) > 0 &&
			  !schedtune_prefer_high_cap(p) &&
			   p->prio <= DEFAULT_PRIO ?
			sched_capacity_margin_up_boosted[task_cpu(p)] :
			sched_capacity_margin_up[task_cpu(p)];

	return capacity * 1024 > uclamp_task(p) * margin;
}
#else
static inline int task_fits_capacity(struct task_struct *p, long capacity)
{
	return fits_capacity(uclamp_task_util(p), capacity);
}
#endif

static inline bool task_fits_max(struct task_struct *p, int cpu)
{
#ifdef CONFIG_SCHED_WALT
	unsigned long capacity = capacity_orig_of(cpu);
	unsigned long max_capacity = cpu_rq(cpu)->rd->max_cpu_capacity.val;
	unsigned long task_boost = per_task_boost(p);

	if (capacity == max_capacity)
		return true;

	if (is_min_capacity_cpu(cpu)) {
		if (task_boost_policy(p) == SCHED_BOOST_ON_BIG ||
			task_boost > 0 ||
			walt_should_kick_upmigrate(p, cpu))
			return false;
	} else { /* mid cap cpu */
		if (task_boost > TASK_BOOST_ON_MID)
			return false;
	}

	return task_fits_capacity(p, capacity, cpu);
#else
	return false;
#endif
}

static inline bool task_demand_fits(struct task_struct *p, int cpu)
{
#ifdef CONFIG_SCHED_WALT
	unsigned long capacity = capacity_orig_of(cpu);
	unsigned long max_capacity = cpu_rq(cpu)->rd->max_cpu_capacity.val;

	if (capacity == max_capacity)
		return true;

	return task_fits_capacity(p, capacity, cpu);
#else
	return false;
#endif
}

struct find_best_target_env {
	int placement_boost;
	int need_idle;
	int fastpath;
	int start_cpu;
	int skip_cpu;
	bool is_rtg;
	bool boosted;
	bool strict_max;
};

static inline void adjust_cpus_for_packing(struct task_struct *p,
			int *target_cpu, int *best_idle_cpu,
			int shallowest_idle_cstate,
			struct find_best_target_env *fbt_env,
			bool boosted)
{
	unsigned long tutil, estimated_capacity;

	if (*best_idle_cpu == -1 || *target_cpu == -1)
		return;

	if (prefer_spread_on_idle(*best_idle_cpu, false))
		fbt_env->need_idle |= 2;

#ifdef CONFIG_SCHED_WALT
	if (task_rtg_high_prio(p) && walt_nr_rtg_high_prio(*target_cpu) > 0) {
		*target_cpu = -1;
		return;
	}
#endif

	if (fbt_env->need_idle || task_placement_boost_enabled(p) || boosted ||
		shallowest_idle_cstate <= 0) {
		*target_cpu = -1;
		return;
	}

	if (task_in_cum_window_demand(cpu_rq(*target_cpu), p))
		tutil = 0;
	else
		tutil = task_util(p);

	estimated_capacity = cpu_util_cum(*target_cpu, tutil);
	estimated_capacity = add_capacity_margin(estimated_capacity,
							*target_cpu);

	/*
	 * If there is only one active CPU and it is already above its current
	 * capacity, avoid placing additional task on the CPU.
	 */
	if (estimated_capacity > capacity_curr_of(*target_cpu)) {
		*target_cpu = -1;
		return;
	}

	if (fbt_env->is_rtg)
		*best_idle_cpu = -1;
}

static inline void update_misfit_status(struct task_struct *p, struct rq *rq)
{
	if (!static_branch_unlikely(&sched_asym_cpucapacity))
		return;

	if (!p || p->nr_cpus_allowed == 1) {
		rq->misfit_task_load = 0;
		return;
	}
#ifdef CONFIG_SCHED_WALT
	if (task_fits_max(p, cpu_of(rq))) {
#else
        if (task_fits_capacity(p, cpu_of(rq))) {
#endif
		rq->misfit_task_load = 0;
		return;
	}

	/*
	 * Make sure that misfit_task_load will not be null even if
	 * task_h_load() returns 0.
	 */
	rq->misfit_task_load = max_t(unsigned long, task_h_load(p), 1);
}

#else /* CONFIG_SMP */

#define UPDATE_TG	0x0
#define SKIP_AGE_LOAD	0x0
#define DO_ATTACH	0x0

static inline void update_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se, int not_used1)
{
	cfs_rq_util_change(cfs_rq, 0);
}

static inline void remove_entity_load_avg(struct sched_entity *se) {}

static inline void
attach_entity_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se, int flags) {}
static inline void
detach_entity_load_avg(struct cfs_rq *cfs_rq, struct sched_entity *se) {}

static inline int idle_balance(struct rq *rq, struct rq_flags *rf)
{
	return 0;
}

static inline void
util_est_enqueue(struct cfs_rq *cfs_rq, struct task_struct *p) {}

static inline void
util_est_dequeue(struct cfs_rq *cfs_rq, struct task_struct *p,
		 bool task_sleep) {}
static inline void update_misfit_status(struct task_struct *p, struct rq *rq) {}

#endif /* CONFIG_SMP */

static void
place_entity(struct cfs_rq *cfs_rq, struct sched_entity *se, int flags)
{
	u64 vslice, vruntime = avg_vruntime(cfs_rq);
	s64 lag = 0;

	se->slice = sysctl_sched_base_slice;
	vslice = calc_delta_fair(se->slice, se);

	/*
	 * Due to how V is constructed as the weighted average of entities,
	 * adding tasks with positive lag, or removing tasks with negative lag
	 * will move 'time' backwards, this can screw around with the lag of
	 * other tasks.
	 *
	 * EEVDF: placement strategy #1 / #2
	 */
	if (sched_feat(PLACE_LAG) && cfs_rq->nr_running) {
		struct sched_entity *curr = cfs_rq->curr;
		unsigned long load;

		lag = se->vlag;

		/*
		 * If we want to place a task and preserve lag, we have to
		 * consider the effect of the new entity on the weighted
		 * average and compensate for this, otherwise lag can quickly
		 * evaporate.
		 *
		 * Lag is defined as:
		 *
		 *   lag_i = S - s_i = w_i * (V - v_i)
		 *
		 * To avoid the 'w_i' term all over the place, we only track
		 * the virtual lag:
		 *
		 *   vl_i = V - v_i <=> v_i = V - vl_i
		 *
		 * And we take V to be the weighted average of all v:
		 *
		 *   V = (\Sum w_j*v_j) / W
		 *
		 * Where W is: \Sum w_j
		 *
		 * Then, the weighted average after adding an entity with lag
		 * vl_i is given by:
		 *
		 *   V' = (\Sum w_j*v_j + w_i*v_i) / (W + w_i)
		 *      = (W*V + w_i*(V - vl_i)) / (W + w_i)
		 *      = (W*V + w_i*V - w_i*vl_i) / (W + w_i)
		 *      = (V*(W + w_i) - w_i*l) / (W + w_i)
		 *      = V - w_i*vl_i / (W + w_i)
		 *
		 * And the actual lag after adding an entity with vl_i is:
		 *
		 *   vl'_i = V' - v_i
		 *         = V - w_i*vl_i / (W + w_i) - (V - vl_i)
		 *         = vl_i - w_i*vl_i / (W + w_i)
		 *
		 * Which is strictly less than vl_i. So in order to preserve lag
		 * we should inflate the lag before placement such that the
		 * effective lag after placement comes out right.
		 *
		 * As such, invert the above relation for vl'_i to get the vl_i
		 * we need to use such that the lag after placement is the lag
		 * we computed before dequeue.
		 *
		 *   vl'_i = vl_i - w_i*vl_i / (W + w_i)
		 *         = ((W + w_i)*vl_i - w_i*vl_i) / (W + w_i)
		 *
		 *   (W + w_i)*vl'_i = (W + w_i)*vl_i - w_i*vl_i
		 *                   = W*vl_i
		 *
		 *   vl_i = (W + w_i)*vl'_i / W
		 */
		load = cfs_rq->avg_load;
		if (curr && curr->on_rq)
			load += scale_load_down(curr->load.weight);

		lag *= load + scale_load_down(se->load.weight);
		if (WARN_ON_ONCE(!load))
			load = 1;
		lag = div_s64(lag, load);
	}

	se->vruntime = vruntime - lag;

	/*
	 * When joining the competition; the exisiting tasks will be,
	 * on average, halfway through their slice, as such start tasks
	 * off with half a slice to ease into the competition.
	 */
	if (sched_feat(PLACE_DEADLINE_INITIAL) && (flags & ENQUEUE_INITIAL))
		vslice /= 2;

	/*
	 * EEVDF: vd_i = ve_i + r_i/w_i
	 */
	se->deadline = se->vruntime + vslice;
}

static void check_enqueue_throttle(struct cfs_rq *cfs_rq);

static inline void check_schedstat_required(void)
{
#ifdef CONFIG_SCHEDSTATS
	if (schedstat_enabled())
		return;

	/* Force schedstat enabled if a dependent tracepoint is active */
	if (trace_sched_stat_wait_enabled()    ||
			trace_sched_stat_sleep_enabled()   ||
			trace_sched_stat_iowait_enabled()  ||
			trace_sched_stat_blocked_enabled() ||
			trace_sched_stat_runtime_enabled())  {
		printk_deferred_once("Scheduler tracepoints stat_sleep, stat_iowait, "
			     "stat_blocked and stat_runtime require the "
			     "kernel parameter schedstats=enable or "
			     "kernel.sched_schedstats=1\n");
	}
#endif
}


static void
enqueue_entity(struct cfs_rq *cfs_rq, struct sched_entity *se, int flags)
{
	bool curr = cfs_rq->curr == se;

	/*
	 * If we're the current task, we must renormalise before calling
	 * update_curr().
	 */
	if (curr)
		place_entity(cfs_rq, se, flags);

	update_curr(cfs_rq);

	/*
	 * When enqueuing a sched_entity, we must:
	 *   - Update loads to have both entity and cfs_rq synced with now.
	 *   - Add its load to cfs_rq->runnable_avg
	 *   - For group_entity, update its weight to reflect the new share of
	 *     its group cfs_rq
	 *   - Add its new weight to cfs_rq->load.weight
	 */
	update_load_avg(cfs_rq, se, UPDATE_TG | DO_ATTACH);
	enqueue_runnable_load_avg(cfs_rq, se);
	/*
	 * XXX update_load_avg() above will have attached us to the pelt sum;
	 * but update_cfs_group() here will re-adjust the weight and have to
	 * undo/redo all that. Seems wasteful.
	 */
	update_cfs_group(se);

	/*
	 * XXX now that the entity has been re-weighted, and it's lag adjusted,
	 * we can place the entity.
	 */
	if (!curr)
		place_entity(cfs_rq, se, flags);

	account_entity_enqueue(cfs_rq, se);

	if (flags & ENQUEUE_WAKEUP)
		place_entity(cfs_rq, se, 0);
	/* Entity has migrated, no longer consider this task hot */
	if (flags & ENQUEUE_MIGRATED)
		se->exec_start = 0;

	check_schedstat_required();
	update_stats_enqueue(cfs_rq, se, flags);
	if (!curr)
		__enqueue_entity(cfs_rq, se);
	se->on_rq = 1;

	if (cfs_rq->nr_running == 1) {
		list_add_leaf_cfs_rq(cfs_rq);
		check_enqueue_throttle(cfs_rq);
	}
}

static void __clear_buddies_next(struct sched_entity *se)
{
	for_each_sched_entity(se) {
		struct cfs_rq *cfs_rq = cfs_rq_of(se);
		if (cfs_rq->next != se)
			break;

		cfs_rq->next = NULL;
	}
}

static void clear_buddies(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	if (cfs_rq->next == se)
		__clear_buddies_next(se);
}

static __always_inline void return_cfs_rq_runtime(struct cfs_rq *cfs_rq);

static void
dequeue_entity(struct cfs_rq *cfs_rq, struct sched_entity *se, int flags)
{
	/*
	 * Update run-time statistics of the 'current'.
	 */
	update_curr(cfs_rq);

	/*
	 * When dequeuing a sched_entity, we must:
	 *   - Update loads to have both entity and cfs_rq synced with now.
	 *   - Substract its load from the cfs_rq->runnable_avg.
	 *   - Substract its previous weight from cfs_rq->load.weight.
	 *   - For group entity, update its weight to reflect the new share
	 *     of its group cfs_rq.
	 */
	update_load_avg(cfs_rq, se, UPDATE_TG);
	dequeue_runnable_load_avg(cfs_rq, se);

	update_stats_dequeue(cfs_rq, se, flags);

	clear_buddies(cfs_rq, se);

	update_entity_lag(cfs_rq, se);
	if (se != cfs_rq->curr)
		__dequeue_entity(cfs_rq, se);
	se->on_rq = 0;
	account_entity_dequeue(cfs_rq, se);

	/* return excess runtime on last dequeue */
	return_cfs_rq_runtime(cfs_rq);

	update_cfs_group(se);

	/*
	 * Now advance min_vruntime if @se was the entity holding it back,
	 * except when: DEQUEUE_SAVE && !DEQUEUE_MOVE, in this case we'll be
	 * put back on, and if we advance min_vruntime, we'll be placed back
	 * further than we started -- ie. we'll be penalized.
	 */
	if ((flags & (DEQUEUE_SAVE | DEQUEUE_MOVE)) != DEQUEUE_SAVE)
		update_min_vruntime(cfs_rq);
}

static void
set_next_entity(struct cfs_rq *cfs_rq, struct sched_entity *se)
{
	clear_buddies(cfs_rq, se);

	/* 'current' is not kept within the tree. */
	if (se->on_rq) {
		/*
		 * Any task has to be enqueued before it get to execute on
		 * a CPU. So account for the time it spent waiting on the
		 * runqueue.
		 */
		update_stats_wait_end(cfs_rq, se);
		__dequeue_entity(cfs_rq, se);
		update_load_avg(cfs_rq, se, UPDATE_TG);
		/*
		 * HACK, stash a copy of deadline at the point of pick in vlag,
		 * which isn't used until dequeue.
		 */
		se->vlag = se->deadline;
	}

	update_stats_curr_start(cfs_rq, se);
	cfs_rq->curr = se;

	/*
	 * Track our maximum slice length, if the CPU's load is at
	 * least twice that of our own weight (i.e. dont track it
	 * when there are only lesser-weight tasks around):
	 */
	if (schedstat_enabled() && rq_of(cfs_rq)->load.weight >= 2*se->load.weight) {
		schedstat_set(se->statistics.slice_max,
			max((u64)schedstat_val(se->statistics.slice_max),
			    se->sum_exec_runtime - se->prev_sum_exec_runtime));
	}

	se->prev_sum_exec_runtime = se->sum_exec_runtime;
}

/*
 * Pick the next process, keeping these things in mind, in this order:
 * 1) keep things fair between processes/task groups
 * 2) pick the "next" process, since someone really wants that to run
 * 3) pick the "last" process, for cache locality
 * 4) do not run the "skip" process, if something else is available
 */
static struct sched_entity *
pick_next_entity(struct cfs_rq *cfs_rq, struct sched_entity *curr)
{
	/*
	 * Enabling NEXT_BUDDY will affect latency but not fairness.
	 */
	if (sched_feat(NEXT_BUDDY) &&
	    cfs_rq->next && entity_eligible(cfs_rq, cfs_rq->next))
		return cfs_rq->next;

	return pick_eevdf(cfs_rq);
}

static bool check_cfs_rq_runtime(struct cfs_rq *cfs_rq);

static void put_prev_entity(struct cfs_rq *cfs_rq, struct sched_entity *prev)
{
	/*
	 * If still on the runqueue then deactivate_task()
	 * was not called and update_curr() has to be done:
	 */
	if (prev->on_rq)
		update_curr(cfs_rq);

	/* throttle cfs_rqs exceeding runtime */
	check_cfs_rq_runtime(cfs_rq);

	if (prev->on_rq) {
		update_stats_wait_start(cfs_rq, prev);
		/* Put 'current' back into the tree. */
		__enqueue_entity(cfs_rq, prev);
		/* in !on_rq case, update occurred at dequeue */
		update_load_avg(cfs_rq, prev, 0);
	}
	cfs_rq->curr = NULL;
}

static void
entity_tick(struct cfs_rq *cfs_rq, struct sched_entity *curr, int queued)
{
	/*
	 * Update run-time statistics of the 'current'.
	 */
	update_curr(cfs_rq);

	/*
	 * Ensure that runnable average is periodically updated.
	 */
	update_load_avg(cfs_rq, curr, UPDATE_TG);
	update_cfs_group(curr);

#ifdef CONFIG_SCHED_HRTICK
	/*
	 * queued ticks are scheduled to match the slice, so don't bother
	 * validating it and just reschedule.
	 */
	if (queued) {
		resched_curr(rq_of(cfs_rq));
		return;
	}
	/*
	 * don't let the period tick interfere with the hrtick preemption
	 */
	if (!sched_feat(DOUBLE_TICK) &&
			hrtimer_active(&rq_of(cfs_rq)->hrtick_timer))
		return;
#endif
}


/**************************************************
 * CFS bandwidth control machinery
 */

#ifdef CONFIG_CFS_BANDWIDTH

#ifdef CONFIG_JUMP_LABEL
static struct static_key __cfs_bandwidth_used;

static inline bool cfs_bandwidth_used(void)
{
	return static_key_false(&__cfs_bandwidth_used);
}

void cfs_bandwidth_usage_inc(void)
{
	static_key_slow_inc_cpuslocked(&__cfs_bandwidth_used);
}

void cfs_bandwidth_usage_dec(void)
{
	static_key_slow_dec_cpuslocked(&__cfs_bandwidth_used);
}
#else /* CONFIG_JUMP_LABEL */
static bool cfs_bandwidth_used(void)
{
	return true;
}

void cfs_bandwidth_usage_inc(void) {}
void cfs_bandwidth_usage_dec(void) {}
#endif /* CONFIG_JUMP_LABEL */

/*
 * default period for cfs group bandwidth.
 * default: 0.1s, units: nanoseconds
 */
static inline u64 default_cfs_period(void)
{
	return 100000000ULL;
}

static inline u64 sched_cfs_bandwidth_slice(void)
{
	return (u64)sysctl_sched_cfs_bandwidth_slice * NSEC_PER_USEC;
}

/*
 * Replenish runtime according to assigned quota. We use sched_clock_cpu
 * directly instead of rq->clock to avoid adding additional synchronization
 * around rq->lock.
 *
 * requires cfs_b->lock
 */
void __refill_cfs_bandwidth_runtime(struct cfs_bandwidth *cfs_b)
{
	if (cfs_b->quota != RUNTIME_INF)
		cfs_b->runtime = cfs_b->quota;
}

static inline struct cfs_bandwidth *tg_cfs_bandwidth(struct task_group *tg)
{
	return &tg->cfs_bandwidth;
}

/* rq->task_clock normalized against any time this cfs_rq has spent throttled */
static inline u64 cfs_rq_clock_task(struct cfs_rq *cfs_rq)
{
	if (unlikely(cfs_rq->throttle_count))
		return cfs_rq->throttled_clock_pelt - cfs_rq->throttled_clock_pelt_time;

	return rq_clock_task(rq_of(cfs_rq)) - cfs_rq->throttled_clock_pelt_time;
}

/* returns 0 on failure to allocate runtime */
static int assign_cfs_rq_runtime(struct cfs_rq *cfs_rq)
{
	struct task_group *tg = cfs_rq->tg;
	struct cfs_bandwidth *cfs_b = tg_cfs_bandwidth(tg);
	u64 amount = 0, min_amount;

	/* note: this is a positive sum as runtime_remaining <= 0 */
	min_amount = sched_cfs_bandwidth_slice() - cfs_rq->runtime_remaining;

	raw_spin_lock(&cfs_b->lock);
	if (cfs_b->quota == RUNTIME_INF)
		amount = min_amount;
	else {
		start_cfs_bandwidth(cfs_b);

		if (cfs_b->runtime > 0) {
			amount = min(cfs_b->runtime, min_amount);
			cfs_b->runtime -= amount;
			cfs_b->idle = 0;
		}
	}
	raw_spin_unlock(&cfs_b->lock);

	cfs_rq->runtime_remaining += amount;

	return cfs_rq->runtime_remaining > 0;
}

static void __account_cfs_rq_runtime(struct cfs_rq *cfs_rq, u64 delta_exec)
{
	/* dock delta_exec before expiring quota (as it could span periods) */
	cfs_rq->runtime_remaining -= delta_exec;

	if (likely(cfs_rq->runtime_remaining > 0))
		return;

	if (cfs_rq->throttled)
		return;
	/*
	 * if we're unable to extend our runtime we resched so that the active
	 * hierarchy can be throttled
	 */
	if (!assign_cfs_rq_runtime(cfs_rq) && likely(cfs_rq->curr))
		resched_curr(rq_of(cfs_rq));
}

static __always_inline
void account_cfs_rq_runtime(struct cfs_rq *cfs_rq, u64 delta_exec)
{
	if (!cfs_bandwidth_used() || !cfs_rq->runtime_enabled)
		return;

	__account_cfs_rq_runtime(cfs_rq, delta_exec);
}

static inline int cfs_rq_throttled(struct cfs_rq *cfs_rq)
{
	return cfs_bandwidth_used() && cfs_rq->throttled;
}

/* check whether cfs_rq, or any parent, is throttled */
static inline int throttled_hierarchy(struct cfs_rq *cfs_rq)
{
	return cfs_bandwidth_used() && cfs_rq->throttle_count;
}

/*
 * Ensure that neither of the group entities corresponding to src_cpu or
 * dest_cpu are members of a throttled hierarchy when performing group
 * load-balance operations.
 */
static inline int throttled_lb_pair(struct task_group *tg,
				    int src_cpu, int dest_cpu)
{
	struct cfs_rq *src_cfs_rq, *dest_cfs_rq;

	src_cfs_rq = tg->cfs_rq[src_cpu];
	dest_cfs_rq = tg->cfs_rq[dest_cpu];

	return throttled_hierarchy(src_cfs_rq) ||
	       throttled_hierarchy(dest_cfs_rq);
}

static int tg_unthrottle_up(struct task_group *tg, void *data)
{
	struct rq *rq = data;
	struct cfs_rq *cfs_rq = tg->cfs_rq[cpu_of(rq)];

	cfs_rq->throttle_count--;
	if (!cfs_rq->throttle_count) {
		/* adjust cfs_rq_clock_task() */
		cfs_rq->throttled_clock_pelt_time += rq_clock_pelt(rq) -
					     cfs_rq->throttled_clock_pelt;

		/* Add cfs_rq with already running entity in the list */
		if (cfs_rq->nr_running >= 1)
			list_add_leaf_cfs_rq(cfs_rq);
	}

	return 0;
}

static int tg_throttle_down(struct task_group *tg, void *data)
{
	struct rq *rq = data;
	struct cfs_rq *cfs_rq = tg->cfs_rq[cpu_of(rq)];

	/* group is entering throttled state, stop time */
	if (!cfs_rq->throttle_count) {
		cfs_rq->throttled_clock_pelt = rq_clock_pelt(rq);
		list_del_leaf_cfs_rq(cfs_rq);
	}
	cfs_rq->throttle_count++;

	return 0;
}

static void throttle_cfs_rq(struct cfs_rq *cfs_rq)
{
	struct rq *rq = rq_of(cfs_rq);
	struct cfs_bandwidth *cfs_b = tg_cfs_bandwidth(cfs_rq->tg);
	struct sched_entity *se;
	long task_delta, idle_task_delta, dequeue = 1;
	bool empty;

	se = cfs_rq->tg->se[cpu_of(rq_of(cfs_rq))];

	/* freeze hierarchy runnable averages while throttled */
	rcu_read_lock();
	walk_tg_tree_from(cfs_rq->tg, tg_throttle_down, tg_nop, (void *)rq);
	rcu_read_unlock();

	task_delta = cfs_rq->h_nr_running;
	idle_task_delta = cfs_rq->idle_h_nr_running;
	for_each_sched_entity(se) {
		struct cfs_rq *qcfs_rq = cfs_rq_of(se);
		/* throttled entity or throttle-on-deactivate */
		if (!se->on_rq)
			break;

		if (dequeue)
			dequeue_entity(qcfs_rq, se, DEQUEUE_SLEEP);
		qcfs_rq->h_nr_running -= task_delta;
                qcfs_rq->idle_h_nr_running -= idle_task_delta;
		walt_dec_throttled_cfs_rq_stats(&qcfs_rq->walt_stats, cfs_rq);

		if (qcfs_rq->load.weight)
			dequeue = 0;
	}

	if (!se) {
		sub_nr_running(rq, task_delta);
		walt_dec_throttled_cfs_rq_stats(&rq->walt_stats, cfs_rq);
	}

	cfs_rq->throttled = 1;
	cfs_rq->throttled_clock = rq_clock(rq);
	raw_spin_lock(&cfs_b->lock);
	empty = list_empty(&cfs_b->throttled_cfs_rq);

	/*
	 * Add to the _head_ of the list, so that an already-started
	 * distribute_cfs_runtime will not see us. If disribute_cfs_runtime is
	 * not running add to the tail so that later runqueues don't get starved.
	 */
	if (cfs_b->distribute_running)
		list_add_rcu(&cfs_rq->throttled_list, &cfs_b->throttled_cfs_rq);
	else
		list_add_tail_rcu(&cfs_rq->throttled_list, &cfs_b->throttled_cfs_rq);

	/*
	 * If we're the first throttled task, make sure the bandwidth
	 * timer is running.
	 */
	if (empty)
		start_cfs_bandwidth(cfs_b);

	raw_spin_unlock(&cfs_b->lock);
}

void unthrottle_cfs_rq(struct cfs_rq *cfs_rq)
{
	struct rq *rq = rq_of(cfs_rq);
	struct cfs_bandwidth *cfs_b = tg_cfs_bandwidth(cfs_rq->tg);
	struct sched_entity *se;
	int enqueue = 1;
        long task_delta, idle_task_delta;
	struct cfs_rq *tcfs_rq __maybe_unused = cfs_rq;

	se = cfs_rq->tg->se[cpu_of(rq)];

	cfs_rq->throttled = 0;

	update_rq_clock(rq);

	raw_spin_lock(&cfs_b->lock);
	cfs_b->throttled_time += rq_clock(rq) - cfs_rq->throttled_clock;
	list_del_rcu(&cfs_rq->throttled_list);
	raw_spin_unlock(&cfs_b->lock);

	/* update hierarchical throttle state */
	walk_tg_tree_from(cfs_rq->tg, tg_nop, tg_unthrottle_up, (void *)rq);

	if (!cfs_rq->load.weight)
		return;

	task_delta = cfs_rq->h_nr_running;
	idle_task_delta = cfs_rq->idle_h_nr_running;
	for_each_sched_entity(se) {
		if (se->on_rq)
			enqueue = 0;

		cfs_rq = cfs_rq_of(se);
		if (enqueue)
			enqueue_entity(cfs_rq, se, ENQUEUE_WAKEUP);
		cfs_rq->h_nr_running += task_delta;
                cfs_rq->idle_h_nr_running += idle_task_delta;
		walt_inc_throttled_cfs_rq_stats(&cfs_rq->walt_stats, tcfs_rq);

		if (cfs_rq_throttled(cfs_rq))
			break;
	}

	assert_list_leaf_cfs_rq(rq);

	if (!se) {
		add_nr_running(rq, task_delta);
		walt_inc_throttled_cfs_rq_stats(&rq->walt_stats, tcfs_rq);
	}

	/* Determine whether we need to wake up potentially idle CPU: */
	if (rq->curr == rq->idle && rq->cfs.nr_running)
		resched_curr(rq);
}

static u64 distribute_cfs_runtime(struct cfs_bandwidth *cfs_b, u64 remaining)
{
	struct cfs_rq *cfs_rq;
	u64 runtime;
	u64 starting_runtime = remaining;

	rcu_read_lock();
	list_for_each_entry_rcu(cfs_rq, &cfs_b->throttled_cfs_rq,
				throttled_list) {
		struct rq *rq = rq_of(cfs_rq);
		struct rq_flags rf;

		rq_lock(rq, &rf);
		if (!cfs_rq_throttled(cfs_rq))
			goto next;

		/* By the above check, this should never be true */
		SCHED_WARN_ON(cfs_rq->runtime_remaining > 0);

		runtime = -cfs_rq->runtime_remaining + 1;
		if (runtime > remaining)
			runtime = remaining;
		remaining -= runtime;

		cfs_rq->runtime_remaining += runtime;

		/* we check whether we're throttled above */
		if (cfs_rq->runtime_remaining > 0)
			unthrottle_cfs_rq(cfs_rq);

next:
		rq_unlock(rq, &rf);

		if (!remaining)
			break;
	}
	rcu_read_unlock();

	return starting_runtime - remaining;
}

/*
 * Responsible for refilling a task_group's bandwidth and unthrottling its
 * cfs_rqs as appropriate. If there has been no activity within the last
 * period the timer is deactivated until scheduling resumes; cfs_b->idle is
 * used to track this state.
 */
static int do_sched_cfs_period_timer(struct cfs_bandwidth *cfs_b, int overrun)
{
	u64 runtime;
	int throttled;

	/* no need to continue the timer with no bandwidth constraint */
	if (cfs_b->quota == RUNTIME_INF)
		goto out_deactivate;

	throttled = !list_empty(&cfs_b->throttled_cfs_rq);
	cfs_b->nr_periods += overrun;

	/*
	 * idle depends on !throttled (for the case of a large deficit), and if
	 * we're going inactive then everything else can be deferred
	 */
	if (cfs_b->idle && !throttled)
		goto out_deactivate;

	__refill_cfs_bandwidth_runtime(cfs_b);

	if (!throttled) {
		/* mark as potentially idle for the upcoming period */
		cfs_b->idle = 1;
		return 0;
	}

	/* account preceding periods in which throttling occurred */
	cfs_b->nr_throttled += overrun;

	/*
	 * This check is repeated as we are holding onto the new bandwidth while
	 * we unthrottle. This can potentially race with an unthrottled group
	 * trying to acquire new bandwidth from the global pool. This can result
	 * in us over-using our runtime if it is all used during this loop, but
	 * only by limited amounts in that extreme case.
	 */
	while (throttled && cfs_b->runtime > 0 && !cfs_b->distribute_running) {
		runtime = cfs_b->runtime;
		cfs_b->distribute_running = 1;
		raw_spin_unlock(&cfs_b->lock);
		/* we can't nest cfs_b->lock while distributing bandwidth */
		runtime = distribute_cfs_runtime(cfs_b, runtime);
		raw_spin_lock(&cfs_b->lock);

		cfs_b->distribute_running = 0;
		throttled = !list_empty(&cfs_b->throttled_cfs_rq);

		cfs_b->runtime -= min(runtime, cfs_b->runtime);
	}

	/*
	 * While we are ensured activity in the period following an
	 * unthrottle, this also covers the case in which the new bandwidth is
	 * insufficient to cover the existing bandwidth deficit.  (Forcing the
	 * timer to remain active while there are any throttled entities.)
	 */
	cfs_b->idle = 0;

	return 0;

out_deactivate:
	return 1;
}

/* a cfs_rq won't donate quota below this amount */
static const u64 min_cfs_rq_runtime = 1 * NSEC_PER_MSEC;
/* minimum remaining period time to redistribute slack quota */
static const u64 min_bandwidth_expiration = 2 * NSEC_PER_MSEC;
/* how long we wait to gather additional slack before distributing */
static const u64 cfs_bandwidth_slack_period = 5 * NSEC_PER_MSEC;

/*
 * Are we near the end of the current quota period?
 *
 * Requires cfs_b->lock for hrtimer_expires_remaining to be safe against the
 * hrtimer base being cleared by hrtimer_start. In the case of
 * migrate_hrtimers, base is never cleared, so we are fine.
 */
static int runtime_refresh_within(struct cfs_bandwidth *cfs_b, u64 min_expire)
{
	struct hrtimer *refresh_timer = &cfs_b->period_timer;
	u64 remaining;

	/* if the call-back is running a quota refresh is already occurring */
	if (hrtimer_callback_running(refresh_timer))
		return 1;

	/* is a quota refresh about to occur? */
	remaining = ktime_to_ns(hrtimer_expires_remaining(refresh_timer));
	if (remaining < min_expire)
		return 1;

	return 0;
}

static void start_cfs_slack_bandwidth(struct cfs_bandwidth *cfs_b)
{
	u64 min_left = cfs_bandwidth_slack_period + min_bandwidth_expiration;

	/* if there's a quota refresh soon don't bother with slack */
	if (runtime_refresh_within(cfs_b, min_left))
		return;

	hrtimer_start(&cfs_b->slack_timer,
			ns_to_ktime(cfs_bandwidth_slack_period),
			HRTIMER_MODE_REL);
}

/* we know any runtime found here is valid as update_curr() precedes return */
static void __return_cfs_rq_runtime(struct cfs_rq *cfs_rq)
{
	struct cfs_bandwidth *cfs_b = tg_cfs_bandwidth(cfs_rq->tg);
	s64 slack_runtime = cfs_rq->runtime_remaining - min_cfs_rq_runtime;

	if (slack_runtime <= 0)
		return;

	raw_spin_lock(&cfs_b->lock);
	if (cfs_b->quota != RUNTIME_INF) {
		cfs_b->runtime += slack_runtime;

		/* we are under rq->lock, defer unthrottling using a timer */
		if (cfs_b->runtime > sched_cfs_bandwidth_slice() &&
		    !list_empty(&cfs_b->throttled_cfs_rq))
			start_cfs_slack_bandwidth(cfs_b);
	}
	raw_spin_unlock(&cfs_b->lock);

	/* even if it's not valid for return we don't want to try again */
	cfs_rq->runtime_remaining -= slack_runtime;
}

static __always_inline void return_cfs_rq_runtime(struct cfs_rq *cfs_rq)
{
	if (!cfs_bandwidth_used())
		return;

	if (!cfs_rq->runtime_enabled || cfs_rq->nr_running)
		return;

	__return_cfs_rq_runtime(cfs_rq);
}

/*
 * This is done with a timer (instead of inline with bandwidth return) since
 * it's necessary to juggle rq->locks to unthrottle their respective cfs_rqs.
 */
static void do_sched_cfs_slack_timer(struct cfs_bandwidth *cfs_b)
{
	u64 runtime = 0, slice = sched_cfs_bandwidth_slice();

	/* confirm we're still not at a refresh boundary */
	raw_spin_lock(&cfs_b->lock);
	if (cfs_b->distribute_running) {
		raw_spin_unlock(&cfs_b->lock);
		return;
	}

	if (runtime_refresh_within(cfs_b, min_bandwidth_expiration)) {
		raw_spin_unlock(&cfs_b->lock);
		return;
	}

	if (cfs_b->quota != RUNTIME_INF && cfs_b->runtime > slice)
		runtime = cfs_b->runtime;

	if (runtime)
		cfs_b->distribute_running = 1;

	raw_spin_unlock(&cfs_b->lock);

	if (!runtime)
		return;

	runtime = distribute_cfs_runtime(cfs_b, runtime);

	raw_spin_lock(&cfs_b->lock);
	cfs_b->runtime -= min(runtime, cfs_b->runtime);
	cfs_b->distribute_running = 0;
	raw_spin_unlock(&cfs_b->lock);
}

/*
 * When a group wakes up we want to make sure that its quota is not already
 * expired/exceeded, otherwise it may be allowed to steal additional ticks of
 * runtime as update_curr() throttling can not not trigger until it's on-rq.
 */
static void check_enqueue_throttle(struct cfs_rq *cfs_rq)
{
	if (!cfs_bandwidth_used())
		return;

	/* an active group must be handled by the update_curr()->put() path */
	if (!cfs_rq->runtime_enabled || cfs_rq->curr)
		return;

	/* ensure the group is not already throttled */
	if (cfs_rq_throttled(cfs_rq))
		return;

	/* update runtime allocation */
	account_cfs_rq_runtime(cfs_rq, 0);
	if (cfs_rq->runtime_remaining <= 0)
		throttle_cfs_rq(cfs_rq);
}

static void sync_throttle(struct task_group *tg, int cpu)
{
	struct cfs_rq *pcfs_rq, *cfs_rq;

	if (!cfs_bandwidth_used())
		return;

	if (!tg->parent)
		return;

	cfs_rq = tg->cfs_rq[cpu];
	pcfs_rq = tg->parent->cfs_rq[cpu];

	cfs_rq->throttle_count = pcfs_rq->throttle_count;
	cfs_rq->throttled_clock_pelt = rq_clock_pelt(cpu_rq(cpu));
}

/* conditionally throttle active cfs_rq's from put_prev_entity() */
static bool check_cfs_rq_runtime(struct cfs_rq *cfs_rq)
{
	if (!cfs_bandwidth_used())
		return false;

	if (likely(!cfs_rq->runtime_enabled || cfs_rq->runtime_remaining > 0))
		return false;

	/*
	 * it's possible for a throttled entity to be forced into a running
	 * state (e.g. set_curr_task), in this case we're finished.
	 */
	if (cfs_rq_throttled(cfs_rq))
		return true;

	throttle_cfs_rq(cfs_rq);
	return true;
}

static enum hrtimer_restart sched_cfs_slack_timer(struct hrtimer *timer)
{
	struct cfs_bandwidth *cfs_b =
		container_of(timer, struct cfs_bandwidth, slack_timer);

	do_sched_cfs_slack_timer(cfs_b);

	return HRTIMER_NORESTART;
}

extern const u64 max_cfs_quota_period;

static enum hrtimer_restart sched_cfs_period_timer(struct hrtimer *timer)
{
	struct cfs_bandwidth *cfs_b =
		container_of(timer, struct cfs_bandwidth, period_timer);
	int overrun;
	int idle = 0;
	int count = 0;

	raw_spin_lock(&cfs_b->lock);
	for (;;) {
		overrun = hrtimer_forward_now(timer, cfs_b->period);
		if (!overrun)
			break;

		if (++count > 3) {
			u64 new, old = ktime_to_ns(cfs_b->period);

			/*
			 * Grow period by a factor of 2 to avoid losing precision.
			 * Precision loss in the quota/period ratio can cause __cfs_schedulable
			 * to fail.
			 */
			new = old * 2;
			if (new < max_cfs_quota_period) {
				cfs_b->period = ns_to_ktime(new);
				cfs_b->quota *= 2;

				pr_warn_ratelimited(
	"cfs_period_timer[cpu%d]: period too short, scaling up (new cfs_period_us = %lld, cfs_quota_us = %lld)\n",
					smp_processor_id(),
					div_u64(new, NSEC_PER_USEC),
					div_u64(cfs_b->quota, NSEC_PER_USEC));
			} else {
				pr_warn_ratelimited(
	"cfs_period_timer[cpu%d]: period too short, but cannot scale up without losing precision (cfs_period_us = %lld, cfs_quota_us = %lld)\n",
					smp_processor_id(),
					div_u64(old, NSEC_PER_USEC),
					div_u64(cfs_b->quota, NSEC_PER_USEC));
			}

			/* reset count so we don't come right back in here */
			count = 0;
		}

		idle = do_sched_cfs_period_timer(cfs_b, overrun);
	}
	if (idle)
		cfs_b->period_active = 0;
	raw_spin_unlock(&cfs_b->lock);

	return idle ? HRTIMER_NORESTART : HRTIMER_RESTART;
}

void init_cfs_bandwidth(struct cfs_bandwidth *cfs_b)
{
	raw_spin_lock_init(&cfs_b->lock);
	cfs_b->runtime = 0;
	cfs_b->quota = RUNTIME_INF;
	cfs_b->period = ns_to_ktime(default_cfs_period());

	INIT_LIST_HEAD(&cfs_b->throttled_cfs_rq);
	hrtimer_init(&cfs_b->period_timer, CLOCK_MONOTONIC, HRTIMER_MODE_ABS_PINNED);
	cfs_b->period_timer.function = sched_cfs_period_timer;
	hrtimer_init(&cfs_b->slack_timer, CLOCK_MONOTONIC, HRTIMER_MODE_REL);
	cfs_b->slack_timer.function = sched_cfs_slack_timer;
	cfs_b->distribute_running = 0;
}

static void init_cfs_rq_runtime(struct cfs_rq *cfs_rq)
{
	cfs_rq->runtime_enabled = 0;
	INIT_LIST_HEAD(&cfs_rq->throttled_list);
	walt_init_cfs_rq_stats(cfs_rq);
}

void start_cfs_bandwidth(struct cfs_bandwidth *cfs_b)
{
	lockdep_assert_held(&cfs_b->lock);

	if (cfs_b->period_active)
		return;

	cfs_b->period_active = 1;
	hrtimer_forward_now(&cfs_b->period_timer, cfs_b->period);
	hrtimer_start_expires(&cfs_b->period_timer, HRTIMER_MODE_ABS_PINNED);
}

static void destroy_cfs_bandwidth(struct cfs_bandwidth *cfs_b)
{
	/* init_cfs_bandwidth() was not called */
	if (!cfs_b->throttled_cfs_rq.next)
		return;

	hrtimer_cancel(&cfs_b->period_timer);
	hrtimer_cancel(&cfs_b->slack_timer);
}

/*
 * Both these CPU hotplug callbacks race against unregister_fair_sched_group()
 *
 * The race is harmless, since modifying bandwidth settings of unhooked group
 * bits doesn't do much.
 */

/* cpu online calback */
static void __maybe_unused update_runtime_enabled(struct rq *rq)
{
	struct task_group *tg;

	lockdep_assert_held(&rq->lock);

	rcu_read_lock();
	list_for_each_entry_rcu(tg, &task_groups, list) {
		struct cfs_bandwidth *cfs_b = &tg->cfs_bandwidth;
		struct cfs_rq *cfs_rq = tg->cfs_rq[cpu_of(rq)];

		raw_spin_lock(&cfs_b->lock);
		cfs_rq->runtime_enabled = cfs_b->quota != RUNTIME_INF;
		raw_spin_unlock(&cfs_b->lock);
	}
	rcu_read_unlock();
}

/* cpu offline callback */
static void __maybe_unused unthrottle_offline_cfs_rqs(struct rq *rq)
{
	struct task_group *tg;

	lockdep_assert_held(&rq->lock);

	rcu_read_lock();
	list_for_each_entry_rcu(tg, &task_groups, list) {
		struct cfs_rq *cfs_rq = tg->cfs_rq[cpu_of(rq)];

		if (!cfs_rq->runtime_enabled)
			continue;

		/*
		 * clock_task is not advancing so we just need to make sure
		 * there's some valid quota amount
		 */
		cfs_rq->runtime_remaining = 1;
		/*
		 * Offline rq is schedulable till CPU is completely disabled
		 * in take_cpu_down(), so we prevent new cfs throttling here.
		 */
		cfs_rq->runtime_enabled = 0;

		if (cfs_rq_throttled(cfs_rq))
			unthrottle_cfs_rq(cfs_rq);
	}
	rcu_read_unlock();
}

#else /* CONFIG_CFS_BANDWIDTH */

static inline bool cfs_bandwidth_used(void)
{
	return false;
}

static inline u64 cfs_rq_clock_task(struct cfs_rq *cfs_rq)
{
	return rq_clock_task(rq_of(cfs_rq));
}

static void account_cfs_rq_runtime(struct cfs_rq *cfs_rq, u64 delta_exec) {}
static bool check_cfs_rq_runtime(struct cfs_rq *cfs_rq) { return false; }
static void check_enqueue_throttle(struct cfs_rq *cfs_rq) {}
static inline void sync_throttle(struct task_group *tg, int cpu) {}
static __always_inline void return_cfs_rq_runtime(struct cfs_rq *cfs_rq) {}

static inline int cfs_rq_throttled(struct cfs_rq *cfs_rq)
{
	return 0;
}

static inline int throttled_hierarchy(struct cfs_rq *cfs_rq)
{
	return 0;
}

static inline int throttled_lb_pair(struct task_group *tg,
				    int src_cpu, int dest_cpu)
{
	return 0;
}

void init_cfs_bandwidth(struct cfs_bandwidth *cfs_b) {}

#ifdef CONFIG_FAIR_GROUP_SCHED
static void init_cfs_rq_runtime(struct cfs_rq *cfs_rq) {}
#endif

static inline struct cfs_bandwidth *tg_cfs_bandwidth(struct task_group *tg)
{
	return NULL;
}
static inline void destroy_cfs_bandwidth(struct cfs_bandwidth *cfs_b) {}
static inline void update_runtime_enabled(struct rq *rq) {}
static inline void unthrottle_offline_cfs_rqs(struct rq *rq) {}

#endif /* CONFIG_CFS_BANDWIDTH */

/**************************************************
 * CFS operations on tasks:
 */

#ifdef CONFIG_SCHED_HRTICK
static void hrtick_start_fair(struct rq *rq, struct task_struct *p)
{
	struct sched_entity *se = &p->se;

	SCHED_WARN_ON(task_rq(p) != rq);

	if (rq->cfs.h_nr_running > 1) {
		u64 ran = se->sum_exec_runtime - se->prev_sum_exec_runtime;
		u64 slice = se->slice;
		s64 delta = slice - ran;

		if (delta < 0) {
			if (rq->curr == p)
				resched_curr(rq);
			return;
		}
		hrtick_start(rq, delta);
	}
}

/*
 * called from enqueue/dequeue and updates the hrtick when the
 * current task is from our class and nr_running is low enough
 * to matter.
 */
static void hrtick_update(struct rq *rq)
{
	struct task_struct *curr = rq->curr;

	if (!hrtick_enabled(rq) || curr->sched_class != &fair_sched_class)
		return;

	hrtick_start_fair(rq, curr);
}
#else /* !CONFIG_SCHED_HRTICK */
static inline void
hrtick_start_fair(struct rq *rq, struct task_struct *p)
{
}

static inline void hrtick_update(struct rq *rq)
{
}
#endif

static unsigned long capacity_of(int cpu);
#ifdef CONFIG_SCHED_WALT
bool __cpu_overutilized(int cpu, int delta)
{
	return (capacity_orig_of(cpu) * 1024) <
		((cpu_util(cpu) + delta) * sched_capacity_margin_up[cpu]);
}

bool cpu_overutilized(int cpu)
{
	return __cpu_overutilized(cpu, 0);
}
#else
bool cpu_overutilized(int cpu)
{
	return !fits_capacity(cpu_util(cpu), capacity_of(cpu));
}
#endif

#ifdef CONFIG_SCHED_WALT
static bool sd_overutilized(struct sched_domain *sd)
{
	return sd->shared->overutilized;
}

static void set_sd_overutilized(struct sched_domain *sd)
{
	trace_sched_overutilized(sd, sd->shared->overutilized, true);
	sd->shared->overutilized = true;
}

static void clear_sd_overutilized(struct sched_domain *sd)
{
	trace_sched_overutilized(sd, sd->shared->overutilized, false);
	sd->shared->overutilized = false;
}
#endif
#ifdef CONFIG_SCHED_WALT
static inline void update_overutilized_status(struct rq *rq)
{
#ifdef CONFIG_SCHED_WALT
	struct sched_domain *sd;

	rcu_read_lock();
	sd = rcu_dereference(rq->sd);
	if (sd && !sd_overutilized(sd) &&
	    cpu_overutilized(rq->cpu))
		set_sd_overutilized(sd);
	rcu_read_unlock();
#else
	if (!READ_ONCE(rq->rd->overutilized) && cpu_overutilized(rq->cpu)) {
		WRITE_ONCE(rq->rd->overutilized, SG_OVERUTILIZED);
	}
#endif
}
#else
static inline void update_overutilized_status(struct rq *rq) { }
#endif

/*
 * The enqueue_task method is called before nr_running is
 * increased. Here we update the fair scheduling stats and
 * then put the task into the rbtree:
 */
static void
enqueue_task_fair(struct rq *rq, struct task_struct *p, int flags)
{
	struct cfs_rq *cfs_rq;
	struct sched_entity *se = &p->se;
        int idle_h_nr_running = task_has_idle_policy(p);
	bool prefer_idle = sched_feat(EAS_PREFER_IDLE) ?
				(schedtune_prefer_idle(p) > 0) : 0;


	/*
	 * The code below (indirectly) updates schedutil which looks at
	 * the cfs_rq utilization to select a frequency.
	 * Let's add the task's estimated utilization to the cfs_rq's
	 * estimated utilization, before we update schedutil.
	 */
	util_est_enqueue(&rq->cfs, p);

	/*
	 * The code below (indirectly) updates schedutil which looks at
	 * the cfs_rq utilization to select a frequency.
	 * Let's update schedtune here to ensure the boost value of the
	 * current task is accounted for in the selection of the OPP.
	 *
	 * We do it also in the case where we enqueue a throttled task;
	 * we could argue that a throttled task should not boost a CPU,
	 * however:
	 * a) properly implementing CPU boosting considering throttled
	 *    tasks will increase a lot the complexity of the solution
	 * b) it's not easy to quantify the benefits introduced by
	 *    such a more complex solution.
	 * Thus, for the time being we go for the simple solution and boost
	 * also for throttled RQs.
	 */
	schedtune_enqueue_task(p, cpu_of(rq));

#ifdef CONFIG_SCHED_WALT
	p->misfit = !task_fits_max(p, rq->cpu);
#endif
	/*
	 * If in_iowait is set, the code below may not trigger any cpufreq
	 * utilization updates, so do it here explicitly with the IOWAIT flag
	 * passed.
	 */
	if (p->in_iowait && prefer_idle)
		cpufreq_update_util(rq, SCHED_CPUFREQ_IOWAIT);

	for_each_sched_entity(se) {
		if (se->on_rq)
			break;
		cfs_rq = cfs_rq_of(se);
		enqueue_entity(cfs_rq, se, flags);

		cfs_rq->h_nr_running++;
                cfs_rq->idle_h_nr_running += idle_h_nr_running;
		walt_inc_cfs_rq_stats(cfs_rq, p);

		/* end evaluation on encountering a throttled cfs_rq */
		if (cfs_rq_throttled(cfs_rq))
			goto enqueue_throttle;

		flags = ENQUEUE_WAKEUP;
	}

	for_each_sched_entity(se) {
		cfs_rq = cfs_rq_of(se);
		walt_inc_cfs_rq_stats(cfs_rq, p);

		update_load_avg(cfs_rq, se, UPDATE_TG);
		update_cfs_group(se);

		cfs_rq->h_nr_running++;
		cfs_rq->idle_h_nr_running += idle_h_nr_running;

		/* end evaluation on encountering a throttled cfs_rq */
		if (cfs_rq_throttled(cfs_rq))
			goto enqueue_throttle;
	}

enqueue_throttle:
	if (!se) {
		add_nr_running(rq, 1);
		inc_rq_walt_stats(rq, p);
		/*
		 * Since new tasks are assigned an initial util_avg equal to
		 * half of the spare capacity of their CPU, tiny tasks have the
		 * ability to cross the overutilized threshold, which will
		 * result in the load balancer ruining all the task placement
		 * done by EAS. As a way to mitigate that effect, do not account
		 * for the first enqueue operation of new tasks during the
		 * overutilized flag detection.
		 *
		 * A better way of solving this problem would be to wait for
		 * the PELT signals of tasks to converge before taking them
		 * into account, but that is not straightforward to implement,
		 * and the following generally works well enough in practice.
		 *
		 * If the task prefers idle cpu, and it also is the first
		 * task enqueued in this runqueue, then we don't check
		 * overutilized. Hopefully the cpu util will be back to
		 * normal before next overutilized check.
		 */
		if ((flags & ENQUEUE_WAKEUP) && !(prefer_idle && rq->nr_running == 1))
			update_overutilized_status(rq);
	}

	if (cfs_bandwidth_used()) {
		/*
		 * When bandwidth control is enabled; the cfs_rq_throttled()
		 * breaks in the above iteration can result in incomplete
		 * leaf list maintenance, resulting in triggering the assertion
		 * below.
		 */
		for_each_sched_entity(se) {
			cfs_rq = cfs_rq_of(se);

			if (list_add_leaf_cfs_rq(cfs_rq))
				break;
		}
	}

	assert_list_leaf_cfs_rq(rq);

	hrtick_update(rq);
}

static void set_next_buddy(struct sched_entity *se);

/*
 * The dequeue_task method is called before nr_running is
 * decreased. We remove the task from the rbtree and
 * update the fair scheduling stats:
 */
static void dequeue_task_fair(struct rq *rq, struct task_struct *p, int flags)
{
	struct cfs_rq *cfs_rq;
	struct sched_entity *se = &p->se;
	int task_sleep = flags & DEQUEUE_SLEEP;
	int idle_h_nr_running = task_has_idle_policy(p);

	/*
	 * The code below (indirectly) updates schedutil which looks at
	 * the cfs_rq utilization to select a frequency.
	 * Let's update schedtune here to ensure the boost value of the
	 * current task is not more accounted for in the selection of the OPP.
	 */
	schedtune_dequeue_task(p, cpu_of(rq));

	for_each_sched_entity(se) {
		cfs_rq = cfs_rq_of(se);
		dequeue_entity(cfs_rq, se, flags);

		cfs_rq->h_nr_running--;
                cfs_rq->idle_h_nr_running -= idle_h_nr_running;
		walt_dec_cfs_rq_stats(cfs_rq, p);

		/* end evaluation on encountering a throttled cfs_rq */
		if (cfs_rq_throttled(cfs_rq))
			goto dequeue_throttle;

		/* Don't dequeue parent if it has other entities besides us */
		if (cfs_rq->load.weight) {
			/* Avoid re-evaluating load for this entity: */
			se = parent_entity(se);
			/*
			 * Bias pick_next to pick a task from this cfs_rq, as
			 * p is sleeping when it is within its sched_slice.
			 */
			if (task_sleep && se && !throttled_hierarchy(cfs_rq))
				set_next_buddy(se);
			break;
		}
		flags |= DEQUEUE_SLEEP;
	}

	for_each_sched_entity(se) {
		cfs_rq = cfs_rq_of(se);
		walt_dec_cfs_rq_stats(cfs_rq, p);

		update_load_avg(cfs_rq, se, UPDATE_TG);
		update_cfs_group(se);

		cfs_rq->h_nr_running--;
		cfs_rq->idle_h_nr_running -= idle_h_nr_running;

		/* end evaluation on encountering a throttled cfs_rq */
		if (cfs_rq_throttled(cfs_rq))
			goto dequeue_throttle;

	}

dequeue_throttle:
	if (!se) {
		sub_nr_running(rq, 1);
		dec_rq_walt_stats(rq, p);
	}

	util_est_dequeue(&rq->cfs, p, task_sleep);
	hrtick_update(rq);
}

#ifdef CONFIG_SMP

/* Working cpumask for: load_balance, load_balance_newidle. */
DEFINE_PER_CPU(cpumask_var_t, load_balance_mask);
DEFINE_PER_CPU(cpumask_var_t, select_idle_mask);

#ifdef CONFIG_NO_HZ_COMMON
/*
 * per rq 'load' arrray crap; XXX kill this.
 */

/*
 * The exact cpuload calculated at every tick would be:
 *
 *   load' = (1 - 1/2^i) * load + (1/2^i) * cur_load
 *
 * If a CPU misses updates for n ticks (as it was idle) and update gets
 * called on the n+1-th tick when CPU may be busy, then we have:
 *
 *   load_n   = (1 - 1/2^i)^n * load_0
 *   load_n+1 = (1 - 1/2^i)   * load_n + (1/2^i) * cur_load
 *
 * decay_load_missed() below does efficient calculation of
 *
 *   load' = (1 - 1/2^i)^n * load
 *
 * Because x^(n+m) := x^n * x^m we can decompose any x^n in power-of-2 factors.
 * This allows us to precompute the above in said factors, thereby allowing the
 * reduction of an arbitrary n in O(log_2 n) steps. (See also
 * fixed_power_int())
 *
 * The calculation is approximated on a 128 point scale.
 */
#define DEGRADE_SHIFT		7

static const u8 degrade_zero_ticks[CPU_LOAD_IDX_MAX] = {0, 8, 32, 64, 128};
static const u8 degrade_factor[CPU_LOAD_IDX_MAX][DEGRADE_SHIFT + 1] = {
	{   0,   0,  0,  0,  0,  0, 0, 0 },
	{  64,  32,  8,  0,  0,  0, 0, 0 },
	{  96,  72, 40, 12,  1,  0, 0, 0 },
	{ 112,  98, 75, 43, 15,  1, 0, 0 },
	{ 120, 112, 98, 76, 45, 16, 2, 0 }
};

/*
 * Update cpu_load for any missed ticks, due to tickless idle. The backlog
 * would be when CPU is idle and so we just decay the old load without
 * adding any new load.
 */
static unsigned long
decay_load_missed(unsigned long load, unsigned long missed_updates, int idx)
{
	int j = 0;

	if (!missed_updates)
		return load;

	if (missed_updates >= degrade_zero_ticks[idx])
		return 0;

	if (idx == 1)
		return load >> missed_updates;

	while (missed_updates) {
		if (missed_updates % 2)
			load = (load * degrade_factor[idx][j]) >> DEGRADE_SHIFT;

		missed_updates >>= 1;
		j++;
	}
	return load;
}

static struct {
	cpumask_var_t idle_cpus_mask;
	atomic_t nr_cpus;
	int has_blocked;		/* Idle CPUS has blocked load */
	unsigned long next_balance;     /* in jiffy units */
	unsigned long next_blocked;	/* Next update of blocked load in jiffies */
} nohz ____cacheline_aligned;

#endif /* CONFIG_NO_HZ_COMMON */

/**
 * __cpu_load_update - update the rq->cpu_load[] statistics
 * @this_rq: The rq to update statistics for
 * @this_load: The current load
 * @pending_updates: The number of missed updates
 *
 * Update rq->cpu_load[] statistics. This function is usually called every
 * scheduler tick (TICK_NSEC).
 *
 * This function computes a decaying average:
 *
 *   load[i]' = (1 - 1/2^i) * load[i] + (1/2^i) * load
 *
 * Because of NOHZ it might not get called on every tick which gives need for
 * the @pending_updates argument.
 *
 *   load[i]_n = (1 - 1/2^i) * load[i]_n-1 + (1/2^i) * load_n-1
 *             = A * load[i]_n-1 + B ; A := (1 - 1/2^i), B := (1/2^i) * load
 *             = A * (A * load[i]_n-2 + B) + B
 *             = A * (A * (A * load[i]_n-3 + B) + B) + B
 *             = A^3 * load[i]_n-3 + (A^2 + A + 1) * B
 *             = A^n * load[i]_0 + (A^(n-1) + A^(n-2) + ... + 1) * B
 *             = A^n * load[i]_0 + ((1 - A^n) / (1 - A)) * B
 *             = (1 - 1/2^i)^n * (load[i]_0 - load) + load
 *
 * In the above we've assumed load_n := load, which is true for NOHZ_FULL as
 * any change in load would have resulted in the tick being turned back on.
 *
 * For regular NOHZ, this reduces to:
 *
 *   load[i]_n = (1 - 1/2^i)^n * load[i]_0
 *
 * see decay_load_misses(). For NOHZ_FULL we get to subtract and add the extra
 * term.
 */
static void cpu_load_update(struct rq *this_rq, unsigned long this_load,
			    unsigned long pending_updates)
{
	unsigned long __maybe_unused tickless_load = this_rq->cpu_load[0];
	int i, scale;

	this_rq->nr_load_updates++;

	/* Update our load: */
	this_rq->cpu_load[0] = this_load; /* Fasttrack for idx 0 */
	for (i = 1, scale = 2; i < CPU_LOAD_IDX_MAX; i++, scale += scale) {
		unsigned long old_load, new_load;

		/* scale is effectively 1 << i now, and >> i divides by scale */

		old_load = this_rq->cpu_load[i];
#ifdef CONFIG_NO_HZ_COMMON
		old_load = decay_load_missed(old_load, pending_updates - 1, i);
		if (tickless_load) {
			old_load -= decay_load_missed(tickless_load, pending_updates - 1, i);
			/*
			 * old_load can never be a negative value because a
			 * decayed tickless_load cannot be greater than the
			 * original tickless_load.
			 */
			old_load += tickless_load;
		}
#endif
		new_load = this_load;
		/*
		 * Round up the averaging division if load is increasing. This
		 * prevents us from getting stuck on 9 if the load is 10, for
		 * example.
		 */
		if (new_load > old_load)
			new_load += scale - 1;

		this_rq->cpu_load[i] = (old_load * (scale - 1) + new_load) >> i;
	}
}

/* Used instead of source_load when we know the type == 0 */
static unsigned long weighted_cpuload(struct rq *rq)
{
	return cfs_rq_runnable_load_avg(&rq->cfs);
}

#ifdef CONFIG_NO_HZ_COMMON
/*
 * There is no sane way to deal with nohz on smp when using jiffies because the
 * CPU doing the jiffies update might drift wrt the CPU doing the jiffy reading
 * causing off-by-one errors in observed deltas; {0,2} instead of {1,1}.
 *
 * Therefore we need to avoid the delta approach from the regular tick when
 * possible since that would seriously skew the load calculation. This is why we
 * use cpu_load_update_periodic() for CPUs out of nohz. However we'll rely on
 * jiffies deltas for updates happening while in nohz mode (idle ticks, idle
 * loop exit, nohz_idle_balance, nohz full exit...)
 *
 * This means we might still be one tick off for nohz periods.
 */

static void cpu_load_update_nohz(struct rq *this_rq,
				 unsigned long curr_jiffies,
				 unsigned long load)
{
	unsigned long pending_updates;

	pending_updates = curr_jiffies - this_rq->last_load_update_tick;
	if (pending_updates) {
		this_rq->last_load_update_tick = curr_jiffies;
		/*
		 * In the regular NOHZ case, we were idle, this means load 0.
		 * In the NOHZ_FULL case, we were non-idle, we should consider
		 * its weighted load.
		 */
		cpu_load_update(this_rq, load, pending_updates);
	}
}

/*
 * Called from nohz_idle_balance() to update the load ratings before doing the
 * idle balance.
 */
static void cpu_load_update_idle(struct rq *this_rq)
{
	/*
	 * bail if there's load or we're actually up-to-date.
	 */
	if (weighted_cpuload(this_rq))
		return;

	cpu_load_update_nohz(this_rq, READ_ONCE(jiffies), 0);
}

/*
 * Record CPU load on nohz entry so we know the tickless load to account
 * on nohz exit. cpu_load[0] happens then to be updated more frequently
 * than other cpu_load[idx] but it should be fine as cpu_load readers
 * shouldn't rely into synchronized cpu_load[*] updates.
 */
void cpu_load_update_nohz_start(void)
{
	struct rq *this_rq = this_rq();

	/*
	 * This is all lockless but should be fine. If weighted_cpuload changes
	 * concurrently we'll exit nohz. And cpu_load write can race with
	 * cpu_load_update_idle() but both updater would be writing the same.
	 */
	this_rq->cpu_load[0] = weighted_cpuload(this_rq);
}

/*
 * Account the tickless load in the end of a nohz frame.
 */
void cpu_load_update_nohz_stop(void)
{
	unsigned long curr_jiffies = READ_ONCE(jiffies);
	struct rq *this_rq = this_rq();
	unsigned long load;
	struct rq_flags rf;

	if (curr_jiffies == this_rq->last_load_update_tick)
		return;

	load = weighted_cpuload(this_rq);
	rq_lock(this_rq, &rf);
	update_rq_clock(this_rq);
	cpu_load_update_nohz(this_rq, curr_jiffies, load);
	rq_unlock(this_rq, &rf);
}
#else /* !CONFIG_NO_HZ_COMMON */
static inline void cpu_load_update_nohz(struct rq *this_rq,
					unsigned long curr_jiffies,
					unsigned long load) { }
#endif /* CONFIG_NO_HZ_COMMON */

static void cpu_load_update_periodic(struct rq *this_rq, unsigned long load)
{
#ifdef CONFIG_NO_HZ_COMMON
	/* See the mess around cpu_load_update_nohz(). */
	this_rq->last_load_update_tick = READ_ONCE(jiffies);
#endif
	cpu_load_update(this_rq, load, 1);
}

/*
 * Called from scheduler_tick()
 */
void cpu_load_update_active(struct rq *this_rq)
{
	unsigned long load = weighted_cpuload(this_rq);

	if (tick_nohz_tick_stopped())
		cpu_load_update_nohz(this_rq, READ_ONCE(jiffies), load);
	else
		cpu_load_update_periodic(this_rq, load);
}

/*
 * Return a low guess at the load of a migration-source CPU weighted
 * according to the scheduling class and "nice" value.
 *
 * We want to under-estimate the load of migration sources, to
 * balance conservatively.
 */
static unsigned long source_load(int cpu, int type)
{
	struct rq *rq = cpu_rq(cpu);
	unsigned long total = weighted_cpuload(rq);

	if (type == 0 || !sched_feat(LB_BIAS))
		return total;

	return min(rq->cpu_load[type-1], total);
}

/*
 * Return a high guess at the load of a migration-target CPU weighted
 * according to the scheduling class and "nice" value.
 */
static unsigned long target_load(int cpu, int type)
{
	struct rq *rq = cpu_rq(cpu);
	unsigned long total = weighted_cpuload(rq);

	if (type == 0 || !sched_feat(LB_BIAS))
		return total;

	return max(rq->cpu_load[type-1], total);
}

static unsigned long cpu_avg_load_per_task(int cpu)
{
	struct rq *rq = cpu_rq(cpu);
	unsigned long nr_running = READ_ONCE(rq->cfs.h_nr_running);
	unsigned long load_avg = weighted_cpuload(rq);

	if (nr_running)
		return load_avg / nr_running;

	return 0;
}

static void record_wakee(struct task_struct *p)
{
	/*
	 * Only decay a single time; tasks that have less then 1 wakeup per
	 * jiffy will not have built up many flips.
	 */
	if (time_after(jiffies, current->wakee_flip_decay_ts + HZ)) {
		current->wakee_flips >>= 1;
		current->wakee_flip_decay_ts = jiffies;
	}

	if (current->last_wakee != p) {
		current->last_wakee = p;
		current->wakee_flips++;
	}
}

/*
 * Detect M:N waker/wakee relationships via a switching-frequency heuristic.
 *
 * A waker of many should wake a different task than the one last awakened
 * at a frequency roughly N times higher than one of its wakees.
 *
 * In order to determine whether we should let the load spread vs consolidating
 * to shared cache, we look for a minimum 'flip' frequency of llc_size in one
 * partner, and a factor of lls_size higher frequency in the other.
 *
 * With both conditions met, we can be relatively sure that the relationship is
 * non-monogamous, with partner count exceeding socket size.
 *
 * Waker/wakee being client/server, worker/dispatcher, interrupt source or
 * whatever is irrelevant, spread criteria is apparent partner count exceeds
 * socket size.
 */
static int wake_wide(struct task_struct *p, int sibling_count_hint)
{
	unsigned int master = current->wakee_flips;
	unsigned int slave = p->wakee_flips;
	int llc_size = this_cpu_read(sd_llc_size);

	if (sibling_count_hint >= llc_size)
		return 1;

	if (master < slave)
		swap(master, slave);
	if (slave < llc_size || master < slave * llc_size)
		return 0;
	return 1;
}

/*
 * The purpose of wake_affine() is to quickly determine on which CPU we can run
 * soonest. For the purpose of speed we only consider the waking and previous
 * CPU.
 *
 * wake_affine_idle() - only considers 'now', it check if the waking CPU is
 *			cache-affine and is (or	will be) idle.
 *
 * wake_affine_weight() - considers the weight to reflect the average
 *			  scheduling latency of the CPUs. This seems to work
 *			  for the overloaded case.
 */
static int
wake_affine_idle(int this_cpu, int prev_cpu, int sync)
{
	/*
	 * If this_cpu is idle, it implies the wakeup is from interrupt
	 * context. Only allow the move if cache is shared. Otherwise an
	 * interrupt intensive workload could force all tasks onto one
	 * node depending on the IO topology or IRQ affinity settings.
	 *
	 * If the prev_cpu is idle and cache affine then avoid a migration.
	 * There is no guarantee that the cache hot data from an interrupt
	 * is more important than cache hot data on the prev_cpu and from
	 * a cpufreq perspective, it's better to have higher utilisation
	 * on one CPU.
	 */
	if (available_idle_cpu(this_cpu) && cpus_share_cache(this_cpu, prev_cpu))
		return available_idle_cpu(prev_cpu) ? prev_cpu : this_cpu;

	if (sync && cpu_rq(this_cpu)->nr_running == 1)
		return this_cpu;

	return nr_cpumask_bits;
}

static int
wake_affine_weight(struct sched_domain *sd, struct task_struct *p,
		   int this_cpu, int prev_cpu, int sync)
{
	s64 this_eff_load, prev_eff_load;
	unsigned long task_load;

	this_eff_load = target_load(this_cpu, sd->wake_idx);

	if (sync) {
		unsigned long current_load = task_h_load(current);

		if (current_load > this_eff_load)
			return this_cpu;

		this_eff_load -= current_load;
	}

	task_load = task_h_load(p);

	this_eff_load += task_load;
	if (sched_feat(WA_BIAS))
		this_eff_load *= 100;
	this_eff_load *= capacity_of(prev_cpu);

	prev_eff_load = source_load(prev_cpu, sd->wake_idx);
	prev_eff_load -= task_load;
	if (sched_feat(WA_BIAS))
		prev_eff_load *= 100 + (sd->imbalance_pct - 100) / 2;
	prev_eff_load *= capacity_of(this_cpu);

	/*
	 * If sync, adjust the weight of prev_eff_load such that if
	 * prev_eff == this_eff that select_idle_sibling() will consider
	 * stacking the wakee on top of the waker if no other CPU is
	 * idle.
	 */
	if (sync)
		prev_eff_load += 1;

	return this_eff_load < prev_eff_load ? this_cpu : nr_cpumask_bits;
}

static int wake_affine(struct sched_domain *sd, struct task_struct *p,
		       int this_cpu, int prev_cpu, int sync)
{
	int target = nr_cpumask_bits;

	if (sched_feat(WA_IDLE))
		target = wake_affine_idle(this_cpu, prev_cpu, sync);

	if (sched_feat(WA_WEIGHT) && target == nr_cpumask_bits)
		target = wake_affine_weight(sd, p, this_cpu, prev_cpu, sync);

	schedstat_inc(p->se.statistics.nr_wakeups_affine_attempts);
	if (target == nr_cpumask_bits)
		return prev_cpu;

	schedstat_inc(sd->ttwu_move_affine);
	schedstat_inc(p->se.statistics.nr_wakeups_affine);
	return target;
}

#ifdef CONFIG_SCHED_TUNE
struct reciprocal_value schedtune_spc_rdiv;

static long
schedtune_margin(unsigned long signal, long boost, long capacity)
{
	long long margin = 0;

	/*
	 * Signal proportional compensation (SPC)
	 *
	 * The Boost (B) value is used to compute a Margin (M) which is
	 * proportional to the complement of the original Signal (S):
	 *   M = B * (capacity - S)
	 * The obtained M could be used by the caller to "boost" S.
	 */
	if (boost >= 0) {
		if (capacity > signal) {
			margin  = capacity - signal;
			margin *= boost;
		}
	} else
		margin = -signal * boost;

	margin  = reciprocal_divide(margin, schedtune_spc_rdiv);

	if (boost < 0)
		margin *= -1;
	return margin;
}

inline long
schedtune_cpu_margin_with(unsigned long util, int cpu, struct task_struct *p)
{
	int boost = schedtune_cpu_boost_with(cpu, p);
	long margin;

	if (boost == 0)
		margin = 0;
	else
		margin = schedtune_margin(util, boost, capacity_orig_of(cpu));

	trace_sched_boost_cpu(cpu, util, margin);

	return margin;
}

long schedtune_task_margin(struct task_struct *task)
{
	int boost = schedtune_task_boost(task);
	unsigned long util;
	long margin;

	if (boost == 0)
		return 0;

	util = task_util_est(task);
	margin = schedtune_margin(util, boost, SCHED_CAPACITY_SCALE);

	return margin;
}

#else /* CONFIG_SCHED_TUNE */

inline long
schedtune_cpu_margin_with(unsigned long util, int cpu, struct task_struct *p)
{
	return 0;
}

#endif /* CONFIG_SCHED_TUNE */

static unsigned long cpu_util_without(int cpu, struct task_struct *p);

static unsigned long capacity_spare_without(int cpu, struct task_struct *p)
{
	return max_t(long, capacity_of(cpu) - cpu_util_without(cpu, p), 0);
}

/*
 * find_idlest_group finds and returns the least busy CPU group within the
 * domain.
 *
 * Assumes p is allowed on at least one CPU in sd.
 */
static struct sched_group *
find_idlest_group(struct sched_domain *sd, struct task_struct *p,
		  int this_cpu, int sd_flag)
{
	struct sched_group *idlest = NULL, *group = sd->groups;
	struct sched_group *most_spare_sg = NULL;
	unsigned long min_runnable_load = ULONG_MAX;
	unsigned long this_runnable_load = ULONG_MAX;
	unsigned long min_avg_load = ULONG_MAX, this_avg_load = ULONG_MAX;
	unsigned long most_spare = 0, this_spare = 0;
	int load_idx = sd->forkexec_idx;
	int imbalance_scale = 100 + (sd->imbalance_pct-100)/2;
	unsigned long imbalance = scale_load_down(NICE_0_LOAD) *
				(sd->imbalance_pct-100) / 100;

	if (sd_flag & SD_BALANCE_WAKE)
		load_idx = sd->wake_idx;

	do {
		unsigned long load, avg_load, runnable_load;
		unsigned long spare_cap, max_spare_cap;
		int local_group;
		int i;

		/* Skip over this group if it has no CPUs allowed */
		if (!cpumask_intersects(sched_group_span(group),
					&p->cpus_allowed))
			continue;

		local_group = cpumask_test_cpu(this_cpu,
					       sched_group_span(group));

		/*
		 * Tally up the load of all CPUs in the group and find
		 * the group containing the CPU with most spare capacity.
		 */
		avg_load = 0;
		runnable_load = 0;
		max_spare_cap = 0;

		for_each_cpu(i, sched_group_span(group)) {
			/* Bias balancing toward CPUs of our domain */
			if (local_group)
				load = source_load(i, load_idx);
			else
				load = target_load(i, load_idx);

			runnable_load += load;

			avg_load += cfs_rq_load_avg(&cpu_rq(i)->cfs);

			spare_cap = capacity_spare_without(i, p);

			if (spare_cap > max_spare_cap)
				max_spare_cap = spare_cap;
		}

		/* Adjust by relative CPU capacity of the group */
		avg_load = (avg_load * SCHED_CAPACITY_SCALE) /
					group->sgc->capacity;
		runnable_load = (runnable_load * SCHED_CAPACITY_SCALE) /
					group->sgc->capacity;

		if (local_group) {
			this_runnable_load = runnable_load;
			this_avg_load = avg_load;
			this_spare = max_spare_cap;
		} else {
			if (min_runnable_load > (runnable_load + imbalance)) {
				/*
				 * The runnable load is significantly smaller
				 * so we can pick this new CPU:
				 */
				min_runnable_load = runnable_load;
				min_avg_load = avg_load;
				idlest = group;
			} else if ((runnable_load < (min_runnable_load + imbalance)) &&
				   (100*min_avg_load > imbalance_scale*avg_load)) {
				/*
				 * The runnable loads are close so take the
				 * blocked load into account through avg_load:
				 */
				min_avg_load = avg_load;
				idlest = group;
			}

			if (most_spare < max_spare_cap) {
				most_spare = max_spare_cap;
				most_spare_sg = group;
			}
		}
	} while (group = group->next, group != sd->groups);

	/*
	 * The cross-over point between using spare capacity or least load
	 * is too conservative for high utilization tasks on partially
	 * utilized systems if we require spare_capacity > task_util(p),
	 * so we allow for some task stuffing by using
	 * spare_capacity > task_util(p)/2.
	 *
	 * Spare capacity can't be used for fork because the utilization has
	 * not been set yet, we must first select a rq to compute the initial
	 * utilization.
	 */
	if (sd_flag & SD_BALANCE_FORK)
		goto skip_spare;

	if (this_spare > task_util(p) / 2 &&
	    imbalance_scale*this_spare > 100*most_spare)
		return NULL;

	if (most_spare > task_util(p) / 2)
		return most_spare_sg;

skip_spare:
	if (!idlest)
		return NULL;

	/*
	 * When comparing groups across NUMA domains, it's possible for the
	 * local domain to be very lightly loaded relative to the remote
	 * domains but "imbalance" skews the comparison making remote CPUs
	 * look much more favourable. When considering cross-domain, add
	 * imbalance to the runnable load on the remote node and consider
	 * staying local.
	 */
	if ((sd->flags & SD_NUMA) &&
	    min_runnable_load + imbalance >= this_runnable_load)
		return NULL;

	if (min_runnable_load > (this_runnable_load + imbalance))
		return NULL;

	if ((this_runnable_load < (min_runnable_load + imbalance)) &&
	     (100*this_avg_load < imbalance_scale*min_avg_load))
		return NULL;

	return idlest;
}

/*
 * find_idlest_group_cpu - find the idlest CPU among the CPUs in the group.
 */
static int
find_idlest_group_cpu(struct sched_group *group, struct task_struct *p, int this_cpu)
{
	unsigned long load, min_load = ULONG_MAX;
	unsigned int min_exit_latency = UINT_MAX;
	u64 latest_idle_timestamp = 0;
	int least_loaded_cpu = this_cpu;
	int shallowest_idle_cpu = -1;
	int i;

	/* Check if we have any choice: */
	if (group->group_weight == 1)
		return cpumask_first(sched_group_span(group));

	/* Traverse only the allowed CPUs */
	for_each_cpu_and(i, sched_group_span(group), &p->cpus_allowed) {
		if (available_idle_cpu(i)) {
			struct rq *rq = cpu_rq(i);
			struct cpuidle_state *idle = idle_get_state(rq);
			if (idle && idle->exit_latency < min_exit_latency) {
				/*
				 * We give priority to a CPU whose idle state
				 * has the smallest exit latency irrespective
				 * of any idle timestamp.
				 */
				min_exit_latency = idle->exit_latency;
				latest_idle_timestamp = rq->idle_stamp;
				shallowest_idle_cpu = i;
			} else if ((!idle || idle->exit_latency == min_exit_latency) &&
				   rq->idle_stamp > latest_idle_timestamp) {
				/*
				 * If equal or no active idle state, then
				 * the most recently idled CPU might have
				 * a warmer cache.
				 */
				latest_idle_timestamp = rq->idle_stamp;
				shallowest_idle_cpu = i;
			}
		} else if (shallowest_idle_cpu == -1) {
			load = weighted_cpuload(cpu_rq(i));
			if (load < min_load) {
				min_load = load;
				least_loaded_cpu = i;
			}
		}
	}

	return shallowest_idle_cpu != -1 ? shallowest_idle_cpu : least_loaded_cpu;
}

static inline int find_idlest_cpu(struct sched_domain *sd, struct task_struct *p,
				  int cpu, int prev_cpu, int sd_flag)
{
	int new_cpu = cpu;

	if (!cpumask_intersects(sched_domain_span(sd), &p->cpus_allowed))
		return prev_cpu;

	/*
	 * We need task's util for capacity_spare_without, sync it up to
	 * prev_cpu's last_update_time.
	 */
	if (!(sd_flag & SD_BALANCE_FORK))
		sync_entity_load_avg(&p->se);

	while (sd) {
		struct sched_group *group;
		struct sched_domain *tmp;
		int weight;

		if (!(sd->flags & sd_flag)) {
			sd = sd->child;
			continue;
		}

		group = find_idlest_group(sd, p, cpu, sd_flag);
		if (!group) {
			sd = sd->child;
			continue;
		}

		new_cpu = find_idlest_group_cpu(group, p, cpu);
		if (new_cpu == cpu) {
			/* Now try balancing at a lower domain level of 'cpu': */
			sd = sd->child;
			continue;
		}

		/* Now try balancing at a lower domain level of 'new_cpu': */
		cpu = new_cpu;
		weight = sd->span_weight;
		sd = NULL;
		for_each_domain(cpu, tmp) {
			if (weight <= tmp->span_weight)
				break;
			if (tmp->flags & sd_flag)
				sd = tmp;
		}
	}

	return new_cpu;
}

#ifdef CONFIG_SCHED_SMT
DEFINE_STATIC_KEY_FALSE(sched_smt_present);
EXPORT_SYMBOL_GPL(sched_smt_present);

static inline void set_idle_cores(int cpu, int val)
{
	struct sched_domain_shared *sds;

	sds = rcu_dereference(per_cpu(sd_llc_shared, cpu));
	if (sds)
		WRITE_ONCE(sds->has_idle_cores, val);
}

static inline bool test_idle_cores(int cpu, bool def)
{
	struct sched_domain_shared *sds;

	if (static_branch_likely(&sched_smt_present)) {
		sds = rcu_dereference(per_cpu(sd_llc_shared, cpu));
		if (sds)
			return READ_ONCE(sds->has_idle_cores);
	}

	return def;
}

/*
 * Scans the local SMT mask to see if the entire core is idle, and records this
 * information in sd_llc_shared->has_idle_cores.
 *
 * Since SMT siblings share all cache levels, inspecting this limited remote
 * state should be fairly cheap.
 */
void __update_idle_core(struct rq *rq)
{
	int core = cpu_of(rq);
	int cpu;

	rcu_read_lock();
	if (test_idle_cores(core, true))
		goto unlock;

	for_each_cpu(cpu, cpu_smt_mask(core)) {
		if (cpu == core)
			continue;

		if (!available_idle_cpu(cpu))
			goto unlock;
	}

	set_idle_cores(core, 1);
unlock:
	rcu_read_unlock();
}

/*
 * Scan the entire LLC domain for idle cores; this dynamically switches off if
 * there are no idle cores left in the system; tracked through
 * sd_llc->shared->has_idle_cores and enabled through update_idle_core() above.
 */
static int select_idle_core(struct task_struct *p, struct sched_domain *sd, int target)
{
	struct cpumask *cpus = this_cpu_cpumask_var_ptr(select_idle_mask);
	int core, cpu;

	if (!static_branch_likely(&sched_smt_present))
		return -1;

	if (!test_idle_cores(target, false))
		return -1;

	cpumask_and(cpus, sched_domain_span(sd), &p->cpus_allowed);

	for_each_cpu_wrap(core, cpus, target) {
		bool idle = true;

		for_each_cpu(cpu, cpu_smt_mask(core)) {
			cpumask_clear_cpu(cpu, cpus);
			if (!available_idle_cpu(cpu))
				idle = false;
				break;
			}
		}
		cpumask_andnot(cpus, cpus, cpu_smt_mask(core));

		if (idle)
			return core;
	}

	/*
	 * Failed to find an idle core; stop looking for one.
	 */
	set_idle_cores(target, 0);

	return -1;
}

/*
 * Scan the local SMT mask for idle CPUs.
 */
static int select_idle_smt(struct task_struct *p, struct sched_domain *sd, int target)
{
	int cpu;

	if (!static_branch_likely(&sched_smt_present))
		return -1;

	for_each_cpu(cpu, cpu_smt_mask(target)) {
		if (!cpumask_test_cpu(cpu, &p->cpus_allowed))
			continue;
		if (available_idle_cpu(cpu))
			return cpu;
	}

	return -1;
}

#else /* CONFIG_SCHED_SMT */

static inline int select_idle_core(struct task_struct *p, struct sched_domain *sd, int target)
{
	return -1;
}

static inline int select_idle_smt(struct task_struct *p, struct sched_domain *sd, int target)
{
	return -1;
}

#endif /* CONFIG_SCHED_SMT */

/*
 * Scan the LLC domain for idle CPUs; this is dynamically regulated by
 * comparing the average scan cost (tracked in sd->avg_scan_cost) against the
 * average idle time for this rq (as found in rq->avg_idle).
 */
static int select_idle_cpu(struct task_struct *p, struct sched_domain *sd, int target)
{
	struct cpumask *cpus = this_cpu_cpumask_var_ptr(select_idle_mask);
	struct sched_domain *this_sd;
	u64 time, cost;
	s64 delta;
	int cpu, nr = INT_MAX;
	int this = smp_processor_id();

	this_sd = rcu_dereference(*this_cpu_ptr(&sd_llc));
	if (!this_sd)
		return -1;
	
	cpumask_and(cpus, sched_domain_span(sd), &p->cpus_allowed);

	if (sched_feat(SIS_PROP)) {
		u64 avg_cost, avg_idle, span_avg;

		/*
		 * Due to large variance we need a large fuzz factor;
		 * hackbench in particularly is sensitive here.
		 */
		avg_idle = this_rq()->avg_idle / 512;
		avg_cost = this_sd->avg_scan_cost + 1;

		span_avg = sd->span_weight * avg_idle;
		if (span_avg > 4*avg_cost)
			nr = div_u64(span_avg, avg_cost);
		else
			nr = 4;

		time = cpu_clock(this);
	}

	for_each_cpu_wrap(cpu, cpus, target) {
		if (!--nr)
			return -1;
		if (!cpumask_test_cpu(cpu, &p->cpus_allowed))
			continue;
	/* Don't need to rebalance while attached to NULL domain */
#ifdef CONFIG_SCHED_WALT
		if (cpu_isolated(cpu))
			continue;
#endif
		if (available_idle_cpu(cpu))
			break;
	}

	if (sched_feat(SIS_PROP)) {
		time = cpu_clock(this) - time;
		cost = this_sd->avg_scan_cost;
		delta = (s64)(time - cost) / 8;
		this_sd->avg_scan_cost += delta;
	}

	return cpu;
}

/*
