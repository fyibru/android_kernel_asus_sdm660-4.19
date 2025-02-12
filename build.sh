error=false
caused=""
C=/workspace/clang/bin
export TZ='Asia/Jakarta'
BUILDDATE=$(date +%Y%m%d)
BUILDTIME=$(date +%H%M%S)

function ret_err(){
        error=true
        caused="error"  
        $error
}

function buildAndUpload(){
        # remove old config
        if [ -f out/.config ]; then
                rm out/.config
        fi
        # build defconfig
        make -j8 O=out ARCH=arm64 CC="ccache $C/clang" CROSS_COMPILE=aarch64-linux-gnu- \
                                  LD=$C/ld.lld AR=$C/llvm-ar NM=$C/llvm-nm STRIP=$C/llvm-strip \
                                  OBJCOPY=$C/llvm-objcopy OBJDUMP=$C/llvm-objdump READELF=$C/llvm-readelf \
                                  HOSTCC=$C/clang HOSTCXX=$C/clang++ HOSTAR=$C/llvm-ar HOSTLD=$C/ld.lld \
                                  asus/X01BD_defconfig

        make -j8 O=out ARCH=arm64 CC="ccache $C/clang" CROSS_COMPILE=aarch64-linux-gnu- \
                                  LD=$C/ld.lld AR=$C/llvm-ar NM=$C/llvm-nm STRIP=$C/llvm-strip \
                                  OBJCOPY=$C/llvm-objcopy OBJDUMP=$C/llvm-objdump READELF=$C/llvm-readelf \
                                  HOSTCC=$C/clang HOSTCXX=$C/clang++ HOSTAR=$C/llvm-ar HOSTLD=$C/ld.lld
}

function update(){
        sudo apt update -y
        sudo apt install cpio flex gcc-aarch64-linux-gnu bc ccache -y
}

update
buildAndUpload