make clean

rm -rf source/auto_tutorials/

if [ "$1" = "wipe" ]; then
    echo "Wiping all .fits files under ./source/tutorials/..."
    find source/tutorials/ \
        -path "source/tutorials/extension" -prune -o \
        -type f -name "*.fits" -delete
fi

if [ "$1" = "wipe" ]; then
    echo "Wiping all .fits.gz files under ./source/tutorials/..."
    find source/tutorials/ \
        -path "source/tutorials/extension" -prune -o \
        -type f -name "*.fits.gz" -delete
fi

make html

# RUN AGAIN TO ENSURE NO "NEW ARGUMENT" WARNINGS + TO ENSURE INCLUDES ARE BUILT
make html
