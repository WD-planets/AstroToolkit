make clean

rm -rf source/auto_tutorials/

if [ "$1" = "wipe" ]; then
    echo "Wiping all .fits and .fits.gz files under ./source/tutorials/..."
    find source/tutorials/ \
        -type f \( -name "*.fits" -o -name "*.fits.gz" \) \
        ! -path "source/tutorials/extension/external_*.fits" \
        ! -path "source/tutorials/extension/external_*.fits.gz" \
        -exec rm -f {} +
fi

make html

# run again to ensure no argument change warnings, and to ensure includes are built
make clean
make html
