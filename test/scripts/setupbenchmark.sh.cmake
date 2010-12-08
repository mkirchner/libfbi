if [ -e "@LIBFBI_BINARY_DIR@/test/testdata.zip" ]
then
    unzip -n "@LIBFBI_BINARY_DIR@/test/testdata.zip" -d "@LIBFBI_BINARY_DIR@/test/"
fi
@LIBFBI_BINARY_DIR@/test/scripts/benchmark.sh 1 9
