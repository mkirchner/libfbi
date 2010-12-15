#!/bin/bash

usage() {
    echo "usage: $0 <fbi|kdtree|clean>"
    exit
}

if [ $# -ne 1 ] ; then
    usage
fi

# unpack benchmark data set
if [ ! -d "@LIBFBI_BINARY_DIR@/test/testdata" ]; then 
    echo "Unpacking benchmark test data."
    unzip -n "@LIBFBI_BINARY_DIR@/test/testdata.zip" \
      -d "@LIBFBI_BINARY_DIR@/test/"
else
    echo "Found benchmark data in @LIBFBI_BINARY_DIR@/test/testdata."
fi

cd "@LIBFBI_BINARY_DIR@/test/testdata"

case $1 in
  clean)
    echo "Deleting all benchmark data."
    rm -rf "@LIBFBI_BINARY_DIR@/test/testdata"
    ;;
  fbi)
    # run the benchmark for libfbi
    echo -e "Points\tTime" > fbi-results.txt
    for i in `seq -f "%02g" 1 9`; do # should be 1-9
        echo "Benchmark dataset ${i}:"
        rm -f ${i}/fbi-results.txt
        for j in {03..13};
        do
            MZ_HIGH=$(echo "(${j} + 1) * 1400 / 14" | bc -l)
            SN_HIGH=$(echo "(${j} + 1) * 2000 / 14" | bc -l)
            echo "... subset (m/z=(0-$MZ_HIGH), sn=(0-$SN_HIGH))"
            @LIBFBI_BINARY_DIR@/examples/example-xic-construction --mzLow=0 \
              --mzHigh=$MZ_HIGH --snLow=0 --snHigh=$SN_HIGH \
              ${i}/pdc >> ${i}/fbi-results.txt
        done;
        cat ${i}/fbi-results.txt >> fbi-results.txt
    done;
    # and generate the graph
    R --slave < @LIBFBI_BINARY_DIR@/test/scripts/plotting.r fbi-results.txt benchmark-fbi.pdf
    ;;
  kdtree)
    # run the benchmark for the kdtree
    echo -e "Points\tTime" > kdtree-results.txt
    for i in `seq -f "%02g" 1 9`; do # should be 1-9
        echo "Benchmark dataset ${i}:"
        rm -f ${i}/kdtree-results.txt
        for j in {03..13};
        do
            MZ_HIGH=$(echo "(${j} + 1) * 1400 / 14" | bc -l)
            SN_HIGH=$(echo "(${j} + 1) * 2000 / 14" | bc -l)
            echo "... subset (m/z=(0-$MZ_HIGH), sn=(0-$SN_HIGH))"
            @LIBFBI_BINARY_DIR@/examples/kdtree-xic-construction --mzLow=0 \
              --mzHigh=$MZ_HIGH --snLow=0 --snHigh=$SN_HIGH \
              ${i}/pdc >> ${i}/kdtree-results.txt
        done;
        cat ${i}/kdtree-results.txt >> kdtree-results.txt
    done;
    R --slave < @LIBFBI_BINARY_DIR@/test/scripts/plotting.r kdtree-results.txt benchmark-kdtree.pdf
    ;;
 *)
    usage
    ;;
esac
