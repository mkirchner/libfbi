#!/bin/bash
cd @LIBFBI_BINARY_DIR@/test/testdata

echo -e "Points\tTime" > results.txt

for i in `seq -f "%02g" $1 $2`;
do
    echo "" > ${i}/results.txt
    for j in {03..13};
    do
        @LIBFBI_BINARY_DIR@/examples/example-xic-construction --mzLow=0 --mzHigh=$(echo "(${j}+ 1) * 1400 / 14" | bc -l) --snLow=0 --snHigh=$(echo "(${j} + 1) * 2000 / 14" | bc -l) ${i}/pdc >> ${i}/results.txt
    done;
    cat ${i}/results.txt >> results.txt
done;

R --slave < @LIBFBI_BINARY_DIR@/test/scripts/plotting.r "results.txt"
