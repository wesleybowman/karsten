#!/bin/bash
gridName="$1"
year="$2"

python makeFirstRun.py $gridName $year

export gridName
qsub RUN_PAR2.sh
wait ${!}
python moveFile $gridName $year 01

for month in {02..12}
do
    python makeRun.py $gridName $year $month
    qsub RUN_PAR2.sh
    wait ${!}
    python moveFile $gridName $year $month
done
