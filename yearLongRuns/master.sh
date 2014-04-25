#!/bin/bash
gridName="$1"
year="$2"

python makeFirstRun.py $gridName $year

export gridName
python myQSUB.py
python moveFile.py $gridName $year 01
wait ${!}

for month in {02..12}
do
    python makeRun.py $gridName $year $month
    python myQSUB.py
    python moveFile.py $gridName $year $month
    wait ${!}
done
