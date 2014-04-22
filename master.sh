#!/bin/bash
gridName="$1"
year="$2"
month="$3"
nextMonth=$((month+1))
echo "0$nextMonth"

python makeFirstRun.py $gridName $year

export gridName
./RUN_PAR2.sh
wait ${!}
python moveFile $gridName $year $month

for number in {02..12}
do
    python makeRun.py $gridName $year $number
    ./RUN_PAR2.sh
    wait ${!}
done
