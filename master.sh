gridName = 'dngrid'
python makeFirstRun.py $gridName 2012

export gridName
./RUN_PAR2.sh
wait ${!}

for number in {02..12}
do
    python makeRun.py $gridName 2012 $number
    ./RUN_PAR2.sh
    wait ${!}
done
