#!/bin/bash
set -a

echo -e "All dates are in the format
year-month-day hour:minute:second
1970-01-30 00:00:00\n"

echo -n "2D or 3D run: "
read run_type

echo -n "Number of Processors: "
read processors

echo -n "Gridname: "
read gridName

echo -n "START_DATE: "
read start_date

echo -n "END_DATE : "
read end_date

echo -n "coldstart or hotstart: "
read startup_type

if [ $startup_type = 'hotstart' ]; then
    echo -n 'STARTUP_FILE'
    read startup_file
    startup_uv_type='set values'
    startup_turb_type='set values'
else
    startup_file='none'
    startup_uv_type='default'
    startup_turb_type='default'
fi

echo -n "EXTSTEP_SECONDS: "
read extstep_seconds


echo -n "RST_FIRST_OUT: "
read rst_first_out

echo -n "RST_OUTPUT_STACK: "
read rst_output_stack

echo -n "NC_FIRST_OUT: "
read nc_first_out

echo -n "NC_OUT_INTERVAL: "
read nc_out_interval

if [ $run_type = '3D' ]; then
    nc_velocity='T'
    exe='fvcom3d_fast'
else
    nc_velocity='F'
    exe='fvcom2d_fast'
fi

echo -n "PROBES_ON (T or F): "
read probes_on

if [ $probes_on = 'T' ]; then
    echo -n "PROBES_NUMBER: "
    read probes_number
    echo -n "PROBES_FILE: "
    read probes_file
else
    probes_number='none'
    probes_file='none'
fi


echo -n "TURBINE_ON (T or F): "
read turbine_on

if [ $turbine_on = 'T' ]; then
    echo -n "TURBINE_FILE: "
    read turbine_file
else
    turbine_file='none'
fi

python makeRun.py
python runpar.py
#RUN_PAR2.sh>>RUNTEST.sh

#gridName="$1"
#year="$2"
#
#python makeFirstRun.py $gridName $year
#
#export gridName
#python myQSUB.py
#python moveFile.py $gridName $year 01
#wait ${!}
#
#for month in {02..12}
#do
#    python makeRun.py $gridName $year $month
#    python myQSUB.py
#    python moveFile.py $gridName $year $month
#    wait ${!}
#done
