#!/bin/bash
HISTFILE=~/.bash_history
set -a

echo -e "All dates are in the format
year-month-day hour:minute:second
1970-01-30 00:00:00\n"

read -ep "2D or 3D run: " run_type
history -s $run_type

read -ep "Number of Processors: " processors
history -s $processors

read -ep "Gridname: " gridName
history -s $gridName

read -ep "START_DATE: " start_date
history -s $start_date

read -ep "END_DATE : " end_date
history -s $end_date

read -ep "coldstart or hotstart: " startup_type
history -s $startup_type

if [ $startup_type = 'hotstart' ]; then
    read -ep 'STARTUP_FILE' startup_file
    history -s $startup_file
    startup_uv_type='set values'
    startup_turb_type='set values'
else
    startup_file='none'
    startup_uv_type='default'
    startup_turb_type='default'
fi

read -ep "EXTSTEP_SECONDS: " extstep_seconds
history -s $extstep_seconds


read -ep "RST_FIRST_OUT: " rst_first_out
history -s $rst_first_out

read -ep "RST_OUTPUT_STACK: " rst_output_stack
history -s $rst_output_stack

read -ep "NC_FIRST_OUT: " nc_first_out
history -s $nc_first_out

read -ep "NC_OUT_INTERVAL: " nc_out_interval
history -s $nc_out_interval

if [ $run_type = '3D' ]; then
    nc_velocity='T'
    exe='fvcom3d_fast'
else
    nc_velocity='F'
    exe='fvcom2d_fast'
fi

read -ep "PROBES_ON (T or F): " probes_on
history -s $probes_on

if [ $probes_on = 'T' ]; then
    read -ep "PROBES_NUMBER: " probes_number
    history -s $probes_number
    read -ep "PROBES_FILE: " probes_file
    history -s $probes_file
else
    probes_number='none'
    probes_file='none'
fi


read -ep "TURBINE_ON (T or F): " turbine_on
history -s $turbine_on

if [ $turbine_on = 'T' ]; then
    read -ep "TURBINE_FILE: " turbine_file
    history -s $turbine_file
else
    turbine_file='none'
fi

python makeRun.py
python runpar.py
