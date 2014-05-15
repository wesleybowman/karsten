import shutil


standardRunDir = './test/'
gridName = 'dngrid'
start_date = '2012-01-01 00:00:00'
end_date = '2012-02-01 00:00:00'
startup_type = 'coldstart'
startup_file = 'none'
startup_uv_type = 'default'
startup_turb_type = 'default'
extstep_seconds = 0.6
rst_first_out = '2012-01-01 00:00:00'
rst_output_stack = 1
nc_first_out = '2012-01-01 00:00:00'
nc_out_interval = '600s'
nc_velocity = 'T'
probes_on = 'F'
probes_number = 75
probes_file = 'none'
turbine_on = 'F'
turbine_file = 'none'
bottom_roughness_minimum = ['0.025', '0.030']


for bottom_roughness in bottom_roughness_minimum:
    top = '''
!================================================================!
_______  _     _  _______  _______  _______  ______     _____
(_______)(_)   (_)(_______)(_______)(_______)(_____ \   (_____)
_____    _     _  _        _     _  _  _  _  _____) )  _  __ _
|  ___)  | |   | || |      | |   | || ||_|| |(_____ (  | |/ /| |
| |       \ \ / / | |_____ | |___| || |   | | _____) )_|   /_| |
|_|        \___/   \______) \_____/ |_|   |_|(______/(_)\_____/
-- Beta Release
!                                                                !
!========DOMAIN DECOMPOSITION USING: METIS 4.0.1 ================!
!======Copyright 1998, Regents of University of Minnesota========!
!                                                                !


&NML_CASE
CASE_TITLE      = '{0}'
TIMEZONE        = 'UTC',
DATE_FORMAT     = 'YMD'
START_DATE      = '{1}'
END_DATE        = '{2}'

&NML_STARTUP
STARTUP_TYPE      = '{3}'
STARTUP_FILE      = '{4}'
STARTUP_UV_TYPE   = '{5}'
STARTUP_TURB_TYPE = '{6}'
STARTUP_TS_TYPE   = 'constant'
STARTUP_T_VALS    = 18
STARTUP_S_VALS    = 35.0
STARTUP_DMAX      =  -10.0
/

&NML_IO
INPUT_DIR       =  './input/'
OUTPUT_DIR      =  './output'
IREPORT         =  720,
VISIT_ALL_VARS  = F,
WAIT_FOR_VISIT  = F,
USE_MPI_IO_MODE = F
/

&NML_INTEGRATION
EXTSTEP_SECONDS =  {7},
ISPLIT          =  1
IRAMP           =  34560
MIN_DEPTH       =  0.5
STATIC_SSH_ADJ  =  0.0
/

&NML_RESTART
RST_ON  = T,
RST_FIRST_OUT      = '{8}'
RST_OUT_INTERVAL   = 'days = 1.0'
RST_OUTPUT_STACK   =           {9}
/

&NML_NETCDF
NC_ON   = T,
NC_FIRST_OUT    = '{10}',
NC_OUT_INTERVAL =  'seconds={11}',
NC_OUTPUT_STACK =  0,
NC_GRID_METRICS = T,
NC_VELOCITY     = {12},
NC_SALT_TEMP    = F,
NC_TURBULENCE   = F,
NC_AVERAGE_VEL  = T,
NC_VERTICAL_VEL = F,
NC_WIND_VEL     = F,
NC_WIND_STRESS  = F,
NC_EVAP_PRECIP  = F,
NC_SURFACE_HEAT = F,
NC_GROUNDWATER = F
/

&NML_NETCDF_AV
NCAV_ON = F,
NCAV_FIRST_OUT  = 'none'
NCAV_OUT_INTERVAL       =  0.0,
NCAV_OUTPUT_STACK       =           0,
NCAV_GRID_METRICS       = F,
NCAV_FILE_DATE  = F,
NCAV_VELOCITY   = F,
NCAV_SALT_TEMP  = F,
NCAV_TURBULENCE = F,
NCAV_AVERAGE_VEL        = F,
NCAV_VERTICAL_VEL       = F,
NCAV_WIND_VEL   = F,
NCAV_WIND_STRESS        = F,
NCAV_EVAP_PRECIP        = F,
NCAV_SURFACE_HEAT       = F,
NCAV_GROUNDWATER        = F,
NCAV_BIO        = F,
NCAV_WQM        = F,
NCAV_VORTICITY  = F
/

&NML_SURFACE_FORCING
WIND_ON = F,
HEATING_ON      = F,
PRECIPITATION_ON        = F,
/

&NML_PHYSICS
HORIZONTAL_MIXING_TYPE          = 'closure'
HORIZONTAL_MIXING_KIND          = 'constant'
HORIZONTAL_MIXING_COEFFICIENT   = 0.3
HORIZONTAL_PRANDTL_NUMBER       = 1.0
VERTICAL_MIXING_TYPE            = 'closure'
VERTICAL_MIXING_COEFFICIENT     = 1.0E-3,
VERTICAL_PRANDTL_NUMBER         = 1.0
BOTTOM_ROUGHNESS_MINIMUM        =  {18}
BOTTOM_ROUGHNESS_LENGTHSCALE    =  0.001
BOTTOM_ROUGHNESS_KIND           = 'constant'
BOTTOM_ROUGHNESS_TYPE           = 'orig'
CONVECTIVE_OVERTURNING          = F,
SCALAR_POSITIVITY_CONTROL       = T,
BAROTROPIC                      = T,
BAROCLINIC_PRESSURE_GRADIENT    = 'sigma levels'
SEA_WATER_DENSITY_FUNCTION      = 'dens2'
RECALCULATE_RHO_MEAN           = F
INTERVAL_RHO_MEAN              = 'seconds=1800.'
TEMPERATURE_ACTIVE              = F,
SALINITY_ACTIVE                 = F,
SURFACE_WAVE_MIXING             = F,
WETTING_DRYING_ON               = T
/

&NML_RIVER_TYPE
RIVER_NUMBER    =           0,
/

&NML_OPEN_BOUNDARY_CONTROL
OBC_ON                      = T,
OBC_NODE_LIST_FILE          = '{0}_obc.dat'
OBC_ELEVATION_FORCING_ON    = T,
OBC_ELEVATION_FILE          = '{0}_el_obc.nc'
OBC_TS_TYPE                 = 3
OBC_TEMP_NUDGING            = F,
OBC_TEMP_FILE               = 'none'
OBC_TEMP_NUDGING_TIMESCALE  =  0.0000000E+00,
OBC_SALT_NUDGING            = F,
OBC_SALT_FILE               = 'none'
OBC_SALT_NUDGING_TIMESCALE  =  0.0000000E+00,
OBC_MEANFLOW                = F,
/

&NML_GRID_COORDINATES
GRID_FILE       = '{0}_grd.dat'
GRID_FILE_UNITS = 'meters'
PROJECTION_REFERENCE  = 'proj=lcc +lon_0=-64.55880 +lat_0=41.78504 +lat_1=39.69152 +lat_2=43.87856'
SIGMA_LEVELS_FILE     = 'sigma.dat'
DEPTH_FILE      = '{0}_dep.dat'
CORIOLIS_FILE   = '{0}_cor.dat'
SPONGE_FILE     = '{0}_spg.dat'
BFRIC_FILE='{0}_bfric.dat'
VVCOE_FILE='{0}_vvcoe.dat'
/

&NML_GROUNDWATER
GROUNDWATER_ON             = F,
GROUNDWATER_FLOW  = 0.0,
GROUNDWATER_FILE           = 'none'
/

&NML_LAG
LAG_PARTICLES_ON        = F,
LAG_START_FILE   = 'none'
LAG_OUT_FILE     = 'none'
LAG_RESTART_FILE = 'none'
LAG_OUT_INTERVAL =  0.000000000000000E+000,
LAG_SCAL_CHOICE  = 'none'
/

&NML_ADDITIONAL_MODELS
DATA_ASSIMILATION       = F,
BIOLOGICAL_MODEL        = F,
SEDIMENT_MODEL  = F,
SEDIMENT_PARAMETER_TYPE = 'constant'
SEDIMENT_MODEL_FILE     = 'generic_sediment.inp'
ICING_MODEL     = F,
ICE_MODEL       = F,
/

&NML_PROBES
PROBES_ON = {13},
PROBES_NUMBER = '{14}',
PROBES_FILE = '{15}',
/

&NML_TURBINE
TURBINE_ON = {16},
TURBINE_FILE = '{17}',
/

&NML_NESTING
NESTING_ON = F
/

&NML_NCNEST
NCNEST_ON = F
/

&NML_BOUNDSCHK
BOUNDSCHK_ON  = F
/

&NML_STATION_TIMESERIES
OUT_STATION_TIMESERIES_ON       = F,
STATION_FILE='NONE'
LOCATION_TYPE='NONE'
OUT_ELEVATION=F,
OUT_VELOCITY_3D=F,
OUT_VELOCITY_2D=F,
OUT_SALT_TEMP =F,
OUT_WIND_VELOCITY=F,
OUT_INTERVAL= 'seconds=1000.0'
/
    ''' .format(gridName, start_date, end_date, startup_type, startup_file,
                startup_uv_type, startup_turb_type, extstep_seconds, rst_first_out,
                rst_output_stack, nc_first_out, nc_out_interval, nc_velocity,
                probes_on, probes_number, probes_file, turbine_on,
                turbine_file, bottom_roughness)


    newFile = '{0}'.format(bottom_roughness)
    shutil.copytree(standardRunDir, newFile)
    outputFile = '{0}/{1}_run.nml'.format(newFile, gridName)

    top = top.split('\n')

    with open(outputFile, 'w') as f:
        for t in top:
            print >> f, t
