import shutil
import sys

arguments = sys.argv

gridName = arguments[1]
year = arguments[2]
month = arguments[3]

fileName = '{0}_0001.nc'
newFile = '{0}_{1}_{2}.nc'
shutil.move(fileName, newFile)

restartFiles = [0001, 0002, 0003, 0004]
for i in restartFiles:
    restartFile = '{0}_restart_{3}.nc'
    newRestartFile = '{0}_restart_{1}_{2}_{3}.nc'
    shutil.move(restartFile, newRestartFile)
