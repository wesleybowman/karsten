import shutil
import sys

arguments = sys.argv

gridName = arguments[1]
year = arguments[2]
month = arguments[3]

fileName = './output/{0}_0001.nc'.format(gridName)
newFile = './output/{0}_{1}_{2}.nc'.format(gridName, year, month)
shutil.move(fileName, newFile)

fileName = './output/{0}_restart_0005.nc'.format(gridName)
newFile = './input/{0}_restart_0005.nc'.format(gridName)
shutil.move(fileName, newFile)


restartFiles = [0001, 0002, 0003, 0004]
for i in restartFiles:
    restartFile = './output/{0}_restart_{3}.nc'.format(gridName, i)
    newRestartFile = './output/{0}_restart_{1}_{2}_{3}.nc'.format(gridName, year, month, i)
    shutil.move(restartFile, newRestartFile)
