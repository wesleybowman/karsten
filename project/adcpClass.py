import numpy as np
import scipy.io as sio


class ADCP:

    def __init__(self, filename):
        self.mat = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)


        self.lat = self.mat['lat']
        self.lon = self.mat['lon']

        self.data = self.mat['data']
        '''
        keys for data:
            bins, north_vel, east_vel, vert_vel, dir_vel, mag_signed_vel,
            ucross, ualong
        '''

        self.pressure = self.mat['pres']
        '''
        keys for pres:
            surf
        '''

        self.time = self.mat['time']
        '''
        keys for time:
            mtime
        '''

if __name__ == '__main__':
    filename = 'Flow_GP-100915-BPa.mat'
    data = ADCP(filename)
