import cPickle as pickle
import scipy.io as sio
import numpy as np


matCoef = sio.loadmat('../speedcoef.m', struct_as_record=False, squeeze_me=True)
#matCoef = sio.loadmat('../elevcoef.m', struct_as_record=False, squeeze_me=True)
mcoef = matCoef['coef']

coef = pickle.load(open('speedcoef.p', 'rb'))
#coef = pickle.load(open('elevcoef.p', 'rb'))

for k in coef.keys():
    mat = eval('mcoef.' + k)
    pyth = coef[k]
    try:
        if not np.allclose(mat, pyth).all():
            print '\n' + k
            print mat
            print pyth
        #print np.allclose(mat, pyth)
    except TypeError:
        if k == 'name':
            names = []
            for i in coef[k]:
                names.append(i.replace(' ', ''))

            mat = eval('mcoef.' + k)
            ind = np.where( mat == names)[0].shape[0]
            size = len(names)
            #if not mat.all() == names:
            if not ind == size:
                print mat
                print pyth
                print mat == pyth

        else:
            for i in coef[k]:

                mat = eval('mcoef.' + k + '.' + i)
                pyth = coef[k][i]
                try:
                    if not np.allclose(mat, pyth).all():
                        print '\n' + k + ':' + i
                        print mat
                        print pyth
                    #print np.allclose(mat, pyth)
                except TypeError:
                    for j in coef[k][i]:
                        mat = eval('mcoef.' + k + '.' + i + '.' + j)
                        pyth = coef[k][i][j]
                        #print '\n' + k + ':' + i + ':' + j
                        #print mat
                        #print pyth
                        if type(pyth) == str:
                            #print mat == pyth
                            if not mat == pyth:
                                print mat
                                print pyth

                        else:
                            if not np.allclose(mat, pyth).all():
                                print '\n' + k + ':' + i + ':' + j
                                print mat
                                print pyth
                            #print np.allclose(mat, pyth)
