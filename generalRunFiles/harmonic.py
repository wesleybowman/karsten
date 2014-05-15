import sys
sys.path.append('/home/wesley/github/')
from utide import ut_solv

coef = ut_solv(time, ua[:,loci[ii]], va[:,loci[ii]], uvnodell[loci[ii],1],
                'auto', Rayleigh[0],'NoTrend','Rmin', 'OLS',
                'NoDiagn', 'LinCI')
