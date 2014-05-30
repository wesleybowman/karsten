time = [734776:1:734795]
y = 15*cos(time)
x = y.*y
lat = 44.86061
coef = ut_solv(time, double(x), [], lat, 'auto','NoTrend','Rmin',1,'OLS','NoDiagn','LinCI');
save elevCoef.mat coef