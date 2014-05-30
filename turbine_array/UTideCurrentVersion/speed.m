time = [734776:1:734795];
time=735604:1/(24*10):735964;
y = 15*cos(time);
x = y.*y;
lat = 44.86061;
coef = ut_solv(time, double(x), double(y), lat, 'auto','NoTrend','Rmin',1,'OLS','NoDiagn','LinCI');
save speedCoef.mat coef