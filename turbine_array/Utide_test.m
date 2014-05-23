clear all

time=735604:1/(24*10):735964;

time_origin=datenum('1899-12-31 12:00:00');

time=time-time(1)+time_origin;

load('ut_constants.mat','const','shallow');

%period=1/const.freq(48)*3600
%period = 44714.16;
amp =3.0;

phase = 23;

lat=45.5;

period=1./const.freq(10:90)*3600;

for jj =1:length(period)

time_series=amp*cos((time-time_origin)*2*pi/(period(jj)/(24*3600))-phase*pi/180);

coef = ut_solv(time,time_series ,[] , lat, 'auto','NoTrend','Rmin',1,'OLS','NoDiagn','LinCI');

amp_err(jj)=amp-coef.A(1);

phase_err(jj)=phase-coef.g(1);

ts_recon=ut_reconstr(time, coef);

err(jj)= sqrt(mean((time_series-ts_recon).^2));


end


figure
plot(amp_err)
hold on
plot(phase_err,'r')
plot(err,'k')
hold off
