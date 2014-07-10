function yd = get_yd(time)
% Gets the year day equivalent of a time in matlab format
%
% Justine McMillan
% Dec 8, 2013

timevec = datevec(time(1));
yd = time - datenum(timevec(1),0,0);