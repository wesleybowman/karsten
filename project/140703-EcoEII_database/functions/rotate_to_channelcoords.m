function [xcross,xalong]=rotate_to_channelcoords(x,y,theta)
% x,y are the coordinates (could be speeds) relative to true
% north -- inputs can be vectors
% theta is the angle between the old coordinate system and the new system
% (measured from x to xcross)
% xcross, xalong are the coordinates cross channel (xcross) and along channel (xalong)
% 
% May 3, 2012

%disp('rotating velocities to be along/across channel')
Theta = theta*pi/180;

xcross = x.*cos(Theta)+y.*sin(Theta);
xalong = -x.*sin(Theta)+y.*cos(Theta);