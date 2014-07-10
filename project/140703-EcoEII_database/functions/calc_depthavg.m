function [avgfield, zout] = calc_depthavg(field,z,varargin)
% Computes the depth_average of a field
% Very simple function that averages over all the rows that are inputed
% Assumes dz is constant
% 
% Inputs:
% field: field to average 
% z: vector containing depths
% varargin
%    1) depth limits (compute average over a subset of the depth)
%    2) dim (dimension that corresponds to z, default = 1)
%
% examples
%   avgfield = calc_depthavg(field,z,[2.1 14],2)
%   avgfield = calc_depthavg(field,z,[],2)
%
% Justine McMillan
% October 2nd, 2012
%
% Changes
% Feb 4th, 2013 -- Found an error in the calculation of H (was
% underestimating H, and therefore overestimating the depth average)


% Set defaults
dim = 1;
zind = 1:length(z);

if nargin >= 3
    depthlims = varargin{1};
    if ~isempty(depthlims)
        zind = find(z>=depthlims(1) & z<=depthlims(2));
    end
end

if nargin==4
    dim = varargin{2};
end

if nargin>4
    error('Too many inputs')
end

%Compute depth average
z = z(zind);
dz = z(2)-z(1);
H = length(z)*dz; %Assumes that all bins are the same size
%H = z(end)-z(1)

if dim == 1
    subfield = field(zind,:);
elseif dim == 2
    subfield = field(:,zind);
end

avgfield = 1/H*nansum(subfield,dim)*dz;
%avgfield = 1/(length(z))*nansum(subfield,dim)
zout = z;

end


