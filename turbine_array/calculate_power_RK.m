%Brian Polagye
%March 6, 2013

%Description: calculate power output time series based on input velocity and turbine
%   parameters

function [P] = calculate_power2(u, t, turbine)

debug = 0;

%operate on speed, not velocity
u = abs(u);

%initialize turbine powre density
P = zeros(size(u));

[nt,nn]=size(u);

turbine.cutin=ones(nt,1)*turbine.cutin;
turbine.rated=ones(nt,1)*turbine.rated;
turbine.Prated=ones(nt,1)*turbine.Prated;
turbine.cutoff=ones(nt,1)*turbine.cutoff;
    
    %calculate power time series - assuming constant efficiency
    P= ((u-turbine.cutin)./(turbine.rated-turbine.cutin)).^3.*turbine.Prated;
    
    
    P(u<turbine.cutin) = 0;
    P(u>turbine.rated) = turbine.Prated(u>turbine.rated);
    P(u>turbine.cutoff) = 0;
    

if debug
    
    %diagnostic plot - turbine power density
    %   plotted in customary units of kW/m2
    figure(1)
    clf
    for i = 1:length(turbine)
        hold all
        plot(t-t(1),P/1000,'-')
    end
    xlabel('Time (days)','fontweight','b')
    ylabel('Turbine power density (kW/m^2)','fontweight','b')
    
end

end