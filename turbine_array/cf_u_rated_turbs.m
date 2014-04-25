function [turbines]=cf_u_rated_turbs(us,turb,u_rated)

[num_t, num_els]=size(us);

%num_rated=1000;
for jj=1:num_els
    jj
    u=us(:,jj);
    max_u=max(u,[],1);
    for ii=1:length(u_rated)
        turbs(ii)=turb;
        turbs(ii).rated=u_rated(ii);
        turbs(ii).Prated=turb.Cp*(1/2*turbs(ii).rho*turbs(ii).area*turbs(ii).rated^3);
        P=calculate_power_RK(u,u,turbs(ii));
        turbs(ii).meanP = mean(P);
        turbs(ii).cf=turbs(ii).meanP/turbs(ii).Prated;
    end
    a=find([turbs.cf]>turb.cf);
    if length(a)>0
    [ur uri]=max([turbs(a).rated]);
    turbines(jj)=turbs(a(uri));
    else
    turbines(jj)=turbs(1);        
    end
        
end

end
