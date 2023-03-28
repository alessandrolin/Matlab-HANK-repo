function [A0,c0,v0] = egm_step(par,R0,F1,Pi0,Y0,T0,bet,dy)

z_g = par.z_g;  

c0  = max(0,(bet*(R0*F1*par.P_z')).^(1/(-par.sig)));

A0  = (par.a_g/R0+c0-z_g'*(Y0-T0)-dy)*Pi0;

for ii = 1:par.N_z
    A0(:,ii) = max(par.a_b,nakeinterp1(A0(:,ii),par.a_g,par.a_g));
end

c0  = dy+par.a_g/Pi0+z_g'*(Y0-T0)-A0/R0;
v0  = c0.^(-par.sig)/Pi0;


