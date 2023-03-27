function [V,D,L,ga] = obtain_micro_pf(dX,dbe,ss,par,T,D)

% obtain variables
Y   = ss.Y+dX(:,1);
Pi  = ss.Pi+dX(:,2);
t   = ss.t+dX(:,4);
R   = ss.R+dX(:,6);
be  = par.bet+dbe;


VV  = nan(par.N_a*par.N_z,T);
V   = ss.v;
ga  = nan(par.N_a*par.N_z,T);
ga_aux = nan(par.N_a,par.N_z,T);

for ind_t = T:-1:1
    [ga_aux(:,:,ind_t),~,V]  = egm_step(par,R(ind_t),...
        V,...
        Pi(ind_t),Y(ind_t),t(ind_t),be(ind_t),0,ss);
    VV(:,ind_t)     = reshape(V,par.N_a*par.N_z,1);
    ga(:,ind_t)     = reshape(ga_aux(:,:,ind_t),par.N_a*par.N_z,1);
end



% Compute the Distribution
DD = NaN(par.N_a*par.N_z,T);



for ind_t=1:T
    [~,indeces,values]    = trans_matrix(ga_aux(:,:,ind_t),par.a_g,par.P_z,1);
    subsparser   = horzcat(indeces,values);
    L = sparse(subsparser(:,1),subsparser(:,2),subsparser(:,3),par.N_a*par.N_z,par.N_a*par.N_z); % generate sparse matrix
    % Simulate
    if ind_t == 1
        DD(:,1) = L'*D; % D is not ss.D, can differ
    else
        DD(:,ind_t) = L'*DD(:,ind_t-1);
    end
end

V  = VV';
ga = ga';
D = DD';






