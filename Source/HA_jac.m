function [J] = HA_jac(dY,dPi,dR,dt,dbe,ss,par,T,epsv)

fprintf('Getting Jacobian HA - shocks - dY:%g\tdPi: %g\tdR: %g\tdt: %g\tdbe:  %g\t',dY,dPi,dR,dt,dbe);

for pm = ["p" "m"]
    if pm == "m"
        dY  = -dY;
        dPi = -dPi;
        dR  = -dR;
        dt  = -dt;
        dbe = -dbe;
    end
    % setup shocks
    Y   = ones(T,1);
    Pi  = ones(T,1)*ss.Pi;
    R   = ones(T,1)*ss.R;
    t   = ones(T,1)*ss.t;
    be  = ones(T,1)*par.bet;

    Y(end)  = Y(end)+dY;
    Pi(end) = Pi(end)+dPi;
    R(end)  = R(end)+dR;
    t(end)  = t(end)+dt;
    be(end) = be(end)+dbe;

    % preallocate matrices
    ga      = nan(par.N_a,par.N_z,T);
    gc      = nan(par.N_a,par.N_z,T);
    v       = nan(par.N_a,par.N_z,T);

% tic
    % Preallocate Cells containing indeces;
    V   = ss.v;

    for ind_t = T:-1:1
        [ga(:,:,ind_t),~,V]  = egm_step(par,R(ind_t),...
            V,...
            Pi(ind_t),Y(ind_t),t(ind_t),be(ind_t),0);
    end
    
    % Compute dD
    dD = NaN(par.N_a*par.N_z,T);
    for ind_t=1:T
       [~,indeces,values]    = trans_matrix(ga(:,:,ind_t),par.a_g,par.P_z,1);
       subsparser   = horzcat(indeces,values);
       L = sparse(subsparser(:,1),subsparser(:,2),subsparser(:,3),par.N_a*par.N_z,par.N_a*par.N_z); % generate sparse matrix
       dD(:,ind_t) = L'*ss.D; % would included -ss.D but cancels out in 66
    end
    
    
    eval(strcat("ga_",pm,"=ga;"))
    eval(strcat("dD_",pm,"=dD;"))
    clearvars ga gv dD v    
    
end

dy = (ga_p-ga_m)/2;
clearvars ga* gc* v*
dD = (dD_p-dD_m)/2;


% % fake news algorithm
F   = zeros(T);
ssga = reshape(ss.ga,par.N_a*par.N_z,1);
for is = 1:T
    F(1,is)     = reshape(dy(:,:,T-is+1),par.N_a*par.N_z,1)'*ss.D;
    F(2,is)     = ssga'*dD(:,T-is+1);
end

%% Build the remainder of F
for is = 1:T
    F(3:end,is) = epsv(1:end-2,:)*dD(:,T-is+1);
end

J = nan(T,T);
JJ = nan(T,T*2-1);
for it = 1:T
    JJ(it,:) = [zeros(1,T-it) F(it,:) zeros(1,it-1)];
end
for it = 2:T
    JJ(it,:) = JJ(it,:)+JJ(it-1,:);
end
for it = 1:T
    J(it,:) = JJ(it,T-it+1:2*T-it);
end






