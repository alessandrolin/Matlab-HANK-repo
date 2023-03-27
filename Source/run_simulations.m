% Pre-invert the Jacobian w.r.t prices, for speed;
invJ  = inv(J_x_pf);

%% PF counterfactual 
% Equilibrium Deviations:
shock = [zeros(param.TTT*(param.nX),1)]+...
            J_be_pf*dbe_cf;
        
X_PF = -invJ*shock; 

X_PF = reshape(X_PF,param.TTT,param.nX);
[V_PF,D_PF,~,GA_PF] = obtain_micro_pf(X_PF,dbe_cf,ss1,param,param.TTT,1,ss1.D); % obtain the value function
% Set Xs up
  
clearvars a b RES err i ind_t D_new it it_1 GA_PF
