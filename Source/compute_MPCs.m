%% compute_MPCs.m

dy = 0.01;
T = 50;

%% One Asset 
% Compute Policy at t=1
[ga,c,~]  = egm_step(param,ss.R,...
    ss.v,...
    ss.Pi,1,ss.t,param.bet,dy,ss);

DD = NaN(param.N_a*param.N_z,T);

[~,indeces,values]    = trans_matrix(ga,param.a_g,param.P_z,1);
subsparser   = horzcat(indeces,values);
L = sparse(subsparser(:,1),subsparser(:,2),subsparser(:,3),param.N_a*param.N_z,param.N_a*param.N_z); % generate sparse matrix
for ind_t=1:T
    % Simulate
    if ind_t == 1
        DD(:,1) = L'*ss.D;  
    else
        DD(:,ind_t) = ss.L'*DD(:,ind_t-1);
    end
end 
    
    
%% Consumption Function
MPC = NaN(1,T+1);
for ind_t=1:T
    % Simulate
    if ind_t == 1
        MPC(ind_t) = (c(:)'*ss1.D-ss1.c)/dy; 
    else
         MPC(ind_t) = (ss1.gc(:)'*DD(:,ind_t-1)-ss1.c)/dy;
    end
end

fprintf('1-year MPC is %3.2f \n', sum(MPC(1:4)))

clear c ga dy ind_t T

