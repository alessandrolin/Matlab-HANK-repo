%% This function calculates the steady state REAL interest rate
%  There is a multilayers (2) loop
%   - update the real rate to clear the asset market
%       - calculate the policy function (with endogenous grid method)
%       - obtain the stationary distribution
%       - obtain the steady state demand for assets  
function [ss1,hist] = steady_state(param,config);
global guess

%%  Preallocate - Check if distribution exists
if isfile(strcat(config.root,'/Store/ss_file.mat'))
    load(strcat(config.root,'/Store/ss_file.mat'))
    guess.R     = ss1.R;
    guess.a     = ss1.ga;  
    param.bet   = ss1.bet;
    guess.c     = ss1.c;
    guess.v     = ss1.v; 
else
    guess.R     = 1.025^0.25;        
    guess.a     = param.a_g*0.95+param.z_g'*0.01;
    guess.c     = param.a_g/param.SPi+param.z_g'*1-param.bbb*(guess.R/param.SPi-1)-max(param.a_b,guess.a)/guess.R;
    guess.v     = exp(-param.sig*guess.c)/param.SPi;  
end

%% Temporary Tolerance assignments
oldtolss = config.tolss;
config.tolss = config.hightol; 

% note: apply your favorite solver to the function below
[ss1,hist] = findss(guess,config,param,1000);
close all
config.tolss = oldtolss;

% Save;
mkdir('Store')
save(strcat(config.root,'/Store/ss_file'),'ss1')
save(strcat(config.root,'/Store/param',date),'param')
 
%% Store
save(strcat(config.root,'/Store/ss_file'),'ss1')
save(strcat(config.root,'/Store/param',date),'param')
