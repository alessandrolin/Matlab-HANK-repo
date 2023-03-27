%%  HANK - 1 asset
%   by Alessandro Lin & Marcel Peruffo
%   This version: 22 Mar 2023

%%  Structure
%   0) Matlab Configuration
%   1) Parameters declaration
%   2) Steady State finding
%   3) Setup Jacobians
%   4) Specify shock
%   5) Run Simulations

%%  0) Matlab Configuration
clear;                      clc;                        close all;
config.root = cd;                    
addpath([config.root '/Source']);    
%delete('Store/ss_file.mat')                         % comment if no change in parameters
%try;delete(strcat('Store/param',date,'.mat'));end   % comment if no change in parameters

%%  1) Parameters declaration
declare_parameters
close

%%  2) Steady State finding
[ss,hist]    = steady_state(param,config);
ss1 = ss;
param.bet    = ss1.bet; 
close all

%%  2.1) Compute MPCs
compute_MPCs 

%%  3) Setup Jacobians  
[J_x_pf,J_be_pf] = setup_jacobian(param,ss,config,param.TTT);
save(strcat(config.root,'/Store/J_x_pf',date),'J_x_pf')
save(strcat(config.root,'/Store/J_be_pf',date),'J_be_pf')

close all

%%  4) Specify Exogenous Shock values
load_exo_shocks 

%%  5) Simulation and plots
run_simulations 
create_plots;
 
