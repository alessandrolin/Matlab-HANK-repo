%%  Preferences
param.bet      = 0.9805;     % discount factor
param.sig      = 1.5;       % intertemporal elas. substitution
param.frisch   = 1;

%%  Phillips curve
param.kap    = 0.01; 

%%  Income Shocks 
param.z_g = [0.492297350129250	1	2.03129267248230]'; % z-grid
param.N_z = size(param.z_g,1);                          % size of z-grid 
param.P_z  = [0.965954808900000	0.0337503822000000	0.000294808900000001
    0.0168751911000000	0.966249617800000	0.0168751911000000
    0.000294808900000001	0.0337503822000000	0.965954808900000]; % z-transition matrix

param.Dzb = param.P_z'^100000;
param.Dzb = param.Dzb(:,1);                             % stationary D_z
param.z_g = param.z_g./(param.z_g'*param.Dzb);          % adjust for average = 1


%% Supply of Assets
param.bbb       = 1;    % Bss/Yss 25 of yearly GDP (3 times quarterly GDP)

%%  Assets
param.a_b = -1;  
param.amax   = (2.5);                               % exp of that
param.N_a = 101;                                    % size of a-grid
param.a_g = logspace(-4,param.amax,param.N_a)';     % create log-spaced a-grid
param.a_g = param.a_g-min(param.a_g)+param.a_b;     % shift for borrowing constraint
[~,idx] = min(abs(param.a_g));                      % Re set a value to be equal to 0
param.a_g(idx) = 0;
clear idx

figure;
subplot(1,2,1); plot(1:param.N_a,param.a_g,1:param.N_a,zeros(1,param.N_a))
title('Grid Points','Fontsize',14)
subplot(1,2,2); plot(param.a_g,zeros(param.N_a,1),'*'); xline(0)
title('Grid Spacing','Fontsize',14)
close

%%  Monetary policy
param.php   = 1.5;                  % monetary policy reaction
param.phy   = 0.125;                % monetary policy reaction
param.phr   = 0;                    % monetary policy reaction
param.SPi   = 1.02^0.25;            % steady state inflation

%%  Aggregate params
param.nX    = 6;                    % number of aggregate variables
param.TTT   = 300;                  % T truncation

%%  Algorithm Tolerance
config.tolss    = 1e-10;            % precise
config.hightol	= 1e-12;             
config.tolV	    = 1e-12;            % precise
