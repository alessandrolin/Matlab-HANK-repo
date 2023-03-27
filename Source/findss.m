%% findss.m
% Finds the steady state up to a specified tolerance
% Equilibrium specifies a beta for a given r

%created b mp 30/3/2022
function[ss,hist] = findss(guess,config,param,maxout_iter)
%   Simplify:
na = param.N_a;
nz = param.N_z;
%   Assign Guesses
guess_R     = guess.R;
guess_a     = guess.a;
guess_c     = guess.c;
guess_v     = guess.v;

%% Initialize Errors;
it1         = 0;
it1_hist    = [[1:1000]',nan(1000,4)]; % iteration, R, A
err1        = 1;
ss.Y        = 1;

%%  Start loop
while 1
    it1 = it1+1;
    errliq    = 1;
    it2     = 0;
    %% Asset Prices
    q = 1;
    %% Find Taxes
    T0 = param.bbb*(1/param.SPi-1/guess_R);
    % Inner loop: given beta and prices, obtain policy functions
    while errliq>config.tolV
        it2 = it2+1;
        [guess_a_t,guess_c,guess_v]  = egm_step(param,guess_R,guess_v,...
            param.SPi,1,T0,param.bet,0,ss);
        errliq = max(max(abs(guess_a-guess_a_t)));
        guess_a = guess_a_t;
    end
    %% One Asset
    if param.N_z ~=1 % HA
        L = trans_matrix(guess_a,param.a_g,param.P_z); % Obtain transition matrix
        D = (L')^10000; D = D(:,1);     % obtain stationary distribution
        ss.c = guess_c(:)'*D;
    else % RA
        L = eye(param.N_a);
        D = zeros(param.N_a,1);
        D(param.a_g ==0) = 1;
    end
    % Asset Market Clearing
    errliq = (D(:,1)'*repmat(param.a_g,param.N_z,1)-param.bbb); % excess savings
    % print results for current iteration
    fprintf('Outer iteration %d:\n',it1)
    %% Ad-hoc updates
    err1 = max(abs([(errliq)]));
    ss.err2 = (err1);
    % Display;
    fprintf('R  = %3.3f Excess savings in B: %3.8f %%.\n', [guess_R,100*errliq/(eps+param.bbb)])
    fprintf('beta = %3.6f .\n', [param.bet])

    % Quit if reached maxout iter (relevant for broyden)
    if it1>=maxout_iter
        break
    end
    if err1<config.hightol
        break
    end
    % Update beta - if too low savings (errliq<0), save more
    param.bet = param.bet*exp(-0.01/50*errliq).*0.9 + ...
        param.bet*0.1;

end

%%  Plot stationary distribution
figure('units','normalized','outerposition',[0.5 0 0.5 0.5])
for ii=1:param.N_z
    subplot(2,param.N_z+1,ii)
    plot(param.a_g,D((ii-1)*param.N_a+1:param.N_a*ii));hold on; xlabel('a-grid'); title(strcat('z=',num2str(param.z_g(ii))))
    subplot(2,param.N_z+1,ii+1+param.N_z)
    plot(D((ii-1)*param.N_a+1:param.N_a*ii));hold on; ylabel('grid')
end
subplot(2,param.N_z+1,param.N_z+1);
yyaxis left; plot(param.a_g,reshape(D,param.N_a,param.N_z)*ones(param.N_z,1));hold on;
yyaxis right; plot(param.a_g,tril(ones(param.N_a))*reshape(D,param.N_a,param.N_z)*ones(param.N_z,1));hold on; xlabel('a-grid'); title('Aggregate')
subplot(2,param.N_z+1,2*(param.N_z+1))
plot(reshape(D,param.N_a,param.N_z)*ones(param.N_z,1));hold on; ylabel('grid')
sgtitle('Stationary Asset Distribution')
sfbr    = 0;
for ii=1:param.N_z
    sfbr = sfbr + D((ii-1)*param.N_a+1);
end
fprintf('Stationary fraction at Borrowing Constraint: %g\t\n',sfbr);

%%  Allocate other SS variables
ss.Pi   = param.SPi;
ss.R    = guess_R;

ss.ga   = guess_a;
ss.gc   = guess_c;
ss.v    = guess_v;
ss.b    = param.bbb;
ss.t    = T0;
ss.D    = D;
ss.Y    = 1;

ss.bet  = param.bet;

ss.L = L;
%% Our Vector of endovars;
ss.X = [ss.Y,ss.Pi,ss.b,ss.t,0,ss.R];

%%  Save history of errors
hist.err    = err1;
hist.hist   = it1_hist;

%% Save errors themselves
ss.err1     = err1;
ss.errliq     = errliq;
ss.err2     = err1;

%% Ensure Aggregate Resource Constraint is met
ss.C   = ss.gc(:)'*D ;
ss.arc = ss.C-ss.Y;

b=1;
