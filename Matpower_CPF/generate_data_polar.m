clc; clear all; close all;
warning('off','all')
%% 

mpc0 = pglib_opf_case118_ieee;
mpc = mpc0;
mpc.gencost(:,5) = mpc.gencost(:,6)/1000;

nsamples = 100;

Pd = mpc0.bus(:,3)/100;
Qd = mpc0.bus(:,4)/100;
nbus = size(mpc.bus,1);
idx_gen_on = find(mpc.gen(:,8)); % Indices of generators that are online
idx_P = find(mpc.bus(:,3)); % Indices for active power demand
idx_Q = find(mpc.bus(:,4)); % Indices for reactive power demand

DEMAND = zeros(nsamples, 2 * nbus);
DISPATCH_AC = zeros(nsamples, 2 * length(idx_gen_on));
DISPATCH_DC = zeros(nsamples, length(idx_gen_on));

correlation = 0.7; % Demand correlation 0.7
Sig = eye(nbus);
Sig(Sig==0) = correlation; 
mu = ones(nbus,1); % Mean of the distribution
lb = -0.3*ones(nbus,1); % Lower bound
ub = 0.3*ones(nbus,1); % Upper bound

random_scale = mvrandn_andre(lb, ub, Sig, nsamples);
random_scale = mu + random_scale;

histogram(random_scale)
correlation = 0.7; % Demand correlation 0.7
Sig = eye(nbus);
Sig(Sig==0) = correlation; 
mu = ones(nbus,1); % Mean of the distribution
lb = -0.15*ones(nbus,1); % Lower bound
ub = 0*ones(nbus,1); % Upper bound

random_Q = mvrandn_andre(lb, ub, Sig, nsamples);
random_Q = mu + random_Q;

P = (mpc.bus(:,3).*random_scale)';
Q = (mpc.bus(:,4).*random_Q)';
%% 

RandDemand = struct();

RandDemand.P = P;
RandDemand.Q = Q;
save('filename.mat', 'RandDemand', '-v7.3');
%% 

load("RandPQ.mat")
%% 

data = pglib_opf_case118_ieee;
mpc = mpc0;

nbranch = size(data.branch,1);
nbus = size(data.bus,1); 
ngen = size(data.gen,1);

data.gencost(:,5) = data.gencost(:,6)/1000;

% Initialization of variables
nsamples = 4000;
timeAC = zeros(1, nsamples);
timeDC = zeros(1, nsamples);
flag = zeros(1, nsamples);
DEMAND = zeros(nsamples, 2 * size(data.bus, 1));
DISPATCH_AC = zeros(nsamples, 2 * ngen);
DISPATCH_DC = zeros(nsamples, ngen);
COST_AC = zeros(1, nsamples);
COST_DC = zeros(1, nsamples);
DualAC = zeros(nsamples, size(data.bus, 1));
DualDC = zeros(nsamples, size(data.bus, 1));

% Iteration over samples
iter = 1;
while iter <= nsamples
    fprintf('%d\n', iter)
    
    mpc.bus(:, 3) = P(iter, :); % Assign active power demand
    mpc.bus(:, 4) = Q(iter, :); % Assign reactive power demand

    tic;
    s_ac = runopf(mpc, mpoption('opf.start', 3, 'out.all', 0, 'verbose', 0));
    timeAC(iter) = toc; 

    tic; % Start measuring execution time of rundcopf
    s_dc = rundcopf(mpc, mpoption('opf.start', 3, 'out.all', 0, 'verbose', 0));
    timeDC(iter) = toc; % Measure execution time of rundcopf

    if s_ac.success == 1 && s_dc.success == 1
        cost_ac = s_ac.f;
        Pg_ac = s_ac.var.val.Pg;
        Qg_ac = s_ac.var.val.Qg;
        cost_dc = s_dc.f;
        Pg_dc = s_dc.var.val.Pg;
        flag(iter) = s_ac.success;
        Pd = P(iter, :);
        Qd = Q(iter, :);
        DEMAND(iter, :) = [Pd, Qd] / mpc0.baseMVA;
        DISPATCH_AC(iter, :) = [Pg_ac', Qg_ac'];
        DISPATCH_DC(iter, :) = Pg_dc';
        COST_AC(iter) = cost_ac;
        COST_DC(iter) = cost_dc;
        DualAC(iter, :) = s_ac.bus(:, 14);
        DualDC(iter, :) = s_dc.bus(:, 14);
    else
        flag(iter) = 0;
        disp("No convergence")
    end

    iter = iter + 1;
end

%% 

figure
plot(COST_AC, COST_DC, 'o')
hold on
plot(COST_AC, COST_AC)
refline(1, 0)

%% 
indices_eliminar = []; % Indices to remove
Limite = mpc.gen(:,9)'/100;
for i = 1:nsamples
    for k = 1:7
        if DISPATCH_AC(i, k) > Limite(1, k)
            indices_eliminar = [indices_eliminar, i];
        end
    end
end
indices_eliminar = unique(indices_eliminar);
DISPATCH_AC(indices_eliminar, :) = [];
DISPATCH_DC(indices_eliminar, :) = [];
DEMAND(indices_eliminar, :) = [];
COST_AC(:, indices_eliminar) = [];
COST_DC(:, indices_eliminar) = [];

%% 
DEMAND = DEMAND(flag == 1, :);
COST_AC = COST_AC(flag == 1);
COST_DC = COST_DC(flag == 1);
DISPATCH_DC = DISPATCH_DC(flag == 1, :);
DISPATCH_AC = DISPATCH_AC(flag == 1, :);
DualAC = DualAC(flag == 1, :);
DualDC = DualDC(flag == 1, :);
timeAC = timeAC(flag == 1);
timeDC = timeDC(flag == 1);

%% 
DataMATPOWER = struct();
DataMATPOWER.DualAC = DualAC;
DataMATPOWER.DualDC = DualDC;
DataMATPOWER.timeAC = timeAC;
DataMATPOWER.timeDC = timeDC;
DataMATPOWER.CostDC = COST_DC';
DataMATPOWER.CostAC = COST_AC';
DataMATPOWER.DispatchAC = DISPATCH_AC;
DataMATPOWER.DispatchDC = DISPATCH_DC;
DataMATPOWER.Demand = DEMAND;

system_name = 'pglib_opf_case118_ieee';

% Concatenate the file name
file_name = sprintf('Generation_Matpower_%s.mat', system_name);

% Save the data in the .mat file
save(file_name, 'DataMATPOWER');

