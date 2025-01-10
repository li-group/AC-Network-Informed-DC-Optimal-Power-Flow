clc; clear all; close all;
warning('off','all')

system_name = 'pglib_opf_case118_ieee';
data = pglib_opf_case118_ieee;
matpower_data_name = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Matpower_CPF\\Generation_Matpower_%s.mat', system_name);

load(matpower_data_name);

DEMAND = Data1.DataMATPOWER.Demand_Probe;

nbranch = size(data.branch,1);
nbus = size(data.bus,1);

Pd = DEMAND(:,1:nbus)*100;
Qd = DEMAND(:,nbus+1:2*nbus)*100;

%% 
flagAC = zeros(186,1);
flagDC = zeros(186,1);
for l = 155:155
    fprintf('%d\n', l)
    nsamples = 1;
    mpc = data;
    nbus = size(mpc.bus,1);
    idx_gen_on = find(mpc.gen(:,8)); % Indices of generators that are online
    idx_P = find(mpc.bus(:,3)); % Indices for active power demand
    idx_Q = find(mpc.bus(:,4)); % Indices for reactive power demand
    mpc.branch(l,:) = []; % Remove the branch at index l
    
    % Initialize matrices for dispatch data
    DISPATCH_AC = zeros(nsamples, length(idx_gen_on));
    DISPATCH_DC = zeros(nsamples, length(idx_gen_on));
    
    iter = 276;
    while iter <= 276
        fprintf('%d\n', iter)    
        
        % Update demands for active and reactive power
        mpc.bus(:,3) = Pd(iter,:);
        mpc.bus(:,4) = Qd(iter,:);
        
        % Run AC and DC optimal power flow
        s_ac = runopf(mpc, mpoption('out.all', 0, 'opf.start', 3, ...
                                    'verbose', 0));
        s_dc = rundcopf(mpc, mpoption('out.all', 0, 'opf.start', 3, ...
                                      'verbose', 0));
        Pg_ac = s_ac.var.val.Pg; % Active power generation (AC)
        Pg_dc = s_dc.var.val.Pg; % Active power generation (DC)

        flagAC(l, iter) = s_ac.success; % Check if AC OPF succeeded
        % flagDC(l, iter) = s_dc.success; % Uncomment if needed
        
        DISPATCH_AC(iter,:) = Pg_ac';
        % DISPATCH_DC(iter,:) = Pg_dc'; % Uncomment if needed
        iter = iter + 1;
        
        % Create dynamic field name using sprintf
        field_name = sprintf('Pg_%d', l);
        
        % Save data into the structure with the dynamic field name
        Data1.DataMATPOWER.ContingenciaAC.(field_name) = DISPATCH_AC;
        % Data1.DataMATPOWER.ContingenciaDC.(field_name) = DISPATCH_DC; % Uncomment if needed
    end
end

%% 
% Concatenate the file name
output_file_name = sprintf('Generation_Matpower_%s.mat', system_name);

% Save the data into the .mat file
save(output_file_name, 'Data1');

%% 
% Example cost data processing for different test cases
data1 = pglib_opf_case30_ieee;
data2 = pglib_opf_case57_ieee;
data3 = pglib_opf_case118_ieee;
c1 = data1.gencost(:,6)*100;
c2 = data2.gencost(:,6)*100;
c3 = data3.gencost(:,6)*100;

BaseCost = Data1.DataMATPOWER.ContingenciaAC.Pg_94 * c3;

%% 

% Export data to CSV files (examples, uncomment as needed)
% writematrix(DEMAND, 'DEMAND_PROBE.csv');
% writematrix(DISPATCH_AC(:,1:7), 'DISPATCH_AC_PROBE.csv');
% writematrix(DISPATCH_DC, 'DISPATCH_DC_PROBE.csv');
% writematrix(COST_AC', 'COST_AC_PROBE.csv');
% writematrix(COST_DC', 'COST_DC_PROBE.csv');
% writematrix(X, 'X.csv');
% writematrix(Y, 'Y.csv');
% writematrix(Z, 'Z.csv');
