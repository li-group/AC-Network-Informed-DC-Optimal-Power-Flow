clc; clear all; close all;
warning('off','all')
% system_name = 'pglib_opf_case2742_goc';
% system_name = 'pglib_opf_case14_ieee';
% system_name = 'pglib_opf_case30_ieee';
% system_name = 'pglib_opf_case118_ieee';
% system_name = 'pglib_opf_case118_ieeeq';
% system_name = 'pglib_opf_case300_ieee';
system_name = "pglib_opf_case1354_pegase";
% Concatenate the file name
file_name = sprintf('Generation_Matpower_%s.mat', system_name);
% file_name = sprintf('Load.mat');
% mpc0 = pglib_opf_case57_ieee;
% mpc0 = pglib_opf_case14_ieee;
% mpc0 = pglib_opf_case30_ieee;
% mpc0 =  pglib_opf_case118_ieee;

%% 
mpc0 =  pglib_opf_case1354_pegase;
mpci = ext2int(mpc0);
% mpci.branch(94,:)=[]; % Uncomment if you want to remove specific branches
[BBUS, BF, PBUSINJ, PFINJ] = makeBdc(mpci);

BBUS = full(BBUS);
BF = full(BF);
PBUSINJ = PBUSINJ/mpc0.baseMVA;
PFINJ = PFINJ/mpc0.baseMVA;
% save("MatrixBDC.mat","PFINJ","PBUSINJ","BF","BBUS")

%% 
load(file_name)
DataMATPOWER.BBUS = BBUS;
DataMATPOWER.BF = BF;
DataMATPOWER.PBUSINJ = PBUSINJ;
DataMATPOWER.PFINJ = PFINJ;

% RandDemand.BBUS = BBUS;
% RandDemand.BF = BF;
% RandDemand.PBUSINJ = PBUSINJ;
% RandDemand.PFINJ = PFINJ;

%% 

% save("Load.mat",'RandDemand', '-v7.3')
% Save the data in the .mat file
save(file_name, 'DataMATPOWER');

%% 
mpc0 = pglib_opf_case57_ieee;
s_ac = rundcopf(mpc0, mpoption('out.all', 0, 'verbose', 0));
DUAL = s_ac.bus(:,14);
Pg = s_ac.gen(:,2);
Cost = s_ac.gencost(:,6);

A = DUAL(s_ac.gen(:,1)).*Pg;
B = Pg.*Cost;
A - B
