clear all ,clc

% Contingencias_DC=load('Contingencias_DC.mat')
% Contingencias_PDC=load('Contingencias_PDC.mat')
% Contingencias_Real=load('Contingencias_Real.mat')
% Contingencias_DC=load('Contingencias_DC.mat');

% Define los nombres de los sistemas
system_names = {'pglib_opf_case30_ieee','pglib_opf_case57_ieee','pglib_opf_case118_ieee'};
num_systems = numel(system_names);

% Define los nombres de las variables
variable_names = {'Matpower_Data', 'Generator_Data', 'pDC_Model', 'DC_Model','Contingencias_AC','Contingencias_DC','Contingencias_pDC','DC_Model_unconstrain','pDC_Model_unconstrain',};

% Crea una estructura para almacenar los datos
all_data = struct();

% Rellena la estructura con los datos de cada sistema
for i = 1:num_systems
    system_name = system_names{i};
    nombre_datos_matpower = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Matpower_CPF\\Generation_Matpower_%s.mat', system_name);
    nombre_datos_mpec = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\NeuralNetwork\\Generation_MPCE_%s.mat', system_name);
    nombre_datos_pDCOPF = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos\\pDC_Model_%s.mat', system_name);
    nombre_datos_DCOPF = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos\\DC_Model_%s.mat', system_name);
    nombre_datos_ContigenciasAC = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos\\Contingencias_AC_%s.mat', system_name);
    nombre_datos_ContigenciasDC = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos\\Contingencias_DC_%s.mat', system_name);
    nombre_datos_ContigenciaspDC = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos\\Contingencias_pDC_%s.mat', system_name);
    nombre_datos_pDCOPF_unconstrain = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos\\pDC_Model_unconstrain_%s.mat', system_name);
    nombre_datos_DCOPF_unconstrain = sprintf('C:\\Users\\andre\\OneDrive\\Escritorio\\CodeDepurado300\\Modelos\\DC_Model_unconstrain_%s.mat', system_name);
    

    all_data(i).Matpower_Data = load(nombre_datos_matpower);
    all_data(i).Generation_Data = load(nombre_datos_mpec);
    all_data(i).ContingenciaAC_Data = load(nombre_datos_ContigenciasAC);
    all_data(i).ContingenciaDC_Data = load(nombre_datos_ContigenciasDC);
    all_data(i).ContingenciapDC_Data = load(nombre_datos_ContigenciaspDC);
    all_data(i).pDC_Model = load(nombre_datos_pDCOPF);
    all_data(i).DC_Model = load(nombre_datos_DCOPF);
    all_data(i).pDC_Model_unconstrain = load(nombre_datos_pDCOPF_unconstrain);
    all_data(i).DC_Model_unconstrain = load(nombre_datos_DCOPF_unconstrain);
end

% Accede a los datos de un sistema específico
% system_index = 1; % Cambia este índice según el sistema que desees acceder
% Matpower_Data = all_data(system_index).Matpower_Data;
% Generator_Data = all_data(system_index).Generator_Data;
% pDC_Model = all_data(system_index).pDC_Model;
% DC_Model = all_data(system_index).DC_Model;
% Feasibility = all_data(system_index).Feasibility;


set(groot, 'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultAxesFontSize',14);
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultLegendFontSize',16);
set(groot,'defaultLineLineWidth',1);

%% 

scatter(all_data(2).Matpower_Data.DataMATPOWER.DualAC_Probe,all_data(2).DC_Model.Dual/100)
hold on 
%scatter(all_data(2).Matpower_Data.DataMATPOWER.DualAC_Probe,all_data(2).pDC_Model.Dual/100)
refline(1,0)
%% TABLES ZONE

%Tablas para error PG 
experimento1Pg = abs(all_data(1).pDC_Model.Pg-all_data(1).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:6));
experimento2Pg = abs(all_data(2).pDC_Model.Pg-all_data(2).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:7));
experimento3Pg = abs(all_data(3).pDC_Model.Pg-all_data(3).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:54));
experimento4Pg = abs(all_data(1).DC_Model.Pg-all_data(1).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:6));
experimento5Pg = abs(all_data(2).DC_Model.Pg-all_data(2).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:7));
experimento6Pg = abs(all_data(3).DC_Model.Pg-all_data(3).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:54));

experimento1Pg = reshape(experimento1Pg, [], 1);
experimento2Pg = reshape(experimento2Pg, [], 1);
experimento3Pg = reshape(experimento3Pg, [], 1);
experimento4Pg = reshape(experimento4Pg, [], 1);
experimento5Pg = reshape(experimento5Pg, [], 1);
experimento6Pg = reshape(experimento6Pg, [], 1);




% % Calcular mínimo, máximo y media de cada experimento
% minimo_exp1 = sprintf('%.4f',min(min(experimento1Pg)));
% minimo_exp2 = sprintf('%.4f',min(min(experimento2Pg)));
% minimo_exp3 = sprintf('%.4f',min(min(experimento3Pg)));
% minimo_exp4 = sprintf('%.4f',min(min(experimento4Pg)));
% minimo_exp5 = sprintf('%.4f',min(min(experimento5Pg)));
% minimo_exp6 = sprintf('%.4f',min(min(experimento6Pg)));
% 
% sprintf('%.4f',mean(mean(experimento1Pg)));
% media_exp1 = sprintf('%.4f',mean(mean(experimento1Pg)));
% media_exp2 = sprintf('%.4f',mean(mean(experimento2Pg)));
% media_exp3 = sprintf('%.4f',mean(mean(experimento3Pg)));
% media_exp4 = sprintf('%.4f',mean(mean(experimento4Pg)));
% media_exp5 = sprintf('%.4f',mean(mean(experimento5Pg)));
% media_exp6 = sprintf('%.4f',mean(mean(experimento6Pg)));
% 
% maximo_exp1 = sprintf('%.4f',max(max(experimento1Pg)));
% maximo_exp2 = sprintf('%.4f',max(max(experimento2Pg)));
% maximo_exp3 = sprintf('%.4f',max(max(experimento3Pg)));
% maximo_exp4 = sprintf('%.4f',max(max(experimento4Pg)));
% maximo_exp5 = sprintf('%.4f',max(max(experimento5Pg)));
% maximo_exp6 = sprintf('%.4f',max(max(experimento6Pg)));


Analysis_Pg = table([minimo_exp1; minimo_exp4; media_exp1;media_exp4; maximo_exp1; maximo_exp4], ...
                    [minimo_exp2; minimo_exp5; media_exp2;media_exp5; maximo_exp2; maximo_exp5], ...
                    [minimo_exp3; minimo_exp6; media_exp3;media_exp6; maximo_exp3; maximo_exp6], ...
                    'RowNames', {'Minimum pDC', 'Minimum DC','Mean pDC','Mean DC','Maximum pDC', 'Maximum DC'}, ...
                    'VariableNames', {'30 IEEE bus', '57 IEEE bus', '118 IEEE bus'});


%% 

%Tablas para error Costos 
experimento1Cost = all_data(1).pDC_Model.ErrorCost.Absolute_Error;
experimento2Cost = all_data(2).pDC_Model.ErrorCost.Absolute_Error;
experimento3Cost = all_data(3).pDC_Model.ErrorCost.Absolute_Error;
experimento4Cost = all_data(1).DC_Model.ErrorCost.Absolute_Error;
experimento5Cost = all_data(2).DC_Model.ErrorCost.Absolute_Error;
experimento6Cost = all_data(3).DC_Model.ErrorCost.Absolute_Error;
% Calcular mínimo, máximo y media de cada experimento
minimo_exp1 = sprintf('%.4e',min(experimento1Cost));
minimo_exp2 = sprintf('%.4e',min(experimento2Cost));
minimo_exp3 = sprintf('%.4e',min(experimento3Cost));
minimo_exp4 = sprintf('%.4e',min(experimento4Cost));
minimo_exp5 = sprintf('%.4e',min(experimento5Cost));
minimo_exp6 = sprintf('%.4e',min(experimento6Cost));

%minimo_exp3 = sprintf('%.4e',min(experimento3Cost));

maximo_exp1 = sprintf('%.4e',max(experimento1Cost));
maximo_exp2 = sprintf('%.4e',max(experimento2Cost));
maximo_exp3 = sprintf('%.4e',max(experimento3Cost));
maximo_exp4 = sprintf('%.4e',max(experimento4Cost));
maximo_exp5 = sprintf('%.4e',max(experimento5Cost));
maximo_exp6 = sprintf('%.4e',max(experimento6Cost));

%maximo_exp3 = sprintf('%.4e',max(experimento3Cost));

media_exp1 = sprintf('%.4e',mean(experimento1Cost));
media_exp2 = sprintf('%.4e',mean(experimento2Cost));
media_exp3 = sprintf('%.4e',mean(experimento3Cost));
media_exp4 = sprintf('%.4e',mean(experimento4Cost));
media_exp5 = sprintf('%.4e',mean(experimento5Cost));
media_exp6 = sprintf('%.4e',mean(experimento6Cost));

%media_exp3 = sprintf('%.4e',mean(experimento3Cost));

Analysis_Cost = table([minimo_exp1; minimo_exp4; media_exp1; media_exp4; maximo_exp1; maximo_exp4], ...
                      [minimo_exp2; minimo_exp5; media_exp2; media_exp5; maximo_exp2; maximo_exp5], ...
                      [minimo_exp3; minimo_exp6; media_exp3; media_exp6; maximo_exp3; maximo_exp6], ...
                      'RowNames', {'Minimum pDC', 'Minimum DC', 'Mean pDC','Mean DC', 'Maximum pDC', 'Maximum DC'}, ...
                      'VariableNames', {'30 IEEE bus', '57 IEEE bus', '118 IEEE bus'});

%% ANALISIS REVENUE ADEQUENCY

%Tablas 
experimento1RA = all_data(1).pDC_Model.RA_Result;
experimento2RA = all_data(2).pDC_Model.RA_Result;
experimento3RA = all_data(3).pDC_Model.RA_Result;
experimento4RA = all_data(1).pDC_Model_unconstrain.RA_Result;
experimento5RA = all_data(2).pDC_Model_unconstrain.RA_Result;
experimento6RA = all_data(3).pDC_Model_unconstrain.RA_Result;
% Calcular mínimo, máximo y media de cada experimento
minimo_exp1 = round(min(experimento1RA),4);
minimo_exp2 = round(min(experimento2RA),4);
minimo_exp3 = round(min(experimento3RA),4);
minimo_exp4 = round(min(experimento4RA),4);
minimo_exp5 = round(min(experimento5RA),4);
minimo_exp6 = round(min(experimento6RA),4);

%minimo_exp3 = round(min(experimento3RA));

maximo_exp1 = round(max(experimento1RA),4);
maximo_exp2 = round(max(experimento2RA),4);
maximo_exp3 = round(max(experimento3RA),4);
maximo_exp4 = round(max(experimento4RA),4);
maximo_exp5 = round(max(experimento5RA),4);
maximo_exp6 = round(max(experimento6RA),4);

%maximo_exp3 = round(max(experimento3RA));

media_exp1 = round(mean(experimento1RA),4);
media_exp2 = round(mean(experimento2RA),4);
media_exp3 = round(mean(experimento3RA),4);
media_exp4 = round(mean(experimento4RA),4);
media_exp5 = round(mean(experimento5RA),4);
media_exp6 = round(mean(experimento6RA),4);

%media_exp3 = sprintf('%.4e',mean(experimento3RA));

Analysis_RA = table([minimo_exp1; minimo_exp4; media_exp1; media_exp4; maximo_exp1; maximo_exp4], ...
                      [minimo_exp2; minimo_exp5; media_exp2; media_exp5; maximo_exp2; maximo_exp5], ...
                      [minimo_exp3; minimo_exp6; media_exp3; media_exp6; maximo_exp3; maximo_exp6], ...
                      'RowNames', {'Minimum pDC', 'Minimum DC', 'Mean pDC','Mean DC', 'Maximum pDC', 'Maximum DC'}, ...
                      'VariableNames', {'30 IEEE bus', '57 IEEE bus', '118 IEEE bus'});
experimento1RA = size(all_data(1).pDC_Model.RA_OK)
experimento2RA = size(all_data(2).pDC_Model.RA_OK)
experimento3RA = size(all_data(3).pDC_Model.RA_OK)
experimento4RA = size(all_data(1).pDC_Model_unconstrain.RA_OK)
experimento5RA = size(all_data(2).pDC_Model_unconstrain.RA_OK)
experimento6RA = size(all_data(3).pDC_Model_unconstrain.RA_OK)
%% ANALISIS COST RECOVERY

%Tablas 
experimento1CR = all_data(1).pDC_Model.CR_Result;
experimento2CR = all_data(2).pDC_Model.CR_Result;
experimento3CR = all_data(3).pDC_Model.CR_Result;
experimento4CR = all_data(1).pDC_Model_unconstrain.CR_Result;
experimento5CR = all_data(2).pDC_Model_unconstrain.CR_Result;
experimento6CR = all_data(3).pDC_Model_unconstrain.CR_Result;
% Calcular mínimo, máximo y media de cada experimento

experimento1CR = reshape(experimento1CR, [], 1);
experimento2CR = reshape(experimento2CR, [], 1);
experimento3CR = reshape(experimento3CR, [], 1);
experimento4CR = reshape(experimento4CR, [], 1);
experimento5CR = reshape(experimento5CR, [], 1);
experimento6CR = reshape(experimento6CR, [], 1);


minimo_exp1 = round(min(experimento1CR),4);
minimo_exp2 = round(min(experimento2CR),4);
minimo_exp3 = round(min(experimento3CR),4);
minimo_exp4 = round(min(experimento4CR),4);
minimo_exp5 = round(min(experimento5CR),4);
minimo_exp6 = round(min(experimento6CR),4);

%minimo_exp3 = round(min(experimento3CR));

maximo_exp1 = round(max(experimento1CR),4);
maximo_exp2 = round(max(experimento2CR),4);
maximo_exp3 = round(max(experimento3CR),4);
maximo_exp4 = round(max(experimento4CR),4);
maximo_exp5 = round(max(experimento5CR),4);
maximo_exp6 = round(max(experimento6CR),4);

%maximo_exp3 = round(max(experimento3CR));

media_exp1 = round(mean(experimento1CR),4);
media_exp2 = round(mean(experimento2CR),4);
media_exp3 = round(mean(experimento3CR),4);
media_exp4 = round(mean(experimento4CR),4);
media_exp5 = round(mean(experimento5CR),4);
media_exp6 = round(mean(experimento6CR),4);

%media_exp3 = sprintf('%.4e',mean(experimento3RA));

Analysis_CR = table([minimo_exp1; minimo_exp4; media_exp1; media_exp4; maximo_exp1; maximo_exp4], ...
                      [minimo_exp2; minimo_exp5; media_exp2; media_exp5; maximo_exp2; maximo_exp5], ...
                      [minimo_exp3; minimo_exp6; media_exp3; media_exp6; maximo_exp3; maximo_exp6], ...
                      'RowNames', {'Minimum pDC', 'Minimum DC', 'Mean pDC','Mean DC', 'Maximum pDC', 'Maximum DC'}, ...
                      'VariableNames', {'30 IEEE bus', '57 IEEE bus', '118 IEEE bus'});
experimento1CR = size(all_data(1).pDC_Model.CR_OK)
experimento2CR = size(all_data(2).pDC_Model.CR_OK)
experimento3CR = size(all_data(3).pDC_Model.CR_OK)
experimento4CR = size(all_data(1).pDC_Model_unconstrain.CR_OK)
experimento5CR = size(all_data(2).pDC_Model_unconstrain.CR_OK)
experimento6CR = size(all_data(3).pDC_Model_unconstrain.CR_OK)



%% Tabla tiempos

time_matpower1 = round(sum(all_data(1).Matpower_Data.DataMATPOWER.timeAC)/60,4);
time_generationbeta1 = round(sum(all_data(1).Generation_Data.Time)/60,4);
time_nn1 = round(478.40/60,4);

time_matpower2 = round(sum(all_data(2).Matpower_Data.DataMATPOWER.timeAC)/60,4);
time_generationbeta2 = round(sum(all_data(2).Generation_Data.Time)/60,4);
time_nn2 = round(569.50/60,4);

time_matpower3 = round(sum(all_data(3).Matpower_Data.DataMATPOWER.timeAC)/60,4);
time_generationbeta3 = round(sum(all_data(3).Generation_Data.Time)/60,4);
time_nn3 = round(867.233/60,4);

total1 = time_matpower1 + time_generationbeta1 + time_nn1;
total2 = time_matpower2 + time_generationbeta2 + time_nn2;
total3 = time_matpower3 + time_generationbeta3 + time_nn3;

Analysis_Time = table([time_matpower1; time_generationbeta1; time_nn1;total1], ...
                      [time_matpower2; time_generationbeta2; time_nn2;total2], ...
                      [time_matpower3; time_generationbeta3; time_nn3; total3], ...
                      'RowNames', {'$t^{MATPOWER}$', '$t^{GENERATION \beta}$', '$t^{NN}$','$t^{TOTAL}$'}, ...
                      'VariableNames', {'30 IEEE bus', '57 IEEE bus', '118 IEEE bus'});


%% 

% Analisis Revenue Adequency (Con/Sin) Restriccion
Constrained = rand(10,1);
Unconstrain = rand(10,1);


% Calcular mínimo, máximo y media de cada experimento
minimo_exp1 = sprintf('%.4e',min(Constrained));
minimo_exp2 = sprintf('%.4e',min(Unconstrain));


maximo_exp1 = sprintf('%.4e',max(Constrained));
maximo_exp2 = sprintf('%.4e',max(Unconstrain));


media_exp1 = sprintf('%.4e',mean(Constrained));
media_exp2 = sprintf('%.4e',mean(Unconstrain));


%Analisys_Constrained = table([minimo_exp1; minimo_exp2; minimo_exp3], ...
                  % [maximo_exp1; maximo_exp2; maximo_exp3], ...
                  % [media_exp1; media_exp2; media_exp3], ...
                  % 'RowNames', {'Unconstrained', 'Constrained'}, ...
                  % 'VariableNames', {'Minimum', 'Maximum', 'Mean'});



%%  POLTS ZONE
data1 = pglib_opf_case30_ieee;
data2 = pglib_opf_case57_ieee;
data3 = pglib_opf_case118_ieee;
c1 = data1.gencost(:,6)*100;
c2 = data2.gencost(:,6)*100;
c3 = data3.gencost(:,6)*100;

BaseCostAC1 = (all_data(1).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:6)*c1)/1000; 
BaseCostAC2 = (all_data(2).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:7)*c2)/1000;  
BaseCostAC3 = (all_data(3).Matpower_Data.Data1.DataMATPOWER.DispatchAC_Probe(:,1:54)*c3)/1000;

BaseCostpDC1 = (all_data(1).pDC_Model.Pg_Feasibility_Cost(:,1:6)*c1)/1000; 
BaseCostpDC2 = (all_data(2).pDC_Model.Pg_Feasibility_Cost(:,1:7)*c2)/1000; 
BaseCostpDC3 = (all_data(3).pDC_Model.Pg_Feasibility_Cost(:,1:54)*c3)/1000; 


BaseCostDC1 = (all_data(1).DC_Model.Pg_Feasibility_Cost(:,1:6)*c1)/1000; 
BaseCostDC2 = (all_data(2).DC_Model.Pg_Feasibility_Cost(:,1:7)*c2)/1000; 
BaseCostDC3 = (all_data(3).DC_Model.Pg_Feasibility_Cost(:,1:54)*c3)/1000; 


ErrorpDC =abs(BaseCostpDC1 - BaseCostAC1)./BaseCostAC1 *100;
ErrorDC = abs(BaseCostDC1 - BaseCostAC1)./BaseCostAC1 *100;

Contador = 0;
for i = 1:1000
    if ErrorpDC(i)<=ErrorDC(i)
        Contador = Contador +1;
    end
end


%% 


% Plot Costs
figure(1)

% Subplot 1
subplot(3,1,1)
scatter(BaseCostAC1, all_data(1).pDC_Model.BaseCost/1000,'o', 'SizeData', 10);
%scatter(BaseCostAC1, BaseCostpDC1,'o', 'SizeData', 10);
hold on
scatter(BaseCostAC1, all_data(1).DC_Model.BaseCost/1000,'o', 'SizeData', 10);
%scatter(BaseCostAC1, BaseCostDC1,'o', 'SizeData', 10);
hline1 = refline(1,0);
hline1.LineWidth = 1;
hline1.Color = 'k';
% Obtener los límites actuales de los ejes
xlim_current = xlim;
ylim_current = ylim;

% Calcular los nuevos límites de los ejes para que cubran toda la línea de referencia
xlim_new = [min(xlim_current(1), ylim_current(1)), max(xlim_current(2), ylim_current(2))];
ylim_new = xlim_new;

% Establecer los nuevos límites de los ejes
xlim(xlim_new);
ylim(ylim_new);
title('30-bus system')
ylabel('$\textrm{(p)DC-OPF model}$');
xlabel('$\textrm{AC-OPF model}$');
grid on;

% Subplot 2
subplot(3,1,2)
scatter(BaseCostAC2, all_data(2).pDC_Model.BaseCost/1000,'o', 'SizeData', 10);
%scatter(BaseCostAC2, BaseCostpDC2,'o', 'SizeData', 10);
hold on
scatter(BaseCostAC2, all_data(2).DC_Model.BaseCost/1000,'o', 'SizeData', 10);
%scatter(BaseCostAC2, BaseCostDC2,'o', 'SizeData', 10);
hline2 = refline(1,0);
hline2.LineWidth = 1;
hline2.Color = 'k';
% Obtener los límites actuales de los ejes
xlim_current = xlim;
ylim_current = ylim;

% Calcular los nuevos límites de los ejes para que cubran toda la línea de referencia
xlim_new = [min(xlim_current(1), ylim_current(1)), max(xlim_current(2), ylim_current(2))];
ylim_new = xlim_new;

% Establecer los nuevos límites de los ejes
xlim(xlim_new);
ylim(ylim_new);
ylabel('$\textrm{(p)DC-OPF model}$');
xlabel('$\textrm{AC-OPF model}$');
title('57-bus system')
sgtitle('Total cost (thousands of \$)', 'Interpreter', 'latex','FontSize',14);
grid on;

% Subplot 3
subplot(3,1,3)
scatter(BaseCostAC3, all_data(3).pDC_Model.BaseCost/1000,'o', 'SizeData', 10);
%scatter(BaseCostAC3, BaseCostpDC3,'o', 'SizeData', 10);
hold on
scatter(BaseCostAC3, all_data(3).DC_Model.BaseCost/1000, 'o', 'SizeData', 10);
%scatter(BaseCostAC3, BaseCostDC3,'o', 'SizeData', 10);
hline3 = refline(1,0);
hline3.LineWidth = 1;
hline3.Color = 'k';
% Obtener los límites actuales de los ejes
xlim_current = xlim;
ylim_current = ylim;

% Calcular los nuevos límites de los ejes para que cubran toda la línea de referencia
xlim_new = [min(xlim_current(1), ylim_current(1)), max(xlim_current(2), ylim_current(2))];
ylim_new = xlim_new;

% Establecer los nuevos límites de los ejes
xlim(xlim_new);
ylim(ylim_new);
ylabel('$\textrm{(p)DC-OPF model}$');
xlabel('$\textrm{AC-OPF model}$');
title('118-bus system')
grid on;

% Subplot 6
exportgraphics(gcf, 'cost_comparison.pdf', 'ContentType', 'vector');
%Cambiar posteriormente el refline

% 5.032111721056281,9
% 28.970914918998279,45
% 80,110


% figure(2)
% subplot(1,3,1)
% plot(BaseCostAC1, all_data(1).DC_Model.BaseCost,'o','Color', [0.5, 0.25, 0]);
% hline4 = refline(1,0);
% hline4.Color = 'b';
% hline4.LineWidth = 1;
% xlabel(' [\$/hr]');
% ylabel('$Cost^{p-DCOPF}$ [\$/hr]');
% title('30-bus')
% grid on;
% 
% 
% % Subplot 2
% subplot(1,3,2)
% plot(BaseCostAC2, all_data(2).DC_Model.BaseCost,'o','Color', [0.5, 0.25, 0]);
% hline5 = refline(1,0);
% hline5.Color = 'b';
% hline5.LineWidth = 1.1;
% xlabel('[\$/hr]');
% title('57-bus')
% grid on;
% 
% 
% % Subplot 3
% subplot(1,3,3)
% plot(BaseCostAC3, all_data(3).DC_Model.BaseCost,'o','Color', [0.5, 0.25, 0]);
% hline6 = refline(1,0);
% hline6.Color = 'b';
% hline6.LineWidth = 1.1;
% 
% xlabel('[\$/hr]');
% title('118-bus')
% grid on;
% hold off;
% exportgraphics(gcf, 'cost_comparisonDC.pdf', 'ContentType', 'vector');


%% COST RECOVERY


experimento1CR = all_data(1).pDC_Model.CR_Result;
experimento2CR = all_data(2).pDC_Model.CR_Result;
experimento3CR = all_data(3).pDC_Model.CR_Result;
experimento4CR = all_data(1).pDC_Model_unconstrain.CR_Result;
experimento5CR = all_data(2).pDC_Model_unconstrain.CR_Result;
experimento6CR = all_data(3).pDC_Model_unconstrain.CR_Result;
% Calcular mínimo, máximo y media de cada experimento

experimento1CR = reshape(experimento1CR, [], 1);
experimento2CR = reshape(experimento2CR, [], 1);
experimento3CR = reshape(experimento3CR, [], 1);
experimento4CR = reshape(experimento4CR, [], 1);
experimento5CR = reshape(experimento5CR, [], 1);
experimento6CR = reshape(experimento6CR, [], 1);

figure(4)

% Primer subplot
subplot(3,1,1)
histogram(experimento1CR/1000);
grid on
hold on
histogram(experimento4CR/1000, 'FaceAlpha', 0.5); % Segundo histograma más transparente
legend('Constrained','Unconstrained')
title('30-$bus$')

% Segundo subplot
subplot(3,1,2)
histogram(experimento2CR/1000);
hold on
histogram(experimento5CR/1000, 'FaceAlpha', 0.5); % Segundo histograma más transparente
ylabel('Frequency');
title('57-$bus$')
legend('Constrained','Unconstrained')
grid on

% Tercer subplot
subplot(3,1,3)
histogram(experimento3CR/1000);
hold on
histogram(experimento6CR/1000, 'FaceAlpha', 0.5); % Segundo histograma más transparente
xlabel('Thousand of ($)');
legend('Constrained','Unconstrained')
title('118-$bus$')
grid on

exportgraphics(gcf, 'CostRecovery.pdf', 'ContentType', 'vector');

%% 


% Etiquetas para los sistemas y las variables
% Plot Errors
figure(3)

boxchart([experimento1Pg,experimento4Pg]);
grid on
set(gca,'TickLabelInterpreter','latex');
set(gca,'XTickLabel',{'p-DCOPF','DCOPF'});

ylabel('Value [p.u.]');

figure(4)
boxchart([experimento2Pg,experimento5Pg]);
grid on
set(gca,'TickLabelInterpreter','latex');
set(gca,'XTickLabel',{'p-DCOPF','DCOPF'});

ylabel('Value [p.u.]');
grid on
figure(5)
boxchart([experimento3Pg,experimento6Pg]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'XTickLabel',{'p-DCOPF','DCOPF'});

grid on
ylabel('Value [p.u.]');



%h = findobj(gca,'Tag','Box');  % Busca las cajas en el gráfico
%legend(h, ); % Leyenda para las cajas

%boxplot([error1(:,2),error2(:,2)], 'Labels', {'p-DCOPF','DCOPF'},'Widths', 0.5);

%boxplot([error1(:,3),error2(:,3)], 'Labels', {'p-DCOPF','DCOPF'},'Widths', 0.5);


%% 
% Generar datos de ejemplo
data1 = randn(100, 1); % Primera columna
data2 = randn(100, 1) + 2; % Segunda columna

% Crear el boxplot
boxplot([data1, data2], 'Labels', {'Columna 1', 'Columna 2'});

% Añadir leyendas
h = findobj(gca,'Tag','Box');  % Busca las cajas en el gráfico
legend(h, {'Columna 1', 'Columna 2'}); % Leyenda para las cajas
hold on;
%scatter(ones(size(data1)), data1, 'r', 'filled'); % Puntos para la primera columna
%scatter(2*ones(size(data2)), data2, 'b', 'filled'); % Puntos para la segunda columna
%hold off;
legend('show'); % Mostrar la segunda leyenda

%% 

% Plot Correlation with each dispatch

experimento1Pg = all_data(1).pDC_Model.Pg-all_data(1).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:6);
experimento2Pg = all_data(2).pDC_Model.Pg-all_data(2).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:7);
experimento3Pg = all_data(3).pDC_Model.Pg-all_data(3).Matpower_Data.Data1.DataMATPOWER.DispatchAC_Probe(:,1:54);
experimento4Pg = all_data(1).DC_Model.Pg-all_data(1).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:6);
experimento5Pg = all_data(2).DC_Model.Pg-all_data(2).Matpower_Data.DataMATPOWER.DispatchAC_Probe(:,1:7);
experimento6Pg = all_data(3).DC_Model.Pg-all_data(3).Matpower_Data.Data1.DataMATPOWER.DispatchAC_Probe(:,1:54);

% Encontrar las columnas con todas las filas iguales a cero
columnas_a_eliminar1 = all(experimento1Pg == 0, 1);
columnas_a_eliminar2 = all(experimento2Pg == 0, 1);
columnas_a_eliminar3 = all(experimento3Pg == 0, 1);
% Eliminar las columnas identificadas
experimento1Pg(:, columnas_a_eliminar1) = [];
experimento2Pg(:, columnas_a_eliminar2) = [];
experimento3Pg(:, columnas_a_eliminar3) = [];
experimento4Pg(:, columnas_a_eliminar1) = [];
experimento5Pg(:, columnas_a_eliminar2) = [];
experimento6Pg(:, columnas_a_eliminar3) = [];

experimento1Pg = reshape(experimento1Pg, [], 1);
experimento2Pg = reshape(experimento2Pg, [], 1);
experimento3Pg = reshape(experimento3Pg, [], 1);
experimento4Pg = reshape(experimento4Pg, [], 1);
experimento5Pg = reshape(experimento5Pg, [], 1);
experimento6Pg = reshape(experimento6Pg, [], 1);

figure(4)

% Subplot 1
subplot(3,1,1)
histogram(experimento1Pg, 'BinWidth', 0.01, 'Normalization', 'probability')
title('30-bus system', 'FontSize', 16)
ylabel('Density', 'FontSize', 17)
xlabel('$$p_k^{\rm g} - p_k^{\rm g,*} \,\textrm{(p.u.)}$$', 'FontSize', 17)
grid on
hold on
histogram(experimento4Pg, 'BinWidth', 0.01, 'Normalization', 'probability')
legend({'pDC-OPF', 'DC-OPF'}, 'Location', 'northoutside', 'Orientation', 'horizontal')
%set(gca, 'YScale', 'log') % Establecer escala logarítmica en el eje y
%ylim([10^-5, 10^0]); % Ajuste de los límites del eje y

% Subplot 2
subplot(3,1,2)
histogram(experimento2Pg, 'BinWidth', 0.5, 'Normalization', 'probability')
title('57-bus system', 'FontSize', 16)
ylabel('Density', 'FontSize', 17)
xlabel('$$p_k^{\rm g} - p_k^{\rm g,*} \,\textrm{(p.u.)}$$', 'FontSize', 17)
grid on
hold on
histogram(experimento5Pg, 'BinWidth', 0.5, 'Normalization', 'probability')
%set(gca, 'YScale', 'log') % Establecer escala logarítmica en el eje y
%ylim([10^-5, 10^0]); % Ajuste de los límites del eje y

% Subplot 3
subplot(3,1,3)
histogram(experimento3Pg, 'BinWidth',0.5, 'Normalization', 'probability')
title('118-bus system', 'FontSize', 16)
ylabel('Density', 'FontSize', 17)
xlabel('$$p_k^{\rm g} - p_k^{\rm g,*} \,\textrm{(p.u.)}$$', 'FontSize', 17)
grid on
hold on
histogram(experimento6Pg, 'BinWidth', 0.5, 'Normalization', 'probability')
%set(gca, 'YScale', 'log') % Establecer escala logarítmica en el eje y
%ylim([10^-5, 10^0]); % Ajuste de los límites del eje y

exportgraphics(gcf, 'RelativePG.pdf', 'ContentType', 'vector');


%CAMBIAR A 18
%% ELECTRICITY MARKETS


experimento1RA = all_data(1).pDC_Model.RA_Result;
experimento2RA = all_data(2).pDC_Model.RA_Result;
experimento3RA = all_data(3).pDC_Model.RA_Result;
experimento4RA = all_data(1).pDC_Model_unconstrain.RA_Result;
experimento5RA = all_data(2).pDC_Model_unconstrain.RA_Result;
experimento6RA = all_data(3).pDC_Model_unconstrain.RA_Result;


% Revenue Adequency
figure(6)
subplot(3,2,1)
boxplot(all_data(1).pDC_Model.RA_Result/1000);

title('30-$bus$')
grid on
subplot(3,2,2)
boxplot(all_data(1).pDC_Model_unconstrain.RA_Result/1000);
title('30-$bus$')
grid on
subplot(3,2,3)
boxplot(all_data(2).pDC_Model.RA_Result/1000);
ylabel('Thousand of ($)');
title('57-$bus$')
grid on
subplot(3,2,4)
boxplot(all_data(2).pDC_Model_unconstrain.RA_Result/1000);
title('57-$bus$')
grid on

subplot(3,2,5)
boxplot(all_data(3).pDC_Model.RA_Result/1000);
xlabel('$pDC-OPF_{\mathcal{C}}$')
title('118-$bus$')
grid on

subplot(3,2,6)
boxplot(all_data(3).pDC_Model_unconstrain.RA_Result/1000);
xlabel('$pDC-OPF_{\mathcal{U}}$')
title('118-$bus$')
grid on
exportgraphics(gcf, 'RevenueAdequency.pdf', 'ContentType', 'vector');
%% 

figure(7)
% Grafica de Factibilidad para modelos p-DC y DC
percentiles = 0:1:100; % Porcentajes del 0 al 100
Porcentiles_pDC1 = prctile(sqrt(all_data(1).pDC_Model.Feasibility/6), percentiles); % Calcular los valores correspondientes a los percentiles
Porcentiles_DC1 = prctile(sqrt(all_data(1).DC_Model.Feasibility/6), percentiles); 
Porcentiles_pDC2 = prctile(sqrt(all_data(2).pDC_Model.Feasibility/7), percentiles); % Calcular los valores correspondientes a los percentiles
Porcentiles_DC2 = prctile(sqrt(all_data(2).DC_Model.Feasibility/7), percentiles); 
Porcentiles_pDC3 = prctile(sqrt(all_data(3).pDC_Model.Feasibility/54), percentiles); % Calcular los valores correspondientes a los percentiles
Porcentiles_DC3 = prctile(sqrt(all_data(3).DC_Model.Feasibility/54), percentiles); 

% Porcentiles_pDC1 = prctile(100*sum(abs(all_data(1).pDC_Model.Pg-all_data(1).pDC_Model.Pg_Feasibility),2)./sum(all_data(1).pDC_Model.Pg_Feasibility,2), percentiles); % Calcular los valores correspondientes a los percentiles
% Porcentiles_DC1 = prctile(100*sum(abs(all_data(1).DC_Model.Pg-all_data(1).DC_Model.Pg_Feasibility),2)./sum(all_data(1).DC_Model.Pg_Feasibility,2), percentiles); 
% Porcentiles_pDC2 = prctile(100*sum(abs(all_data(2).pDC_Model.Pg-all_data(2).pDC_Model.Pg_Feasibility),2)./sum(all_data(2).pDC_Model.Pg_Feasibility,2), percentiles); % Calcular los valores correspondientes a los percentiles
% Porcentiles_DC2 = prctile(100*sum(abs(all_data(2).DC_Model.Pg-all_data(2).DC_Model.Pg_Feasibility),2)./sum(all_data(2).DC_Model.Pg_Feasibility,2), percentiles); 
% Porcentiles_pDC3 = prctile(100*sum(abs(all_data(3).pDC_Model.Pg-all_data(3).pDC_Model.Pg_Feasibility),2)./sum(all_data(3).pDC_Model.Pg_Feasibility,2), percentiles); % Calcular los valores correspondientes a los percentiles
% Porcentiles_DC3 = prctile(100*sum(abs(all_data(3).DC_Model.Pg-all_data(3).DC_Model.Pg_Feasibility),2)./sum(all_data(3).DC_Model.Pg_Feasibility,2), percentiles); 




% Porcentiles_NN = prctile(FeseabilityNN, percentiles); 
h1 = stairs(Porcentiles_pDC1, percentiles, '-', 'LineWidth', 1.5, 'Color', "#0072BD");
hold on
h2 = stairs(Porcentiles_DC1, percentiles, '--', 'LineWidth', 1.5, 'Color', "#0072BD");
h3 = stairs(Porcentiles_pDC2, percentiles, '-', 'LineWidth', 1.5, 'Color', "#D95319");
h4 = stairs(Porcentiles_DC2, percentiles, '--', 'LineWidth', 1.5, 'Color', "#D95319");
h5 = stairs(Porcentiles_pDC3, percentiles, '-', 'LineWidth', 1.5, 'Color', "#77AC30");
h6 = stairs(Porcentiles_DC3, percentiles, '--', 'LineWidth', 1.5, 'Color', "#77AC30");

set(gca, 'XScale', 'log'); % Cambiar la escala del eje X a logarítmica
% stairs(Porcentiles_NN, percentiles, '--');

xlabel('Distance to feasibility (p.u.)', 'FontSize', 19); % Cambiar el tamaño de la fuente del label del eje X
ylabel('Percentage of instances', 'FontSize', 20); % Cambiar el tamaño de la fuente del label del eje Y

legend([h1, h3, h5], {'30-bus system', '57-bus system', '118-bus system'}, 'FontSize', 14); % Ajustar la leyenda y su tamaño de fuente

grid on
exportgraphics(gcf, 'Feasibility.pdf', 'ContentType', 'vector');


%legend('$$30bus^{p-DCOPF}$$','$$30bus^{DCOPF}$$','$$57bus^{p-DCOPF}$$','$$57bus^{DCOPF}$$','$$118bus^{p-DCOPF}$$','$$118bus^{DCOPF}$$')

%% CONTINGENCIAS

% Matriz de errores de costo de contingencias 1
% Preasignación de la matriz de errores
Errores_Contingencia_COSTPDC_1 = zeros(1000, 41); % Ajusta num_filas según sea necesario
Errores_Contingencia_COSTDC_1 = zeros(1000, 1); % Ajusta num_filas según sea necesario
Errores_Contingencia_COSTPDC_2 = zeros(1000, 80); % Ajusta num_filas según sea necesario
Errores_Contingencia_COSTDC_2 = zeros(1000, 1); % Ajusta num_filas según sea necesario
Errores_Contingencia_COSTPDC_3 = zeros(1000, 186); % Ajusta num_filas según sea necesario
Errores_Contingencia_COSTDC_3 = zeros(1000, 1); % Ajusta num_filas según sea necesario

%Flag30 = [0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1];
%Flag57 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0];

%Flag118  = [1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];



for l = 1:41
    % Construir el nombre de la variable
    variable_name = sprintf('all_data(1).ContingenciapDC_Data.ErrorCost_%d.Absolute_Error', l);
    
    % Acceder al valor de la variable usando getfield
    valor_variable = getfield(all_data(1).ContingenciapDC_Data, ['ErrorCost_' num2str(l)], 'Absolute_Error');
    
    % Asignar el valor a la matriz de errores
    Errores_Contingencia_COSTPDC_1(:, l) = valor_variable;
end



% Matriz de errores de costo de contingencias 1
for l = 1:41
    variable_name = sprintf('all_data(1).ContingenciaDC_Data.ErrorCost_%d.Absolute_Error', l);
    Errores_Contingencia_COSTDC_1(:, l) = getfield(all_data(1).ContingenciaDC_Data, ['ErrorCost_' num2str(l)], 'Absolute_Error');
end


% Matriz de errores de costo de contingencias 2 (PDC)
for l = 1:80
    variable_name = sprintf('all_data(2).ContingenciapDC_Data.ErrorCost_%d.Absolute_Error', l);
    Errores_Contingencia_COSTPDC_2(:, l) = getfield(all_data(2).ContingenciapDC_Data, ['ErrorCost_' num2str(l)], 'Absolute_Error');
end

% Matriz de errores de costo de contingencias 2 (DC)
for l = 1:80
    variable_name = sprintf('all_data(2).ContingenciaDC_Data.ErrorCost_%d.Absolute_Error', l);
    Errores_Contingencia_COSTDC_2(:, l) = getfield(all_data(2).ContingenciaDC_Data, ['ErrorCost_' num2str(l)], 'Absolute_Error');
end


% Matriz de errores de costo de contingencias 3 (PDC)
for l = 1:186
    variable_name = sprintf('all_data(3).ContingenciapDC_Data.ErrorCost_%d.Absolute_Error', l);
    Errores_Contingencia_COSTPDC_3(:, l) = getfield(all_data(3).ContingenciapDC_Data, ['ErrorCost_' num2str(l)], 'Absolute_Error');
end


% Matriz de errores de costo de contingencias 3 (DC)
for l = 1:186
    variable_name = sprintf('all_data(3).ContingenciaDC_Data.ErrorCost_%d.Absolute_Error', l);
    Errores_Contingencia_COSTDC_3(:, l) = getfield(all_data(3).ContingenciaDC_Data, ['ErrorCost_' num2str(l)], 'Absolute_Error');
end

%Malas_muestras = [6, 7, 8, 50,94, 98, 124, 125, 178, 179, 184, 185, 186]
%Errores_Contingencia_COSTDC_3(:,Malas_muestras) = [];
%Errores_Contingencia_COSTPDC_3(:,Malas_muestras) = [];

% figure(6)
% subplot(2,1,1)
% boxplot(Errores_Contingencia_COSTPDC_1(:,1:40))
% grid on
% ylabel('$\epsilon^{pDC-OPF}$ [\%]')
% subplot(2,1,2)
% boxplot(Errores_Contingencia_COSTDC_1(:,1:40))
% grid on
% ylabel('$\epsilon^{DC-OPF}$ [\%]')
% exportgraphics(gcf, 'Contingency_30.pdf', 'ContentType', 'vector');
% figure(7)
% subplot(2,1,1)
% boxplot(Errores_Contingencia_COSTPDC_2(:,1:40))
% grid on
% ylabel('$\epsilon^{pDC-OPF}$ [\%]')
% subplot(2,1,2)
% boxplot(Errores_Contingencia_COSTDC_2(:,1:40))
% grid on
% ylabel('$\epsilon^{DC-OPF}$ [\%]')
% exportgraphics(gcf, 'Contingency_57.pdf', 'ContentType', 'vector');
% figure(8)
% subplot(2,1,1)
% boxplot(Errores_Contingencia_COSTPDC_3(:,1:40))
% grid on
% ylabel('$\epsilon^{pDC-OPF}$ [\%]')
% subplot(2,1,2)
% boxplot(Errores_Contingencia_COSTDC_3(:,1:40))
% grid on
% ylabel('$\epsilon^{DC-OPF}$ [\%]')
% exportgraphics(gcf, 'Contingency_118.pdf', 'ContentType', 'vector');
flag30 = readmatrix('flagAC30.csv');
fag57 = readmatrix('flagAC57.csv');
flag118 = readmatrix('flagAC118.csv');

flag30 = reshape(flag30,[],1);
 flag57 = reshape(fag57,[],1);
 flag118 = reshape(flag118,[],1);


mask30 = flag30 == 1;
mask57 = flag57 ==1;
mask118 = flag118 ==1;


x1 = all_data(1).pDC_Model.ErrorCost.Absolute_Error;

%x1(x1>=100)=[];
x2 = reshape(Errores_Contingencia_COSTPDC_1, [], 1);
x2 = x2(mask30);
%x2(x2>=100)=[];
%x3(x3>=100)=[];
x3 =all_data(2).pDC_Model.ErrorCost.Absolute_Error;

x4 = reshape(Errores_Contingencia_COSTPDC_2, [], 1);
x4 = x4(mask57);
indices_a_eliminar = find(x4 == 100);
x4(indices_a_eliminar) = [];

x10 = reshape(Errores_Contingencia_COSTDC_2, [], 1);
x10 = x10(mask57);
x10(indices_a_eliminar) = [];


%x3(indices_a_eliminar) = [];
%x4(x4>=100)=[];
x5 = all_data(3).pDC_Model.ErrorCost.Absolute_Error;
%x5(x5>=100)=[];
x6 = reshape(Errores_Contingencia_COSTPDC_3, [], 1);
x6 = x6(mask118);
indices_a_eliminar = find(x6 == 100);
x6(indices_a_eliminar) = [];
%x6(x6>=100)=[];

x7 = all_data(1).DC_Model.ErrorCost.Absolute_Error;
x8 = reshape(Errores_Contingencia_COSTDC_1, [], 1);
x8 = x8(mask30);

x9 = all_data(2).DC_Model.ErrorCost.Absolute_Error;


x11 = all_data(3).DC_Model.ErrorCost.Absolute_Error;
x12 = reshape(Errores_Contingencia_COSTDC_3, [], 1);
x12 = x12(mask118);
x12(indices_a_eliminar) = [];

figure(9)

% Subplot 1
subplot(3,1,1)
x1 = sort(x1);
n = length(x1);
empirical_cdf = (1:n)' / n;
h1 = plot(x1, empirical_cdf, 'Color', "#0072BD");
hold on

[posicion_fila, posicion_columna]  = find(Errores_Contingencia_COSTDC_1 == max(max(x8)));

x2 = sort(x2);
n = length(x2);
empirical_cdf = (1:n)' / n;
plot(x2, empirical_cdf, '--', 'Color', "#0072BD")

x7 = sort(x7);
n = length(x7);
empirical_cdf = (1:n)' / n;
h2 = plot(x7, empirical_cdf, 'Color', "#D95319");

x8 = sort(x8);
n = length(x8);
empirical_cdf = (1:n)' / n;
plot(x8, empirical_cdf, '--', 'Color', "#D95319")

ylabel('eCDF', 'FontSize', 19)
xlabel('Cost approximation error $(\%)$', 'Interpreter', 'latex', 'FontSize', 19)
title('30-bus system', 'FontSize', 18)
grid on

% Subplot 2
subplot(3,1,2)
x3 = sort(x3);
n = length(x3);
empirical_cdf = (1:n)' / n;
plot(x3, empirical_cdf, 'Color', "#0072BD")
hold on

x4 = sort(x4);
n = length(x4);
empirical_cdf = (1:n)' / n;
plot(x4, empirical_cdf, '--', 'Color', "#0072BD")

x9 = sort(x9);
n = length(x9);
empirical_cdf = (1:n)' / n;
plot(x9, empirical_cdf, 'Color', "#D95319")

x10 = sort(x10);
n = length(x10);
empirical_cdf = (1:n)' / n;
plot(x10, empirical_cdf, '--', 'Color', "#D95319")

ylabel('eCDF', 'FontSize', 19)
xlabel('Cost approximation error $(\%)$', 'Interpreter', 'latex', 'FontSize', 19)
title('57-bus system', 'FontSize', 18)
grid on

% Subplot 3
subplot(3,1,3)
x5 = sort(x5);
n = length(x5);
empirical_cdf = (1:n)' / n;
plot(x5, empirical_cdf, 'Color', "#0072BD")
hold on

x6 = sort(x6);
n = length(x6);
empirical_cdf = (1:n)' / n;
plot(x6, empirical_cdf, '--', 'Color', "#0072BD")

x11 = sort(x11);
n = length(x11);
empirical_cdf = (1:n)' / n;
plot(x11, empirical_cdf, 'Color', "#D95319")

x12 = sort(x12);
n = length(x12);
empirical_cdf = (1:n)' / n;
plot(x12, empirical_cdf, '--', 'Color', "#D95319")

ylabel('eCDF', 'FontSize', 19)
xlabel('Cost approximation error $(\%)$', 'Interpreter', 'latex', 'FontSize', 19)
title('118-bus system', 'FontSize', 18)
grid on

legend([h1, h2], 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 15);


%mu = mean(x5);
% x5= sort(x5);
% n = length(x5);
% empirical_cdf = (1:n)' / n;
% %plot(x5,empirical_cdf)
% boxplot(x5)
% title('118 bus')
% ylabel('Percentage of cost error')
% xlabel('Base Model')
% grid on
% subplot(1,2,2)
% mu = mean(x6);
% x6= sort(x6);
% n = length(x6);
% empirical_cdf = (1:n)' / n;
% %plot(x6,empirical_cdf)
% boxplot(x6)
% 
% 
% title('118 bus')
% xlabel('Topology N-1')
% grid on
exportgraphics(gcf, 'Contingency.pdf', 'ContentType', 'vector');
%% 

percentiles = 98:0.1:100; % Porcentajes del 0 al 100

% Calcular los valores correspondientes a los percentiles
Porcentiles_pDCO30 = prctile(x1, percentiles); 
Porcentiles_pDC30 = prctile(x2, percentiles); 
Porcentiles_DCO30 = prctile(x7, percentiles); 
Porcentiles_DC30 = prctile(x8, percentiles); 
Porcentiles_pDCO57 = prctile(x3, percentiles); 
Porcentiles_pDC57 = prctile(x4, percentiles); 
Porcentiles_DCO57 = prctile(x9, percentiles); 
Porcentiles_DC57 = prctile(x10, percentiles); 
Porcentiles_pDCO118 = prctile(x5, percentiles); 
Porcentiles_pDC118 = prctile(x6, percentiles); 
Porcentiles_DCO118 = prctile(x11, percentiles); 
Porcentiles_DC118 = prctile(x12, percentiles); 

% Truncar los valores a 3 decimales
Porcentiles_pDCO30 = round(Porcentiles_pDCO30, 3);
Porcentiles_pDC30 = round(Porcentiles_pDC30, 3);
Porcentiles_DCO30 = round(Porcentiles_DCO30, 3);
Porcentiles_DC30 = round(Porcentiles_DC30, 3);
Porcentiles_pDCO57 = round(Porcentiles_pDCO57, 3);
Porcentiles_pDC57 = round(Porcentiles_pDC57, 3);
Porcentiles_DCO57 = round(Porcentiles_DCO57, 3);
Porcentiles_DC57 = round(Porcentiles_DC57, 3);
Porcentiles_pDCO118 = round(Porcentiles_pDCO118, 3);
Porcentiles_pDC118 = round(Porcentiles_pDC118, 3);
Porcentiles_DCO118 = round(Porcentiles_DCO118, 3);
Porcentiles_DC118 = round(Porcentiles_DC118, 3);

% Crear las matrices
Matriz_30 = [Porcentiles_pDCO30', Porcentiles_DCO30', Porcentiles_pDC30', Porcentiles_DC30'];
Matriz_57 = [Porcentiles_pDCO57', Porcentiles_DCO57', Porcentiles_pDC57', Porcentiles_DC57'];
Matriz_118 = [Porcentiles_pDCO118', Porcentiles_DCO118', Porcentiles_pDC118', Porcentiles_DC118'];


%% 

figure(9)

% Matriz comparacion de Despachos de contingencias
for l = 1:56
    i=nbranch_depure(l);
    variable_name1 = sprintf('Contingencias_PDC.DISPATCH_Line%d',i);
    variable_name2 = sprintf('Contingencias_DC.DISPATCH_Line%d',i);
    variable_name3 = sprintf('Contingencias_Real.DISPATCH_AC_Line%d',i);

    subplot(7,8, l);
    histogram(sum(eval(variable_name3),2));
    xlabel('Valores');

    hold on
     histogram(sum(eval(variable_name2),2));
    grid on
    set(gca, 'YTickLabel', []); % Removes y-axis tick labels
    set(gca, 'XLabel', []); % Set XLabel to empty []
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]); % Maximize figure window
end


figure(13)

scatter(Duales_Matpower.Duales_Matpower.DualAC,pDC_Model.Dual)
grid on
ylabel("Duales AC ($/hr)")
xlabel("Duales pDC ($/hr)")
refline(1,0)

figure(14)

scatter(Duales_Matpower.Duales_Matpower.DualAC,DC_Model.Dual)
grid on
ylabel("Duales AC ($/hr)")
xlabel("Duales DC ($/hr)")
refline(1,0)

figure(15)

scatter(Duales_Matpower.Duales_Matpower.DualDC,pDC_Model.Dual)
grid on
refline(1,0)
ylabel("Duales DC MATPOWER ($/hr)")
xlabel("Duales pDC ($/hr)")

% Matriz de errores de costo de contingencias
for l = 1:56
    i=nbranch_depure(l);
    variable_name1 = sprintf('Contingencias_PDC.Objective_Feseability%d',i);
    Feseability_Cont_PDC(:,l) =eval(variable_name1);
end
Feseability_Cont_PDC= [FeseabilityPDC,Feseability_Cont_PDC(:,1:3)];
figure(10)
boxplot(Feseability_Cont_PDC)
ylabel("Distance to feasibility")
grid on

%% 
LMPSAC1 = all_data(1).Matpower_Data.DataMATPOWER.DualAC_Probe;          
LMPSpDC1 = all_data(1).pDC_Model.Dual / 100;
LMPSDC1 = all_data(1).Matpower_Data.DataMATPOWER.DualDC_Probe;

LMPSAC2 = all_data(2).Matpower_Data.DataMATPOWER.DualAC_Probe;      
LMPSpDC2 = all_data(2).pDC_Model.Dual / 100;
LMPSDC2 = all_data(2).Matpower_Data.DataMATPOWER.DualDC_Probe;

LMPSAC3 = all_data(3).Matpower_Data.Data1.DataMATPOWER.DualAC_Probe;     
LMPSpDC3 = all_data(3).pDC_Model.Dual / 100;
LMPSDC3 = all_data(3).Matpower_Data.Data1.DataMATPOWER.DualDC_Probe;

figure(16)

% Definir colores
color_pDC_OPF = [0 0.4470 0.7410];  % Azul claro
color_AC_OPF = [0.8500 0.3250 0.0980];  % Rojo anaranjado
color_DC_OPF = [0.9290 0.6940 0.1250];  % Verde

% Subplot 1
figure(4)

% Subplot 1
subplot(3,1,1)
[counts1, edges1] = histcounts(LMPSDC1, 'BinWidth', 1, 'Normalization', 'probability');
edges1 = edges1(1:end-1) + diff(edges1)/2; % Ajustar los bordes al centro de los bins
stairs(edges1, counts1, 'Color', color_DC_OPF, 'LineWidth', 1.5); % Graficar densidad

hold on
[counts2, edges2] = histcounts(LMPSAC1, 'BinWidth', 1, 'Normalization', 'probability');
edges2 = edges2(1:end-1) + diff(edges2)/2;
stairs(edges2, counts2, 'Color', color_AC_OPF, 'LineWidth', 1.5);

[counts3, edges3] = histcounts(LMPSpDC1, 'BinWidth', 1, 'Normalization', 'probability');
edges3 = edges3(1:end-1) + diff(edges3)/2;
stairs(edges3, counts3, 'Color', color_pDC_OPF, 'LineWidth', 1.5);

ylabel('Density', 'FontSize', 20)
xlabel('Locational marginal prices (\$/MW)', 'Interpreter', 'latex', 'FontSize', 20)
legend({'DC-OPF', 'AC-OPF', 'pDC-OPF'}, 'Location', 'northoutside', 'Orientation', 'horizontal')
title('30-bus system', 'FontSize', 19)
%set(gca, 'YScale', 'log')
grid on

% Subplot 2
subplot(3,1,2)
[counts1, edges1] = histcounts(LMPSDC2, 'BinWidth', 1, 'Normalization', 'probability');
edges1 = edges1(1:end-1) + diff(edges1)/2;
stairs(edges1, counts1, 'Color', color_DC_OPF, 'LineWidth', 1.5);

hold on
[counts2, edges2] = histcounts(LMPSAC2, 'BinWidth', 1, 'Normalization', 'probability');
edges2 = edges2(1:end-1) + diff(edges2)/2;
stairs(edges2, counts2, 'Color', color_AC_OPF, 'LineWidth', 1.5);

[counts3, edges3] = histcounts(LMPSpDC2, 'BinWidth', 1, 'Normalization', 'probability');
edges3 = edges3(1:end-1) + diff(edges3)/2;
stairs(edges3, counts3, 'Color', color_pDC_OPF, 'LineWidth', 1.5);

ylabel('Density', 'FontSize', 20)
xlabel('Locational marginal prices (\$/MW)', 'Interpreter', 'latex', 'FontSize', 20)
title('57-bus system', 'FontSize', 19)
%set(gca, 'YScale', 'log')
grid on

% Subplot 3
subplot(3,1,3)
[counts1, edges1] = histcounts(LMPSDC3, 'BinWidth', 1, 'Normalization', 'probability');
edges1 = edges1(1:end-1) + diff(edges1)/2;
stairs(edges1, counts1, 'Color', color_DC_OPF, 'LineWidth', 1.5);

hold on
[counts2, edges2] = histcounts(LMPSAC3, 'BinWidth', 1, 'Normalization', 'probability');
edges2 = edges2(1:end-1) + diff(edges2)/2;
stairs(edges2, counts2, 'Color', color_AC_OPF, 'LineWidth', 1.5);

[counts3, edges3] = histcounts(LMPSpDC3, 'BinWidth', 1, 'Normalization', 'probability');
edges3 = edges3(1:end-1) + diff(edges3)/2;
stairs(edges3, counts3, 'Color', color_pDC_OPF, 'LineWidth', 1.5);

ylabel('Density', 'FontSize', 20)
xlabel('Locational marginal prices (\$/MW)', 'Interpreter', 'latex', 'FontSize', 20)
title('118-bus system', 'FontSize', 19)
%set(gca, 'YScale', 'log')
grid on



% subplot(3,2,3)
% histogram(LMPSAC2,'BinWidth', 4)
% hold on 
% histogram(LMPSpDC2,'BinWidth', 4)
% legend('AC-OPF','pDC-OPF')
% ylabel('Frequency')
% title('57 bus')
% grid on 
% 
% 
% subplot(3,2,4)
% histogram(LMPSDC2,'BinWidth', 4)
% hold on 
% histogram(LMPSpDC2,'BinWidth', 4)
% legend('DC-OPF','pDC-OPF')
% title('57 bus')
% grid on 
% 
% 
% subplot(3,2,5)
% histogram(LMPSAC3,'BinWidth', 4)
% hold on 
% histogram(LMPSpDC3,'BinWidth', 4)
% legend('AC-OPF','pDC-OPF')
% title('118 bus')
% xlabel('Marginal Prices (\$/p.u.)')
% grid on 
% 
% 
% 
% subplot(3,2,6)
% histogram(LMPSDC3,'BinWidth', 4)
% hold on 
% histogram(LMPSpDC3,'BinWidth', 4)
% legend('DC-OPF','pDC-OPF')
% title('118 bus')
% xlabel('Marginal Prices (\$/p.u.)')
% grid on

% % Obtener los límites actuales de los ejes
% xlim_current = xlim;
% ylim_current = ylim;
% 
% % Calcular los nuevos límites de los ejes para que cubran toda la línea de referencia
% xlim_new = [min(xlim_current(1), ylim_current(1)), max(xlim_current(2), ylim_current(2))];
% ylim_new = xlim_new;
% 
% % Establecer los nuevos límites de los ejes
% xlim(xlim_new);
% ylim(ylim_new);
exportgraphics(gcf,'MarginalCosts.pdf', 'ContentType', 'vector');