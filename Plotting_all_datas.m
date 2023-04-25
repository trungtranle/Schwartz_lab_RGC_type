%%
load April21_63cells_68datasets.mat
%%
tic;
R_array = zeros(68,1);
fit_array = zeros(68,1);
tau_array = zeros(68,1);
tau_sd_array = zeros(68,1);
capacitance_array = zeros(68,1);
max_depol_array = zeros(68,1);
min_depol_array = zeros(68,1);
sag_coeff_array = zeros(68,1);

for i = 1:68;
    R = FeatureExtract(results(i));
    R_array(i) = R.resistance.R_MOhm;
    fit_array(i) = R.resistance.AdjustedRsqaured;
 
    tau_array(i) = R.tau.Tau_ms;
    tau_sd_array(i) = R.tau.SD_ms;
    capacitance_array(i) = R.capacitance_pF;
    max_depol_array(i) = max(R.depol_current_level);
    min_depol_array(i) = min(R.depol_current_level);
    sag_coeff_array(i) = table2array(R.sag_coeff);
end

toc
%%
tbl = table([1:68]',{results.cell_type}',R_array, fit_array,tau_array, tau_sd_array, ...
    capacitance_array, max_depol_array, min_depol_array,sag_coeff_array,...
    'VariableNames', {'No','Cell_type','Resistance', 'Rsquared', ...
    'Tau', 'Tau_sd', 'Capacitance', 'Max_Depol_pA', 'Min_Depol_pA', 'Sag_coeff'})

%%
figure;
boxplot(tbl.Resistance(tbl.Rsquared > 0.9), tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Input Resistance (M\Omega)');
box off;

%%
figure;
boxplot(tbl.Tau, tbl.Cell_type);
xlabel('Cell types');
ylabel('Tau (ms)');
box off;
%%
figure;
boxplot(tbl.Capacitance(tbl.Rsquared > 0.9), tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Capacitance (pF)');
box off;
%%

figure;
tiledlayout('vertical');
nexttile;
boxplot(tbl.Max_Depol_pA(tbl.Rsquared > 0.9), tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Max depol current injected (pA)');
title('Raw values');
set(gca,'FontSize', 15, 'XTickLabel', {});


nexttile;
boxplot(tbl.Min_Depol_pA(tbl.Rsquared > 0.9), tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Min depol current injected (pA)');

set(gca,'FontSize', 15);

%%
figure;
tiledlayout('vertical');
nexttile;
boxplot(tbl.Max_Depol_pA(tbl.Rsquared > 0.9) .* tbl.Resistance(tbl.Rsquared > 0.9) , tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Max current injected * Resistance');
title('Normalized by multiple by Resistance');
set(gca,'FontSize', 15, 'XTickLabel', {});

box off;
nexttile;
boxplot(tbl.Min_Depol_pA(tbl.Rsquared > 0.9) .* tbl.Resistance(tbl.Rsquared > 0.9) , tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Min current injected * Resistance');
set(gca,'FontSize', 15);
box off;
%%
figure;
boxplot(tbl.Sag_coeff(tbl.Rsquared > 0.9) , tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Sag coefficient');
set(gca,'FontSize', 15);
