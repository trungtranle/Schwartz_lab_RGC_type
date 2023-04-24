%%
load April21_63cells_68datasets.mat
%%
tic;
R_array = zeros(68,1);
fit_array = zeros(68,1);
tau_array = zeros(68,1);
tau_sd_array = zeros(68,1);
capacitance_array = zeros(68,1);

for i = 1:68;
    R = FeatureExtract(results(i));
    R_array(i) = R.resistance.R_MOhm;
    fit_array(i) = R.resistance.AdjustedRsqaured;
 
    tau_array(i) = R.tau.Tau_ms;
    tau_sd_array(i) = R.tau.SD_ms;
    capacitance_array(i) = R.capacitance_pF;
end

toc
%%
tbl = table([1:68]',{results.cell_type}',R_array, fit_array,tau_array, tau_sd_array, capacitance_array, ...
    'VariableNames', {'No','Cell_type','Resistance', 'Rsquared', 'Tau', 'Tau_sd', 'Capacitance'})

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