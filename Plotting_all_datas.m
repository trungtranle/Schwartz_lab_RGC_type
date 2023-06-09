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
spontaneous_spike = zeros(68,1);
resting_Vm = zeros(68,1);
V_threshold_array = zeros(68,1);
half_width_time_array = zeros(68,1);
first_AP_peak_array = zeros(68,1);
first_AP_trough_array = zeros(68,1);
spike_number = {};
latency_to_first_spike = {};
adaptaion_index = {};
ISI_cv = {};
block_current_level = zeros(68,1);
max_slope_array = zeros(68,1);
max_spikes_array = zeros(68,1);

for cell = [1:57, 59:68]; %58 is weird
    
    R = FeatureExtract(results(cell));
    R_array(cell) = R.resistance.R_MOhm;
    fit_array(cell) = R.resistance.AdjustedRsqaured;
    tau_array(cell) = R.tau.Tau_ms;
    tau_sd_array(cell) = R.tau.SD_ms;
    capacitance_array(cell) = R.capacitance_pF;
    max_depol_array(cell) = max(R.depol_current_level);
    min_depol_array(cell) = min(R.depol_current_level);
    sag_coeff_array(cell) = table2array(R.sag_coeff);
    resting_Vm(cell) = R.resting_Vm;
    spontaneous_spike(cell) = R.spontaneous_spike;
    V_threshold_array(cell) = R.V_threshold;
    half_width_time_array(cell) = R.half_width_time;
    first_AP_peak_array(cell) = R.first_AP_peak(1);
    first_AP_trough_array(cell) = R.first_AP_trough(1);
    spike_number(cell) = {R.spike_number};
    latency_to_first_spike(cell) = {R.latency_to_first_spike};
    adaptaion_index(cell) = {R.adaptaion_index};
    ISI_cv(cell) = {R.ISI_cv};
    block_current_level(cell) = R.block_current_level;
    max_slope_array(cell) = R.max_slope;
    max_spikes_array(cell) = max(R.spike_number);
end

toc

%%
tbl = table({results.cell_type}',R_array, fit_array,tau_array, tau_sd_array, ...
    capacitance_array, max_depol_array, min_depol_array,sag_coeff_array, spontaneous_spike, resting_Vm, ...
    V_threshold_array, half_width_time_array, first_AP_peak_array, first_AP_trough_array, ...
    spike_number', latency_to_first_spike', adaptaion_index', ISI_cv', block_current_level, max_slope_array, max_spikes_array,...
    'VariableNames', {'Cell_type','Resistance', 'Rsquared', ...
    'Tau', 'Tau_sd', 'Capacitance', 'Max_Depol_pA', 'Min_Depol_pA', 'Sag_coeff', 'Spontaneous_spike', 'Resting_Vm', ...
    'V_threshold', 'Half_width_AP', 'First_AP_Peak', 'First_AP_trough', ...
    'Spike_numbers', 'Latency_to_first_spike', 'Adaptation_index', 'ISI_cv', 'Blocked_current_level', 'Max_slope', 'Max_spike_number'})

%%
figure;
boxplot(tbl.Resistance(tbl.Rsquared > 0.9), tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Input Resistance (M\Omega)');tbl
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
%%
figure;
boxplot(tbl.Spontaneous_spike(tbl.Rsquared > 0.9) , tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Spontaneous spike (Hz)');
set(gca,'FontSize', 15);
%%
figure;
boxplot(tbl.Resting_Vm(tbl.Rsquared > 0.9) , tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Resting V_m (mV)');
set(gca,'FontSize', 15);
%%
max(tbl.Min_Depol_pA)
min(tbl.Max_Depol_pA)
%%
figure;
boxplot(tbl.V_threshold(tbl.Rsquared > 0.9) , tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('AP threshold (mV)');
set(gca,'FontSize', 15);
%%
figure;
boxplot(tbl.Half_width_AP(tbl.Rsquared > 0.9) , tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('AP halfwidth (ms)');
set(gca,'FontSize', 15);
%%
figure;
boxplot(tbl.Blocked_current_level(tbl.Rsquared > 0.9) , tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Depol current level APs get blocked (pA)');
set(gca,'FontSize', 15);
%% 
unique_cell_types = countlabels({results.cell_type});
%% SPIKE PER EPOCH
figure;
tiledlayout('flow');

for i=1:height(unique_cell_types)
nexttile; hold on;
    for j=1:68
        ctype = tbl.Cell_type(j);
        if ctype{1} == unique_cell_types.Label(i)
            spike_n =  tbl.Spike_numbers(j);
            plot(spike_n{1}, "o-");
            ylabel("Spike per epoch")
        end

    end
    hold off;
    title(unique_cell_types.Label(i));
hold off;
end

%% LATENCY
figure;
tiledlayout('flow');

for i=1:height(unique_cell_types)
nexttile; hold on;
    for j=1:68
        ctype = tbl.Cell_type(j);
        if ctype{1} == unique_cell_types.Label(i)
            spike_n =  tbl.Latency_to_first_spike(j);
            plot(spike_n{1}, "o-");
            ylabel("Latency to first spike (ms)")
        end

    end
    hold off;
    title(unique_cell_types.Label(i));
hold off;
end
%% ADAPTATION INDEX
figure;
tiledlayout('flow');

for i=1:height(unique_cell_types)
nexttile; hold on;
    for j=1:68
        ctype = tbl.Cell_type(j);
        if ctype{1} == unique_cell_types.Label(i)
            spike_n =  tbl.Adaptation_index(j);
            plot(spike_n{1}, "o-");
            ylabel("Adaptation Index")
        end

    end
    hold off;
    title(unique_cell_types.Label(i));
hold off;
end

%% ISI CV
figure;
tiledlayout('flow');

for i=1:height(unique_cell_types)
nexttile; hold on;
    for j=1:68
        ctype = tbl.Cell_type(j);
        if ctype{1} == unique_cell_types.Label(i)
            spike_n =  tbl.ISI_cv(j);
            plot(spike_n{1} *100, "o-");
            ylabel("ISI CV %")
        end

    end
    hold off;
    title(unique_cell_types.Label(i));
hold off;
end
%%
figure;
boxplot(tbl.Max_slope(tbl.Rsquared > 0.9) , tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Max slope (mV/mS)');
set(gca,'FontSize', 15);
hold on;
swarmchart(tbl(tbl.Rsquared > 0.9, "Cell_type"), tbl.Max_slope(tbl.Rsquared > 0.9))
%%
two_or_more_cells = unique_cell_types(unique_cell_types.Count > 1, :);
included = zeros(68,1);
for j=1:68
    ctype = tbl.Cell_type(j);
    ctype{1}
    if ismember(ctype{1}, two_or_more_cells.Label)
        included(j) = true
   
    end
end
%%
tbl.Included = included;
newtbl = tbl(tbl.Included == 1 & tbl.Rsquared > 0.9 & ~isnan(tbl.Tau), :)
%%
tbl(tbl.Cell_type{:} == 'Bursty suppressed by contrast', :)