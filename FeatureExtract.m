function [result_struct]=FeatureExtract(result_table)
%% init
start_time = result_table.pre_time_ms * 10^-3 * result_table.sample_rate;
end_time = start_time + result_table.stim_time_ms * 10^-3 *result_table.sample_rate;
hyper_current_epoch = find(result_table.inj_current < 0);
depol_current_epoch = find(result_table.inj_current > 0);
hyper_current_level_pA = result_table.inj_current(hyper_current_epoch);
depol_current_level_pA = result_table.inj_current(depol_current_epoch);
hyper_Vm = result_table.example_traces(hyper_current_epoch,:);
hyper_Vm = hyper_Vm';
depol_Vm = result_table.example_traces(depol_current_epoch,:);
depol_Vm = depol_Vm';
time_in_s = linspace(0, size(hyper_Vm,1), size(hyper_Vm,1)) / result_table.sample_rate;
%%

%%
stable_Vm = mean(hyper_Vm(start_time:end_time,:));
R_linear = fitlm(hyper_current_level_pA, stable_Vm');
resistance = struct();
resistance.R_MOhm = R_linear.Coefficients.Estimate('x1') * 1000; % V/I = mV/pA = 10^9(Giga) => Convert to 10^6 (Mega)Ohm
resistance.AdjustedRsqaured = R_linear.Rsquared.Adjusted;
if R_linear.Rsquared.Adjusted < 0.90
warning('Resistance fit kinda sucky. Check if Ih went brrr')
end


% Calculate Tau (ms)
hyper_epoch_less_than_minus50 = find(hyper_current_level_pA > -50); % INJECTED CURRENT LESS HYPERPOLARIZING THAN -50

ft = fittype('a + b*exp(-x*c)', 'independent', 'x'); %One parameter exp fit with asymt to Vinf
tau_array = zeros(length(hyper_epoch_less_than_minus50),1);
Vinf_array = zeros(length(hyper_current_epoch), 1);
for i=1:length(hyper_epoch_less_than_minus50)
    f = fit([time_in_s(start_time:end_time)]', hyper_Vm(start_time:end_time,...
        hyper_epoch_less_than_minus50(i)), ft, 'StartPoint',[-60,10,30]);


    
    tau_array(i) =  f.c;


end


%Return tau struct
tau = struct();
tau.Tau_ms = mean(1./tau_array)*1000;
tau.SD_ms = std(1./tau_array)*1000;
tau.N = size(tau_array, 2);
if tau.SD_ms > 10
    warning('Tau SD too high')
end

%Return Capacitance
capacitance_pF = tau.Tau_ms / resistance.R_MOhm * 100;

%% Sag
min_Vm = min(hyper_Vm(:,hyper_current_level_pA <= -50), [], 1);
fit_sag_peak_vs_stable = fitlm(min_Vm, stable_Vm(hyper_current_level_pA <= -50));
sag_coeff = fit_sag_peak_vs_stable.Coefficients(2,1);

%%

%% Return values
result_struct.resistance = resistance;
result_struct.tau = tau;
result_struct.capacitance_pF = capacitance_pF;
result_struct.depol_current_level = depol_current_level_pA;
result_struct.hyper_current_level = hyper_current_level_pA;
result_struct.resting_Vm = result_table.vrest;
result_struct.sag_coeff = sag_coeff;
end