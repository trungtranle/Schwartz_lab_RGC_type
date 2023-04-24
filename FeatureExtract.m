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
    stable_Vm = mean(hyper_Vm(start_time:end_time,:));
    R_linear = fitlm(hyper_current_level_pA, stable_Vm');
    resistance = struct();
    %if R_linear.Rsquared.Adjusted > 0.90
        resistance.R_MOhm = R_linear.Coefficients.Estimate('x1') * 1000; % V/I = mV/pA = 10^9(Giga) => Convert to 10^6 (Mega)Ohm
        resistance.AdjustedRsqaured = R_linear.Rsquared.Adjusted;
    %else
        %error('Resistance fit kinda sucky. Check if Ih went brrr')
    %end
    %% Return values
    result_struct.resistance = resistance;
  









end