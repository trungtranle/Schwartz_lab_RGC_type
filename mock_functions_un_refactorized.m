% Provided that the result_table was already loaded
result_table = results(2)

%% Set things up
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
%% Plot hyper epoch
figure;
hold on;
plot(time_in_s, hyper_Vm)
legend(string(hyper_current_level_pA))
xlabel('Time (s)')
ylabel('V_m (mV)')
hold off;

%% Plot depol epoch
figure;
hold on;
plot(time_in_s, depol_Vm)
legend(string(depol_current_level_pA))
xlabel('Time (s)')
ylabel('V_m (mV)')
hold off;

%% Calculate R, Tau, C
% Calculate Resistance (MOhm)
stable_Vm = mean(hyper_Vm(start_time:end_time,:));
figure;
scatter(hyper_current_level_pA, stable_Vm)
hold on;
R_linear = fitlm(hyper_current_level_pA, stable_Vm');
plot(R_linear)

%return Resistance struct
resistance = struct();

    resistance.R_MOhm = R_linear.Coefficients.Estimate('x1') * 1000; % V/I = mV/pA = 10^9(Giga) => Convert to 10^6 (Mega)Ohm
    resistance.AdjustedRsqaured = R_linear.Rsquared.Adjusted;
if R_linear.Rsquared.Adjusted > 0.90
    warning('Resistance fit kinda sucky. Check if Ih went brrr')
end
hold off;
% Calculate Tau (ms)
hyper_epoch_less_than_minus50 = find(hyper_current_level_pA > -50); % INJECTED CURRENT LESS HYPERPOLARIZING THAN -50
figure; hold on;
%plot(time_in_s(start_time:end_time), hyper_Vm(start_time:end_time, hyper_epoch_less_than_minus50))
ft = fittype('a + b*exp(-x*c)', 'independent', 'x'); %One parameter exp fit with asymt to Vinf
tau_array = zeros(length(hyper_epoch_less_than_minus50),1);
Vinf_array = zeros(length(hyper_current_epoch), 1);
for i=1:length(hyper_epoch_less_than_minus50)
    f = fit([time_in_s(start_time:end_time)]', hyper_Vm(start_time:end_time,...
        hyper_epoch_less_than_minus50(i)), ft, 'StartPoint',[-60,10,30]);

    plot(f, time_in_s(start_time:end_time), hyper_Vm(start_time:end_time, hyper_epoch_less_than_minus50(i)))

    tau_array(i) =  f.c;
    f.c


end
hold off;

%Return tau struct
tau = struct();
tau.Tau_ms = mean(1./tau_array)*1000;
tau.SD_ms = std(1./tau_array)*1000;
tau.N = size(tau_array, 2);
if tau.SD_ms > 10
    error('Tau SD too high')
end

%Return Capacitance
capacitance_pF = tau.Tau_ms / resistance.R_MOhm * 100;

%% Correlate Calculated Vinf from model to Vinf from averaging (this part does matter anyway)
% figure;
% hold on;
% scatter(Vinf_array, stable_Vm)
% hold off;
% figure();
% hold on;
% lm_test = fitlm(Vinf_array, stable_Vm);
% scatter(Vinf_array, stable_Vm);
% plot(lm_test)
% xlabel('Calculated V_{inf} from single-term exponential fit (mV)', 'Interpreter','tex');
% ylabel("Stable V_{inf} from averaging (mV)", 'Interpreter','tex');
% title('')
% text(mean(Vinf_array), mean(stable_Vm)-20, ['Adjusted R^2 = ' num2str(lm_test.Rsquared.Adjusted)])
% hold off;
% % most cells have R_squared = 1. Check more cells.

%% Calculate Sag ANGLE
min_Vm = min(hyper_Vm(:,hyper_current_level_pA <= -50), [], 1);
figure; hold on;
%plot(min_Vm, stable_Vm(find(hyper_current_level_pA <= -50)))
fit_sag_peak_vs_stable = fitlm(min_Vm, stable_Vm(hyper_current_level_pA <= -50))
plot(fit_sag_peak_vs_stable)
plot(stable_Vm(hyper_current_level_pA <= -50), stable_Vm(hyper_current_level_pA <= -50) , 'g')
hold off;

%%

test_spike = depol_Vm(:,6);

diff_1 = diff(test_spike(:)) * (result_table.sample_rate/1000);
diff_2 = diff(test_spike(:), 2) * (result_table.sample_rate/1000);

figure;
ax1 = subplot(3,1,1);
plot(test_spike(:));
ax2 = subplot(3,1,2);
plot(diff_1);
yline(0.2*max(diff(test_spike())) * (result_table.sample_rate/1000));
ax3 = subplot(3,1,3);
plot(diff_2);
linkaxes([ax1, ax2, ax3], "x"); hold off;

figure;

[pks, locs, w, p] = findpeaks(test_spike,50000, 'MinPeakProminence', 6, 'MinPeakHeight',-10, 'MinPeakDistance', 0.001);
findpeaks(test_spike,50000, 'MinPeakProminence', 6, 'MinPeakHeight',-5, 'MinPeakDistance', 0.001)
hold off;

n_pks_array = zeros(size(depol_current_level_pA,1), 1);

figure;
hold on
for i=[1:size(depol_current_level_pA, 1)]

    test_spike = depol_Vm(:,i);

    peks = findpeaks(test_spike(start_time : end_time),50000, 'MinPeakProminence',6, 'MinPeakHeight', -10, 'MinPeakDistance', 0.001);

    n_pks_array(i) = size(peks, 1);
end

hold off;



%% Does it spike spontaneously
spontaneous_peak_array = zeros(size(depol_current_level_pA,1), 1);


figure; hold on;
for i=1:size(depol_current_level_pA, 1)
    spontaneous_peak_array(i) = size(findpeaks(depol_Vm(1:start_time, i), ...
        "MinPeakProminence", 6, "MinPeakHeight", -10, "MinPeakDistance", result_table.sample_rate*1e-3),1);
    findpeaks(depol_Vm(1:start_time, i),  ...
       "MinPeakProminence", 6, "MinPeakHeight", -10, "MinPeakDistance", result_table.sample_rate*1e-3)
end

spontaneous_spike = (mean(spontaneous_peak_array))/(start_time/result_table.sample_rate); %Hz

%% Find first spike
first_spike = [0 0 0]; %peak loc epoch
if spontaneous_spike ~= 0
    end_time_find = start_time;
    start_time_find = 1;
else 
    end_time_find = end_time;
    start_time_find = start_time;
end
figure;
i = 1;
while first_spike(2) == 0
    findpeaks(depol_Vm(start_time_find:end_time_find, i), ...
       'MinPeakProminence', 6, 'MinPeakHeight', -10, "MinPeakDistance", result_table.sample_rate*1e-3); %peak separation at least 1 ms
    [pks, locs] = findpeaks(depol_Vm(start_time_find:end_time_find, i), ...
        'MinPeakProminence', 6, 'MinPeakHeight', -10, "MinPeakDistance", result_table.sample_rate*1e-3);
    
 
    try
    first_spike = [pks(1) locs(1) i];
    
    end
    i = i+1;

end
hold on;

if spontaneous_spike == 0
    first_spike(2) = start_time_find + first_spike(2);

end
figure; hold on;
findpeaks(depol_Vm(start_time_find:end_time_find, first_spike(3)), ...
        'MinPeakProminence', 6, 'MinPeakHeight', -10, "MinPeakDistance", result_table.sample_rate*1e-3); %peak separation at least 1 ms

Vm_diff_1 = diff(depol_Vm(first_spike(2):(first_spike(2) + (10 * 1e-3 * result_table.sample_rate)), first_spike(3)),1); %take from 10ms
locations = find(Vm_diff_1 == 0);
trough = [0 0 0];
trough(2) = first_spike(2) + locations(1); % trough location;
trough(3) = first_spike(3); %trough epoch;
trough(1) = depol_Vm(trough(2), trough(3)); %trough level mV


plot(trough(2), trough(1), "x", "LineWidth",100)

hold off;

start_time_for_threshold = max(1, (first_spike(2) - (5 * 1e-3 * result_table.sample_rate))); %take from -5ms or the start

Vm_diff_2 = diff(depol_Vm(start_time_for_threshold : first_spike(2), first_spike(3)), 1); 

plot(Vm_diff_2); hold on;

threshold_loc = find(Vm_diff_2 >= max(Vm_diff_2)*0.2,1);
V_threshold = depol_Vm(threshold_loc + 1, first_spike(3)); % + 1 bc diff lost one position

half_height = (first_spike(1) + trough(1)) / 2;
half_width_times = sum(depol_Vm(threshold_loc : trough(2), trough(3)) >= half_height) /result_table.sample_rate * 1e3; %ms



%% spikes and ISIs during current injections

spike_numbers = zeros(length(depol_current_epoch), 1);
latency_to_first_spike = zeros(length(depol_current_epoch), 1);
adaptation_index = zeros(length(depol_current_epoch), 1);
ISI_cv = zeros(length(depol_current_epoch), 1);
blocked = zeros(length(depol_current_epoch), 1);

%figure; hold on;
for epoch=1:length(depol_current_epoch)
    
    %findpeaks(depol_Vm(start_time:end_time, epoch), ...
    %    'MinPeakProminence', 6, 'MinPeakHeight', -5, "MinPeakDistance", result_table.sample_rate*1e-3); 
    [spikes, locs] = findpeaks(depol_Vm(start_time:end_time, epoch), ...
        'MinPeakProminence', 6, 'MinPeakHeight', -5, "MinPeakDistance", result_table.sample_rate*1e-3); %peak separation at least 1 ms
    try
        latency_to_first_spike(epoch) = locs(1) - start_time / result_table.sample_rate*1e-3;
        spike_numbers(epoch) = length(spikes);
        locs_move_1 = [0; locs(1:end-1)];
        ISIs = locs - locs_move_1; ISIs = ISIs(2:end);
        ISI_cv(epoch) = std(ISIs) / mean(ISIs);
        for j=1:(length(ISIs)-1)
           
            adaptation = (ISIs(j+1) - ISIs(j))/ (ISIs(j+1) + ISIs(j));
            adaptation_index(epoch) = adaptation_index(epoch) + adaptation;
            
        end
        adaptation_index(epoch) = adaptation_index(epoch) / (length(ISIs) - 1);
        if locs(end) < (end_time - start_time)/2
            blocked(epoch) = true;
        else
            blocked(epoch) = false;
        end
        
    catch
        latency_to_first_spike(epoch) = end_time - start_time / result_table.sample_rate * 1e-3;
        spike_numbers(epoch) = length(spikes);
    end
end

blocked_epoch = find(blocked == true, 1);
if isempty(blocked_epoch)
    blocked_current_level = 9999; % not blocked yet
else
    blocked_current_level = depol_current_level_pA(blocked_epoch(1));
end


%interpolation
interpolation_range = linspace(depol_current_level_pA(1), depol_current_level_pA(end), 10)
spike_numbers = interp1(depol_current_level_pA, spike_numbers,interpolation_range)
latency_to_first_spike = interp1(depol_current_level_pA, latency_to_first_spike, interpolation_range)
adaptation_index = interp1(depol_current_level_pA, adaptation_index, interpolation_range)
ISI_cv = interp1(depol_current_level_pA, ISI_cv, interpolation_range)