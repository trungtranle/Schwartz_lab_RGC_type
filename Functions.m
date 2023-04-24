function [stable_Vm, R_linear] = Functions(hyper_Vm, start_time, end_time, hyper_current_level_pA)
stable_Vm = mean(hyper_Vm(start_time:end_time,:));
figure;
scatter(hyper_current_level_pA, stable_Vm)
hold on;
R_linear = fitlm(hyper_current_level_pA, stable_Vm');
plot(R_linear)
end