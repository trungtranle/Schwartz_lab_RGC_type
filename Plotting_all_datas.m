%%
load April21_63cells_68datasets.mat
%%
R_array = zeros(68,1);
fit_array = zeros(68,1);
for i = 1:68;
    R = FeatureExtract(results(i));
    R_array(i) = R.resistance.R_MOhm;
    fit_array(i) = R.resistance.AdjustedRsqaured;
end

tbl = table(R_array, fit_array, {results.cell_type}', 'VariableNames', {'Resistance', 'Rsquared', 'Cell_type'})
figure;
boxplot(tbl.Resistance(tbl.Rsquared > 0.9), tbl.Cell_type(tbl.Rsquared > 0.9));
xlabel('Cell types');
ylabel('Input Resistance (M\Omega)');



