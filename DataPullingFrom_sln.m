%% 
sln_cell.CellEvent * sln_cell.AssignType * sln_cell.Cell & proj(data_group)

%%
sln_symphony.ExperimentCell * sln_cell.CellEvent * sln_cell.AssignType * sln_cell.Cell & proj(data_group)

t = proj(sln_symphony.ExperimentCell * sln_cell.CellEvent * sln_cell.AssignType ...
    * sln_cell.Cell & proj(data_group) & "cell_class" == {'RGC'},'file_name' ,'source_id', 'cell_unid', 'cell_type');

results = fetch(sln_results.DatasetMultiPulsevaryCurrent * t, '*');

% 63 cells, 68 datasets.

%% Datacheck
celltypes = {results.cell_type};
tbl_celltype = countlabels(celltypes)
figure();
bar(tbl_celltype.Label, tbl_celltype.Count);
xlabel('Confirmed cell type');
ylabel('Cells');
text(1, 17, "Total = 63 cells, 68 datasets");
set(gca, 'FontSize', 15, 'TickDir', 'out');
set(gcf, "Color", [1 1 1]);
box off;

