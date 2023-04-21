%% 
sln_cell.CellEvent * sln_cell.AssignType * sln_cell.Cell & proj(data_group)

%%
sln_symphony.ExperimentCell * sln_cell.CellEvent * sln_cell.AssignType * sln_cell.Cell & proj(data_group)

t = proj(sln_symphony.ExperimentCell * sln_cell.CellEvent * sln_cell.AssignType ...
    * sln_cell.Cell & proj(data_group) & "cell_class" == {'RGC'},'file_name' ,'source_id', 'cell_unid', 'cell_type');

results = fetch(sln_results.DatasetMultiPulsevaryCurrent * t, '*');

% 63 cells, 68 datasets.

%%
