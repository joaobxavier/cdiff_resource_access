%% Read data files
base_dir = fullfile(package_root(),'code');
addpath(genpath(fullfile(base_dir, 'analysis')));
data_dir = fullfile(package_root(),'input','flow');
%% Load all cell data
mice = [1 6 7 11 12 13];

allMice = [];
for i = 1:length(mice)
    allCellsMX = readtable(project_data_file('processed', 'flow_cytometry', ...
        'ks10_adaptive_csv_files', 'all_adaptive', ...
        sprintf('all_adaptive_Specimen_001_%d.csv', mice(i))));
    allCellsMX.mouse(:) = mice(i);
    if i == 1
        allMice = allCellsMX;
    else
        allMice = [allMice; allCellsMX];
    end
end

%% Add gate flags
allMiceWithGates = allMice;
cellTypes = {'B_cells' 'CD19negTCRbneg' 'NK_cells' 'T_cells' 'CD4pos' 'CD8pos'};

for j = 1:length(cellTypes)
    allCellsOfType = [];
    dirName = fullfile(data_dir, sprintf('ks10_adaptive_csv_files/%s/*.csv', cellTypes{j}));
    d = dir(dirName);
    for i = 1:length(d)
        allCellsMX = readtable(fullfile(data_dir, sprintf('ks10_adaptive_csv_files/%s/%s', cellTypes{j}, d(i).name)));

        aux = strsplit(d(i).name, '_001_');
        aux = strsplit(aux{2}, '_');
        mouse = str2double(aux{1});
        allCellsMX.mouse(:) = mouse;
        if i == 1
            allCellsOfType = allCellsMX;
        else
            allCellsOfType = [allCellsOfType; allCellsMX];
        end
    end
    binaryGate = ismember(allMice, allCellsOfType);
    allMiceWithGates = addvars(allMiceWithGates, binaryGate);
    allMiceWithGates.Properties.VariableNames{end} = cellTypes{j};
end

%% Remove events that did not pass any gate
didNotPassAnyGate = sum(allMiceWithGates{:, 17:22}, 2) == 0;
allGatedCell = allMiceWithGates(didNotPassAnyGate == 0, :);

%% Find cells within fluorescence range of gated cells
fluor = allGatedCell{:, 7:14};
allFlowMetrics = allGatedCell{:, 1:14};

maxValuesAmongGated = max(allFlowMetrics);
minValuesAmongGated = min(allFlowMetrics);
allFlowMetricsAllCells = allMiceWithGates{:, 1:14};
cellsWithinRanges = all(allFlowMetricsAllCells >= minValuesAmongGated & ...
    allFlowMetricsAllCells <= maxValuesAmongGated, 2);
allCellsWithinRange = allMiceWithGates(cellsWithinRanges, :);

fluor = allCellsWithinRange{:, 7:14};
allFlowMetrics = allCellsWithinRange{:, 1:14};

%% Define mouse groups
UI_mice = [1, 6, 7];
Avirulent_mice = [11, 12, 13];

mouseGroup = repmat({'Other'}, height(allCellsWithinRange), 1);
mouseGroup(ismember(allCellsWithinRange.mouse, UI_mice)) = {'UI'};
mouseGroup(ismember(allCellsWithinRange.mouse, Avirulent_mice)) = {'Avirulent'};
mouseGroup = categorical(mouseGroup);
prettyTypeNames = {'B cells', 'CD19^{-}TCR\beta^{-}', 'NK cells', ...
    'T cells', 'CD4^{+}', 'CD8^{+}'};

%% Generate categorical cell type labels
cellTypeVars = allCellsWithinRange.Properties.VariableNames(17:22);

cellLabels = repmat({'Other'}, height(allCellsWithinRange), 1);
for i = 1:length(cellTypeVars)
    idx = allCellsWithinRange{:, cellTypeVars{i}} == 1;
    cellLabels(idx) = cellTypeVars(i);
end
cellLabels = categorical(cellLabels);

%% Compute cell type fractions per mouse
cellTypeNames = allCellsWithinRange.Properties.VariableNames(17:22);
mouseIDs = unique(allCellsWithinRange.mouse);
countCells = table(mouseIDs, 'VariableNames', {'mouse'});

for i = 1:length(mouseIDs)
    mouse = mouseIDs(i);
    dataMouse = allCellsWithinRange(allCellsWithinRange.mouse == mouse, :);
    totalEvents = height(dataMouse);
    fractions = zeros(1, length(cellTypeNames));
    for j = 1:length(cellTypeNames)
        fractions(j) = sum(dataMouse{:, cellTypeNames{j}}) / totalEvents;
    end
    countCells{i, 2:(length(cellTypeNames)+1)} = fractions;
end
countCells.Properties.VariableNames(2:end) = cellTypeNames;

%% Wilcoxon rank-sum tests
UI = [1 6 7];
Avirulent = [11 12 13];

data = countCells{:, 2:end};
mouseIDs = countCells.mouse;

UI_data = data(ismember(mouseIDs, UI), :);
Avirulent_data = data(ismember(mouseIDs, Avirulent), :);

UI_mean = mean(UI_data, 1);
Avirulent_mean = mean(Avirulent_data, 1);
UI_std = std(UI_data, 0, 1);
Avirulent_std = std(Avirulent_data, 0, 1);

groupMeans = [UI_mean; Avirulent_mean];
nCellTypes = length(cellTypeNames);


pvals = nan(1, nCellTypes);
for j = 1:nCellTypes
    pvals(j) = ranksum(UI_data(:,j), Avirulent_data(:,j));
end

fprintf('\nWilcoxon rank-sum p-values (n=3 per group):\n');
for j = 1:nCellTypes
    fprintf('  %s: p = %.4f\n', cellTypeNames{j}, pvals(j));
end


out=getappdata(0,'CdiffVectorOutput');
save(fullfile(out,'flow_adaptive.mat'),'UI_data','Avirulent_data','prettyTypeNames','countCells','pvals');
writetable(countCells,fullfile(package_root(),'results','tables','main','flow_adaptive_fractions.csv'));
