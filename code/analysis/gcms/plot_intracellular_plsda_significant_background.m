%% Descriptive intracellular PLS-DA after a matched-background screen
% Each bacterial sample carrying suffix _i is corrected against Media_i.
% Metabolites enter the PLS-DA only when their corrected grand mean is
% significantly greater than zero after accounting for strain and batch in
% the balanced design and controlling the false-discovery rate across all 18
% measured metabolites. The PLS-DA is a descriptive full-data visualization,
% not an independently validated classifier.

clear;

scriptPath = mfilename('fullpath');
localRoot = fileparts(fileparts(fileparts(scriptPath)));
repoRoot = fileparts(localRoot);
sourceFile = package_input('intra.xlsx');
tableDir = fullfile(repoRoot, 'results', 'tables', 'gcms', ...
    'gcms_intracellular_matched');
figureDir = fullfile(repoRoot, 'results', 'diagnostic_figures', 'gcms', ...
    'gcms_intracellular_matched');
if ~isfolder(tableDir)
    mkdir(tableDir);
end
if ~isfolder(figureDir)
    mkdir(figureDir);
end

strainOrder = ["ST1-75", "ST1-68", "ST1-6", "VPI10463", "R20291"];
strainColors = [
    36, 104, 162;
    60, 157, 86;
    142, 90, 169;
    217, 71, 63;
    228, 154, 47] ./ 255;
batchMarkers = {'o', 's', '^'};

%% Load and matched-blank correct
source = readtable(sourceFile, 'VariableNamingRule', 'preserve');
sampleID = string(source{:, 1});
metabolite = string(source.Properties.VariableNames(2:end));
metabolite = regexprep(metabolite, '_\d+Results$', '');
rawSignal = double(source{:, 2:end});

group = extractBefore(sampleID, '_');
group(group == "VPI") = "VPI10463";
batch = str2double(extractAfter(sampleID, '_'));

expectedGroups = ["Media", strainOrder];
for g = expectedGroups
    assert(sum(group == g) == 3, 'Expected three labeled rows for %s.', g);
    assert(isequal(sort(batch(group == g))', 1:3), ...
        'Expected batches 1, 2, and 3 for %s.', g);
end

bacterialRow = group ~= "Media";
bacterialIndex = find(bacterialRow);
bacterialSample = sampleID(bacterialRow);
bacterialGroup = group(bacterialRow);
bacterialBatch = batch(bacterialRow);
corrected = nan(sum(bacterialRow), numel(metabolite));
for i = 1:numel(bacterialIndex)
    sourceRow = bacterialIndex(i);
    mediaRow = find(group == "Media" & batch == batch(sourceRow), 1);
    corrected(i, :) = rawSignal(sourceRow, :) - rawSignal(mediaRow, :);
end

%% One-sided, FDR-controlled screen for signal above matched background
% In this balanced 5-strain x 3-batch design, the grand mean is orthogonal to
% strain and batch effects. Interaction/residual variation supplies 8 df.
nStrain = numel(strainOrder);
nBatch = numel(unique(bacterialBatch));
residualDf = size(corrected, 1) - nStrain - nBatch + 1;
grandMean = mean(corrected, 1)';
standardError = nan(numel(metabolite), 1);
tStatistic = nan(numel(metabolite), 1);
oneSidedP = nan(numel(metabolite), 1);

for j = 1:numel(metabolite)
    values = corrected(:, j);
    fitted = repmat(mean(values), size(values));
    for g = 1:nStrain
        rows = bacterialGroup == strainOrder(g);
        fitted(rows) = fitted(rows) + mean(values(rows)) - mean(values);
    end
    for b = unique(bacterialBatch, 'stable')'
        rows = bacterialBatch == b;
        fitted(rows) = fitted(rows) + mean(values(rows)) - mean(values);
    end
    residual = values - fitted;
    residualMeanSquare = sum(residual .^ 2) / residualDf;
    standardError(j) = sqrt(residualMeanSquare / numel(values));
    if standardError(j) == 0
        tStatistic(j) = sign(grandMean(j)) * Inf;
    else
        tStatistic(j) = grandMean(j) / standardError(j);
    end
    oneSidedP(j) = 1 - tcdf(tStatistic(j), residualDf);
end

qValue = bhAdjust(oneSidedP);
include = grandMean > 0 & qValue < 0.05;
assert(sum(include) >= 2, ...
    'Fewer than two metabolites passed the above-background screen.');

screenTable = table(metabolite', grandMean, standardError, tStatistic, ...
    repmat(residualDf, numel(metabolite), 1), oneSidedP, qValue, include, ...
    'VariableNames', {'Metabolite', 'MeanBlankCorrectedSignal', ...
    'StandardErrorAdjustedForStrainAndBatch', 'OneSidedT', 'ResidualDf', ...
    'OneSidedP', 'BenjaminiHochbergQ', 'IncludedInPlsDa'});
writetable(screenTable, fullfile(tableDir, ...
    'intracellular_significantly_above_background_screen.csv'));

%% Two-component descriptive PLS-DA
selectedMetabolite = metabolite(include);
selectedCorrected = corrected(:, include);
selectedZ = (selectedCorrected - mean(selectedCorrected, 1)) ./ ...
    std(selectedCorrected, 0, 1);

classLabel = nan(numel(bacterialGroup), 1);
response = zeros(numel(bacterialGroup), nStrain);
for g = 1:nStrain
    rows = bacterialGroup == strainOrder(g);
    classLabel(rows) = g;
    response(rows, g) = 1;
end

[xLoading, ~, xScore, ~, ~, percentVariance, ~, stats] = ...
    plsregress(selectedZ, response, 2);
assert(size(xScore, 2) >= 2, 'Two PLS components were not returned.');

scoreTable = table(bacterialSample, bacterialGroup, bacterialBatch, ...
    xScore(:, 1), xScore(:, 2), 'VariableNames', ...
    {'Sample', 'Strain', 'Batch', 'PLS1Score', 'PLS2Score'});
writetable(scoreTable, fullfile(tableDir, ...
    'intracellular_descriptive_plsda_2d_scores.csv'));

loadingTable = table(selectedMetabolite', xLoading(:, 1), xLoading(:, 2), ...
    stats.W(:, 1), stats.W(:, 2), 'VariableNames', ...
    {'Metabolite', 'PLS1Loading', 'PLS2Loading', 'PLS1Weight', ...
    'PLS2Weight'});
writetable(loadingTable, fullfile(tableDir, ...
    'intracellular_descriptive_plsda_2d_loadings.csv'));

componentTable = table(["PLS1"; "PLS2"], ...
    percentVariance(1, 1:2)', percentVariance(2, 1:2)', ...
    'VariableNames', {'Component', 'XVarianceFraction', ...
    'YVarianceFraction'});
writetable(componentTable, fullfile(tableDir, ...
    'intracellular_descriptive_plsda_component_variance.csv'));

%% Figure
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1080, 900]);
layout = tiledlayout(figureHandle, 1, 1, 'TileSpacing', 'compact', ...
    'Padding', 'loose');

ax = nexttile(layout);
hold(ax, 'on');
xline(ax, 0, 'Color', [0.88, 0.88, 0.88]);
yline(ax, 0, 'Color', [0.88, 0.88, 0.88]);

% Overlay loadings on the score plane. A single scalar preserves loading
% directions while bringing arrow lengths into the observed score range.
scoreRadius = 0.72 * min(max(abs(xScore(:, 1))), max(abs(xScore(:, 2))));
loadingRadius = max(vecnorm(xLoading(:, 1:2), 2, 2));
loadingScale = scoreRadius / loadingRadius;
scaledLoading = xLoading(:, 1:2) .* loadingScale;
for j = 1:numel(selectedMetabolite)
    quiver(ax, 0, 0, scaledLoading(j, 1), scaledLoading(j, 2), 0, ...
        'Color', [0.49, 0.54, 0.58], 'LineWidth', 1.15, ...
        'MaxHeadSize', 0.20);
    text(ax, 1.08 * scaledLoading(j, 1), ...
        1.08 * scaledLoading(j, 2), selectedMetabolite(j), ...
        'FontName', 'Arial', 'FontSize', 8.5, ...
        'Color', [0.24, 0.28, 0.31], ...
        'HorizontalAlignment', horizontalAlignment(scaledLoading(j, 1)), ...
        'VerticalAlignment', verticalAlignment(scaledLoading(j, 2)));
end

% Plot samples above the loading vectors so that the biological observations
% remain visually primary.
for g = 1:nStrain
    for b = 1:nBatch
        row = classLabel == g & bacterialBatch == b;
        scatter(ax, xScore(row, 1), xScore(row, 2), 95, ...
            strainColors(g, :), batchMarkers{b}, 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
    end
    centroid = [mean(xScore(classLabel == g, 1)), ...
        mean(xScore(classLabel == g, 2))];
    scatter(ax, centroid(1), centroid(2), 240, strainColors(g, :), ...
        'o', 'filled', 'MarkerEdgeColor', [0.15, 0.15, 0.15], ...
        'LineWidth', 1.5);
end
axis(ax, 'equal');
box(ax, 'off');
set(ax, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1);
xlabel(ax, sprintf('PLS1 score (%.1f%% X variance)', ...
    100 * percentVariance(1, 1)));
ylabel(ax, sprintf('PLS2 score (%.1f%% X variance)', ...
    100 * percentVariance(1, 2)));
title(ax, 'Samples and metabolite loadings', ...
    'FontWeight', 'bold');

strainLegend = gobjects(nStrain, 1);
for g = 1:nStrain
    strainLegend(g) = scatter(ax, nan, nan, 80, strainColors(g, :), ...
        'o', 'filled', 'MarkerEdgeColor', 'w');
end
legend(ax, strainLegend, strainOrder, 'Location', 'southoutside', ...
    'NumColumns', nStrain, 'Box', 'off');

sgtitle(layout, sprintf(['Intracellular metabolite pools differ modestly ' ...
    'across strains\nDescriptive PLS-DA using %d metabolites above matched ' ...
    'background; batches 1/2/3 = circles/squares/triangles; arrows = ' ...
    'scaled loadings'], sum(include)), ...
    'FontName', 'Arial', 'FontSize', 15, 'FontWeight', 'bold');

outputFile = fullfile(figureDir, ...
    'intracellular_plsda_2d_significant_above_background.png');
exportgraphics(figureHandle, outputFile, 'Resolution', 300);
close(figureHandle);

fprintf(['Included %d of %d metabolites at BH q<0.05 for a one-sided ' ...
    'grand-mean test above matched background.\n'], ...
    sum(include), numel(metabolite));
fprintf('Included metabolites: %s.\n', strjoin(selectedMetabolite, ', '));
fprintf('PLS1 and PLS2 explain %.2f%% and %.2f%% of X variance.\n', ...
    100 * percentVariance(1, 1), 100 * percentVariance(1, 2));
fprintf('Figure: %s\n', outputFile);

%% Local functions
function adjusted = bhAdjust(pValue)
adjusted = nan(size(pValue));
finite = find(isfinite(pValue));
[sortedP, order] = sort(pValue(finite));
m = numel(sortedP);
sortedQ = sortedP .* m ./ (1:m)';
for i = (m - 1):-1:1
    sortedQ(i) = min(sortedQ(i), sortedQ(i + 1));
end
sortedQ = min(sortedQ, 1);
restored = nan(m, 1);
restored(order) = sortedQ;
adjusted(finite) = restored;
end

function alignment = horizontalAlignment(value)
if value < 0
    alignment = 'right';
else
    alignment = 'left';
end
end

function alignment = verticalAlignment(value)
if value < 0
    alignment = 'top';
else
    alignment = 'bottom';
end
end
