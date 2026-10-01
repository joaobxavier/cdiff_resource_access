%% Generate Figure 2: disease-ranked mouse-screen phenotypes
% Shows the complete 21-strain mouse screen as disease-ranked estimates,
% aligns protection scores to the same strain order, and then displays the
% panel-wide association in a square panel with the biologically informative
% ST1-75, ST1-6, ST1-68 and ST1-49 comparisons labeled. The primary composite
% mouse-score analysis is reused without changing its scoring or models.

clear;

scriptDir = fileparts(mfilename('fullpath'));
localRoot = fileparts(fileparts(scriptDir));
repoRoot = fileparts(localRoot);
addpath(genpath(fullfile(localRoot, 'analysis')));
sourceDir = fullfile(repoRoot, 'results', 'tables', 'main');
figureDir = fullfile(repoRoot, 'results', 'figures');
tableDir = fullfile(repoRoot, 'results', 'tables', 'main');
if ~isfolder(figureDir); mkdir(figureDir); end
if ~isfolder(tableDir); mkdir(tableDir); end

scoreFile = fullfile(repoRoot,'results','tables','mouse_weight_survival_reanalysis','primary_terminal_event_zero_scores.csv');
weightFile = project_data_file('processed', 'mouse', 'scores', ...
    'ProtectionScreen_CDI_mouse.csv');
cfuFile = fullfile(repoRoot, 'results', 'tables', 'main', ...
    'figureS1_fecal_cfu_observations.csv');
assert(all(isfile(string({scoreFile, weightFile, cfuFile}))), ...
    'One or more verified source tables for Figure 2 are unavailable.');

scores = readtable(scoreFile, 'TextType', 'string');
assert(height(scores) == 21, 'Expected 21 strains in the mouse screen.');
requiredVariables = ["Strain", "Mono_Disease_Effect", ...
    "Mono_Disease_CI_Lower", "Mono_Disease_CI_Upper", ...
    "Protection_Effect", "Protection_CI_Lower", ...
    "Protection_CI_Upper"];
assert(all(ismember(requiredVariables, string(scores.Properties.VariableNames))), ...
    'Mouse-score table is missing a required estimate or interval column.');
assert(numel(unique(scores.Strain)) == height(scores), ...
    'Expected one mouse-score row per strain.');

weights = readtable(weightFile, 'TextType', 'string');
weights.experiment = string(weights.experiment);
weights.tx = string(weights.tx);
weights.cdiffstrain = string(weights.cdiffstrain);
weights.mouse = string(weights.mouse);
weights.animal_id = weights.experiment + "|" + weights.tx + "|" + ...
    weights.cdiffstrain + "|" + weights.mouse;
weights = weights(weights.experiment ~= "ks65", :);
assert(all(weights.death(ismissing(weights.relweight)) == 1), ...
    'Focused screen weights are missing outside terminal-event rows.');
st16Weights = selectFocusedExperiment(weights, "st1.6", [3, 3, 3]);
st168Weights = selectFocusedExperiment(weights, "st1.68", [3, 3, 4]);

cfu = readtable(cfuFile, 'TextType', 'string');
focalCfuStrains = ["ST1-75", "ST1-6", "ST1-68", "ST1-49"];
focalCfu = cfu(cfu.Context == "Mono-colonization" & ...
    cfu.Phase == "Late" & ismember(cfu.Candidate, focalCfuStrains), :);
expectedScheduled = [57, 11, 11, 8];
expectedMeasured = [52, 11, 11, 8];
expectedDetected = [47, 11, 11, 8];
for i = 1:numel(focalCfuStrains)
    one = focalCfu(focalCfu.Candidate == focalCfuStrains(i), :);
    assert(height(one) == expectedScheduled(i) && ...
        sum(~ismissing(one.cdiffcfu)) == expectedMeasured(i) && ...
        sum(one.Detected == 1) == expectedDetected(i), ...
        'Unexpected late mono-colonization CFU counts for %s.', ...
        focalCfuStrains(i));
end

rankedScores = sortrows(scores, 'Mono_Disease_Effect', 'descend');
rankedScores.Disease_Rank = (1:height(rankedScores))';
assert(all(diff(rankedScores.Mono_Disease_Effect) <= 0), ...
    'Disease ranking is not descending.');

% Preserve the shared score-table export used by supplementary workflows and
% add an explicitly ranked figure-facing table for this visualization.
writetable(scores, fullfile(tableDir, 'figure2_candidate_scores.csv'));
writetable(rankedScores, fullfile(tableDir, ...
    'figure2_ranked_candidate_scores.csv'));
writetable(st16Weights, fullfile(tableDir, ...
    'figure2_st16_focused_weight_observations.csv'));
writetable(st168Weights, fullfile(tableDir, ...
    'figure2_st168_focused_weight_observations.csv'));
writetable(focalCfu, fullfile(tableDir, ...
    'figure2_focal_monocolonization_cfu.csv'));
focusedSummary = [summarizeFocusedGroups(st16Weights, "ST1-6"); ...
    summarizeFocusedGroups(st168Weights, "ST1-68")];
writetable(focusedSummary, fullfile(tableDir, ...
    'figure2_focused_trajectory_summary.csv'));

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.2, 0.2, 12.5, 12.0]);
% Disease-ranked screen with aligned strain highlights and CFU summaries.

axA = axes(fig, 'Position', [0.065, 0.81, 0.92, 0.13]);
plotRankedBars(axA, rankedScores, 'Disease', false);
panelLabel(axA, 'A');

axB = axes(fig, 'Position', [0.065, 0.63, 0.92, 0.13]);
plotRankedBars(axB, rankedScores, 'Protection', true);
panelLabel(axB, 'B');
linkaxes([axA, axB], 'x');

axC = axes(fig, 'Position', [0.065, 0.075, 0.425, 0.45]);
[rho, pValue] = plotScoreAssociation(axC, scores);
pbaspect(axC, [1, 1, 1]);
panelLabel(axC, 'C');

axD = axes(fig, 'Position', [0.56, 0.408, 0.425, 0.13]);
plotFocusedTrajectories(axD, st16Weights, "st1.6", "ST1-6", ...
    [142, 90, 169] / 255, false);
panelLabel(axD, 'D');

axE = axes(fig, 'Position', [0.56, 0.244, 0.425, 0.13]);
plotFocusedTrajectories(axE, st168Weights, "st1.68", "ST1-68", ...
    [213, 94, 0] / 255, true);
panelLabel(axE, 'E');
linkaxes([axD, axE], 'xy');

axF = axes(fig, 'Position', [0.56, 0.03, 0.425, 0.13]);
plotFocalCfu(axF, focalCfu);
panelLabel(axF, 'F');

outputFile = fullfile(figureDir, 'figure2_ranked_mouse_screen.png');
export_vector_asset(fig, outputFile);
exportgraphics(fig, outputFile, 'Resolution', 300);
close(fig);



Source = ["Terminal-event-zero mouse-score table"; ...
    "Panel co-infection mouse table"; ...
    "Audited fecal-CFU observation table"; ...
    "Figure 2 MATLAB generator"];
Path = string({scoreFile; weightFile; cfuFile; ...
    [mfilename('fullpath'), '.m']});
Use = ["Twenty-one strain disease and protection estimates with 95% confidence intervals"; ...
    "Contemporaneous ST1-6 and ST1-68 strain-alone, co-infection, and VPI10463 trajectories"; ...
    "Late strain-alone fecal CFU measurements for ST1-75, ST1-6, ST1-68, and ST1-49"; ...
    "Disease ranking, aligned protection display, association, focused trajectories, and box-and-whisker CFU panel"];
writetable(table(Source, Path, Use), fullfile(tableDir, ...
    'figure2_source_provenance.csv'));

summaryFile = fullfile(tableDir, 'figure2_analysis_summary.txt');
fid = fopen(summaryFile, 'w');
assert(fid >= 0, 'Unable to create Figure 2 summary.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Figure 2: disease-ranked mouse-screen phenotypes\n');
fprintf(fid, 'Disease-ranked strain order:\n%s\n', ...
    strjoin(rankedScores.Strain, ', '));
fprintf(fid, 'Disease/protection Spearman rho %.12g, p %.12g, n=%d.\n', ...
    rho, pValue, height(scores));
for strain = ["ST1-75", "ST1-6", "ST1-68", "ST1-49"]
    row = scores(scores.Strain == strain, :);
    assert(height(row) == 1, 'Expected one mouse-score row for %s.', strain);
    rank = rankedScores.Disease_Rank(rankedScores.Strain == strain);
    fprintf(fid, ['%s disease rank %d; disease %.12g [%.12g, %.12g]; ' ...
        'protection %.12g [%.12g, %.12g].\n'], strain, rank, ...
        row.Mono_Disease_Effect, row.Mono_Disease_CI_Lower, ...
        row.Mono_Disease_CI_Upper, row.Protection_Effect, ...
        row.Protection_CI_Lower, row.Protection_CI_Upper);
end
fprintf(fid, ['ST1-68 and ST1-49 are emphasized as biological exceptions; ' ...
    'this is a narrative designation, not a statistical outlier test.\n']);
for i = 1:height(focusedSummary)
    fprintf(fid, '%s %s (%s): %d mice, %d terminal events.\n', ...
        focusedSummary.Strain(i), focusedSummary.Condition(i), ...
        focusedSummary.Experiment(i), focusedSummary.N_Animals(i), ...
        focusedSummary.N_Terminal_Events(i));
end
for i = 1:numel(focalCfuStrains)
    one = focalCfu(focalCfu.Candidate == focalCfuStrains(i), :);
    measured = ~ismissing(one.cdiffcfu);
    positive = measured & one.Detected == 1;
    positiveQuartiles = prctile(one.cdiffcfu(positive), [25, 50, 75]);
    fprintf(fid, ['%s late strain-alone fecal CFU: %d/%d measured ' ...
        'samples detected; positive median %.6g CFU/g; ' ...
        'positive IQR %.6g to %.6g; positive range %.6g to %.6g.\n'], ...
        focalCfuStrains(i), ...
        sum(positive), sum(measured), ...
        positiveQuartiles(2), positiveQuartiles(1), positiveQuartiles(3), ...
        min(one.cdiffcfu(positive)), max(one.cdiffcfu(positive)));
end
fprintf(fid, 'Figure: %s\n', outputFile);
disp(fileread(summaryFile));

function focused = selectFocusedExperiment(weights, strainCode, expectedN)
mixedCode = strainCode + ".vpi";
experiment = unique(weights.experiment(weights.cdiffstrain == mixedCode));
assert(isscalar(experiment), ...
    'Expected one focused co-infection experiment for %s.', strainCode);
conditions = [strainCode, mixedCode, "vpi"];
focused = weights(weights.experiment == experiment & ...
    ismember(weights.cdiffstrain, conditions), :);
for i = 1:numel(conditions)
    nAnimals = numel(unique(focused.animal_id( ...
        focused.cdiffstrain == conditions(i))));
    assert(nAnimals == expectedN(i), ...
        'Unexpected %s sample size in experiment %s.', ...
        conditions(i), experiment);
end
assert(all(focused.day >= 0 & focused.day <= 7), ...
    'Focused trajectory day lies outside the expected 0-7 interval.');
end

function summary = summarizeFocusedGroups(raw, strainLabel)
candidateCodes = unique(raw.cdiffstrain( ...
    startsWith(raw.cdiffstrain, "st1.") & ...
    ~endsWith(raw.cdiffstrain, ".vpi")));
assert(isscalar(candidateCodes), ...
    'Expected one strain-alone condition in focused trajectory table.');
codes = [candidateCodes, candidateCodes + ".vpi", "vpi"];
Condition = ["Strain alone"; "Strain + VPI10463"; "VPI10463 alone"];
Strain = repmat(strainLabel, 3, 1);
Experiment = repmat(unique(raw.experiment), 3, 1);
N_Animals = zeros(3, 1);
N_Terminal_Events = zeros(3, 1);
for i = 1:3
    one = raw(raw.cdiffstrain == codes(i), :);
    N_Animals(i) = numel(unique(one.animal_id));
    N_Terminal_Events(i) = sum(one.death == 1);
end
summary = table(Strain, Experiment, Condition, N_Animals, ...
    N_Terminal_Events);
end

function plotRankedBars(ax, rankedScores, phenotype, showStrainLabels)
colors = strainColors(rankedScores.Strain);
x = (1:height(rankedScores))';
hold(ax, 'on');

switch phenotype
    case 'Protection'
        estimate = rankedScores.Protection_Effect;
        lower = rankedScores.Protection_CI_Lower;
        upper = rankedScores.Protection_CI_Upper;
        titleText = 'Protection (in the same strain order)';
        yLabelText = {'Co-infection protection score'; '(+ = better)'};
    case 'Disease'
        estimate = rankedScores.Mono_Disease_Effect;
        lower = rankedScores.Mono_Disease_CI_Lower;
        upper = rankedScores.Mono_Disease_CI_Upper;
        titleText = 'Disease potential';
        yLabelText = {'Mono-colonization disease score'; '(+ = worse)'};
    otherwise
        error('Unknown phenotype: %s', phenotype);
end

bars = bar(ax, x, estimate, 0.76, 'FaceColor', 'flat', ...
    'EdgeColor', 'none');
bars.CData = colors;
errorbar(ax, x, estimate, estimate - lower, upper - estimate, ...
    'LineStyle', 'none', 'Color', [0.20, 0.20, 0.20], ...
    'LineWidth', 0.75, 'CapSize', 2.5);
yline(ax, 0, '-', 'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');

interval = [lower; upper; 0];
padding = max(1.5, 0.08 * range(interval));
ylim(ax, [min(interval) - padding, max(interval) + padding]);
xlim(ax, [0.25, height(rankedScores) + 0.75]);
xticks(ax, x);
if showStrainLabels
    xticklabels(ax, rankedScores.Strain);
    xtickangle(ax, 90);
else
    xticklabels(ax, []);
end
ylabel(ax, yLabelText);
title(ax, titleText, 'FontWeight', 'bold');
set(ax, 'FontName', 'Arial', 'FontSize', 8.0, 'Box', 'off', ...
    'TickDir', 'out', 'Layer', 'top');
end

function [rho, pValue] = plotScoreAssociation(ax, scores)
x = scores.Mono_Disease_Effect;
y = scores.Protection_Effect;
[rho, pValue] = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
hold(ax, 'on');

% Draw model-derived intervals before points so every estimate remains visible.
for i = 1:height(scores)
    plot(ax, [scores.Mono_Disease_CI_Lower(i), ...
        scores.Mono_Disease_CI_Upper(i)], [y(i), y(i)], '-', ...
        'Color', [0.86, 0.86, 0.86], 'LineWidth', 0.75, ...
        'HandleVisibility', 'off');
    plot(ax, [x(i), x(i)], [scores.Protection_CI_Lower(i), ...
        scores.Protection_CI_Upper(i)], '-', ...
        'Color', [0.86, 0.86, 0.86], 'LineWidth', 0.75, ...
        'HandleVisibility', 'off');
end

% The fitted line is a descriptive linear guide; the reported association is
% the rank-based Spearman correlation.
coefficients = polyfit(x, y, 1);
xLine = linspace(min(x), max(x), 200);
plot(ax, xLine, polyval(coefficients, xLine), '--', ...
    'Color', [0.30, 0.30, 0.30], 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
xline(ax, 0, ':', 'Color', [0.55, 0.55, 0.55], ...
    'HandleVisibility', 'off');
yline(ax, 0, ':', 'Color', [0.55, 0.55, 0.55], ...
    'HandleVisibility', 'off');

colors = strainColors(scores.Strain);
scatter(ax, x, y, 52, colors, 'filled', 'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.6);

labelStrains = ["ST1-75", "ST1-6", "ST1-68", "ST1-49"];
labelText = ["ST1-75", "ST1-6", ...
    "ST1-68: low disease, no protection", ...
    "ST1-49: protection despite disease potential"];
offset = [0.35, -0.3; -2.4, 0.7; 0.45, -1.1; 0.45, 0.8];
markers = {'o', '^', 's', 'd'};
for i = 1:numel(labelStrains)
    row = find(scores.Strain == labelStrains(i), 1);
    scatter(ax, x(row), y(row), 92, colors(row, :), 'filled', ...
        'Marker', markers{i}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    text(ax, x(row) + offset(i, 1), y(row) + offset(i, 2), ...
        labelText(i), 'FontSize', 8.0, 'FontWeight', 'bold', ...
        'Color', colors(row, :), 'BackgroundColor', 'w', ...
        'Margin', 1.2, 'Clipping', 'on');
end

xValues = [scores.Mono_Disease_CI_Lower; ...
    scores.Mono_Disease_CI_Upper; 0];
yValues = [scores.Protection_CI_Lower; ...
    scores.Protection_CI_Upper; 0];
xPadding = max(1, 0.07 * range(xValues));
yPadding = max(1, 0.07 * range(yValues));
xlim(ax, [min(xValues) - xPadding, max(xValues) + xPadding]);
ylim(ax, [min(yValues) - yPadding, max(yValues) + yPadding]);
xlabel(ax, 'Mono-colonization disease score (+ = worse)');
ylabel(ax, 'Co-infection protection score (+ = better)');
title(ax, 'Disease potential and protection', ...
    'FontWeight', 'bold');
text(ax, 0.985, 0.96, sprintf('Spearman \\rho = %.3f, p = %.3g', ...
    rho, pValue), 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontSize', 9, 'BackgroundColor', 'w', ...
    'Margin', 2);
set(ax, 'FontName', 'Arial', 'FontSize', 9, 'Box', 'off', ...
    'TickDir', 'out', 'Layer', 'top');
end

function plotFocusedTrajectories(ax, raw, strainCode, strainLabel, ...
        strainColor, showXLabels)
conditionCodes = [strainCode, strainCode + ".vpi", "vpi"];
conditionLabels = ["Strain alone", "Strain + VPI10463", ...
    "VPI10463 alone"];
groupColors = [strainColor; strainColor; 0.27, 0.27, 0.27];
lineStyles = {'--', '-', '-'};
days = sort(unique(raw.day));
meanHandles = gobjects(3, 1);
hold(ax, 'on');
for g = 1:3
    group = raw(raw.cdiffstrain == conditionCodes(g), :);
    ids = unique(group.animal_id, 'stable');
    paleFraction = 0.80;
    if g == 2
        paleFraction = 0.70;
    end
    pale = paleFraction * [1, 1, 1] + ...
        (1 - paleFraction) * groupColors(g, :);
    for i = 1:numel(ids)
        one = sortrows(group(group.animal_id == ids(i), :), 'day');
        observed = ~ismissing(one.relweight);
        plot(ax, one.day(observed), one.relweight(observed), '-', ...
            'Color', pale, 'LineWidth', 0.75, ...
            'HandleVisibility', 'off');
        if any(one.death == 1)
            lastObserved = find(observed, 1, 'last');
            eventColor = strainColor;
            if g == 3; eventColor = [0, 0, 0]; end
            plot(ax, one.day(lastObserved), one.relweight(lastObserved), ...
                'x', 'Color', eventColor, 'MarkerSize', 6.5, ...
                'LineWidth', 1.4, 'HandleVisibility', 'off');
        end
    end
    displayedMean = nan(numel(days), 1);
    for d = 1:numel(days)
        displayedMean(d) = mean(group.relweight(group.day == days(d)), ...
            'omitnan');
    end
    meanHandles(g) = plot(ax, days, displayedMean, lineStyles{g}, ...
        'Color', groupColors(g, :), 'LineWidth', 2.3, ...
        'DisplayName', conditionLabels(g));
end
yline(ax, 100, ':', 'Color', [0.58, 0.58, 0.58], ...
    'HandleVisibility', 'off');
xlim(ax, [0, 7]);
ylim(ax, [70, 115]);
xticks(ax, 0:7);
if showXLabels
    xlabel(ax, 'Day after simultaneous challenge');
else
    xticklabels(ax, []);
end
ylabel(ax, 'Body weight (% of baseline)');
if strainLabel == "ST1-6"
    title(ax, 'ST1-6 protects during co-infection', 'FontWeight', 'bold');
    legend(ax, meanHandles, 'Location', 'southwest', 'Box', 'off', ...
        'FontSize', 6.4, 'NumColumns', 3);
else
    title(ax, 'ST1-68 does not protect during co-infection', ...
        'FontWeight', 'bold');
end
set(ax, 'FontName', 'Arial', 'FontSize', 8.2, 'Box', 'off', ...
    'TickDir', 'out', 'Layer', 'top');
end

function plotFocalCfu(ax, cfu)
strains = ["ST1-75", "ST1-6", "ST1-68", "ST1-49"];
colors = strainColors(strains);
hold(ax, 'on');
fill(ax, [0.45, 4.55, 4.55, 0.45], [0.18, 0.18, 1.18, 1.18], ...
    [0.94, 0.94, 0.94], 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
for i = 1:numel(strains)
    one = cfu(cfu.Candidate == strains(i), :);
    measured = ~ismissing(one.cdiffcfu);
    positive = measured & one.Detected == 1;
    nonDetect = measured & one.Detected == 0;
    logPositive = log10(one.cdiffcfu(positive));
    pale = 0.72 * [1, 1, 1] + 0.28 * colors(i, :);
    quartiles = prctile(logPositive, [25, 50, 75]);
    interquartileRange = quartiles(3) - quartiles(1);
    lowerFence = quartiles(1) - 1.5 * interquartileRange;
    upperFence = quartiles(3) + 1.5 * interquartileRange;
    lowerWhisker = min(logPositive(logPositive >= lowerFence));
    upperWhisker = max(logPositive(logPositive <= upperFence));
    patch(ax, i + [-0.27, 0.27, 0.27, -0.27], ...
        [quartiles(1), quartiles(1), quartiles(3), quartiles(3)], pale, ...
        'EdgeColor', colors(i, :), 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    plot(ax, [i, i], [lowerWhisker, quartiles(1)], '-', ...
        'Color', colors(i, :), 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    plot(ax, [i, i], [quartiles(3), upperWhisker], '-', ...
        'Color', colors(i, :), 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    plot(ax, i + [-0.14, 0.14], [lowerWhisker, lowerWhisker], '-', ...
        'Color', colors(i, :), 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    plot(ax, i + [-0.14, 0.14], [upperWhisker, upperWhisker], '-', ...
        'Color', colors(i, :), 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    plot(ax, i + [-0.27, 0.27], [quartiles(2), quartiles(2)], '-', ...
        'Color', colors(i, :), 'LineWidth', 1.6, ...
        'HandleVisibility', 'off');

    positiveJitter = linspace(-0.22, 0.22, numel(logPositive));
    scatter(ax, i + positiveJitter, logPositive, 16, colors(i, :), ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.35, ...
        'HandleVisibility', 'off');
    if any(nonDetect)
        nonDetectJitter = linspace(-0.14, 0.14, sum(nonDetect));
        scatter(ax, i + nonDetectJitter, ...
            repmat(0.68, sum(nonDetect), 1), 24, colors(i, :), 'x', ...
            'LineWidth', 1.1, 'HandleVisibility', 'off');
    end
    text(ax, i, 6.25, sprintf('%d/%d', sum(positive), sum(measured)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 7.3, 'FontWeight', 'bold', 'Color', colors(i, :));
end
xlim(ax, [0.45, 4.55]);
ylim(ax, [0, 9.25]);
xticks(ax, 1:4);
xticklabels(ax, strains);
yticks(ax, [0.68, 4, 5, 6, 7, 8, 9]);
yticklabels(ax, {'ND', '10^4', '10^5', '10^6', '10^7', '10^8', '10^9'});
ylabel(ax, 'Fecal {\it C. difficile} (CFU/g)');
title(ax, 'Fecal recovery after mono-colonization', ...
    'FontWeight', 'bold');
text(ax, 0.98, 0.96, 'Detected / measured', 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 6.8, 'Color', [0.30, 0.30, 0.30]);
set(ax, 'FontName', 'Arial', 'FontSize', 8.2, 'Box', 'off', ...
    'TickDir', 'out', 'Layer', 'top');
end

function colors = strainColors(strains)
colors = repmat([0.64, 0.66, 0.68], numel(strains), 1);
key = ["ST1-75", "ST1-6", "ST1-68", "ST1-49"];
keyColors = [0, 114, 178; 142, 90, 169; 213, 94, 0; 204, 121, 167] / 255;
for i = 1:numel(key)
    colors(strains == key(i), :) = repmat(keyColors(i, :), ...
        sum(strains == key(i)), 1);
end
end

function panelLabel(ax, label)
text(ax, 0.0, 1.065, label, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 15, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
