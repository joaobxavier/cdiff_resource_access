%% Test whether the PM1 amino-acid breadth association is distributed and specific
% Perform two substrate-class sensitivity analyses:
% (1) remove each amino-acid well from the class score one at a time; and
% (2) compare amino-acid breadth with breadth across all remaining PM1 wells.
% No substrate subset is selected from protection outcomes.

clearvars;
close all;

scriptDir = fileparts(mfilename('fullpath'));
localRoot = fileparts(fileparts(scriptDir));
repoRoot = fileparts(localRoot);
addpath(genpath(fullfile(localRoot, 'analysis')));

reviewDir = fullfile(repoRoot, 'results', 'analyses', ...
    'pm1_models');
figureDir = fullfile(reviewDir, 'figures');
tableDir = fullfile(reviewDir, 'tables');
if ~isfolder(figureDir); mkdir(figureDir); end
if ~isfolder(tableDir); mkdir(tableDir); end
set(0, 'DefaultFigureVisible', 'off');

nPermutation = 10000;

traitFile = fullfile(repoRoot, 'results', 'tables', 'main', ...
    'figure4_phylogeny_aligned_traits_and_breadth.csv');
callFile = fullfile(repoRoot, 'results', 'tables', 'main', ...
    'figure4_pm1_call_matrix_plotted_order.csv');
traits = readtable(traitFile, 'TextType', 'string');
calls = readtable(callFile, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
assert(height(traits) == 21 && height(calls) == 95, ...
    'Expected 21 ST1 strains and 95 PM1 substrates.');

candidateCalls = align_call_matrix(calls, traits.Strain);
vpiCalls = logical(calls.VPI10463);
assert(sum(candidateCalls, 'all') + sum(vpiCalls) == 700, ...
    'The audited 700-call PM1 matrix changed unexpectedly.');

aminoRows = calls.Chemical_Class == "amino acid";
assert(sum(aminoRows) == 16, 'Expected 16 amino-acid-class PM1 wells.');
y = traits.Protection_Effect;
disease = traits.Mono_Disease_Effect;
aminoBreadth = sum(candidateCalls(:, aminoRows'), 2);
nonAminoBreadth = sum(candidateCalls(:, ~aminoRows'), 2);
assert(all(aminoBreadth + nonAminoBreadth == traits.Total_Breadth), ...
    'Amino-acid and non-amino-acid breadth do not reconstruct total breadth.');

sigmaProtection = recover_protection_covariance(traits.Strain, y);

scoreResults = evaluate_scalar_scores(y, disease, aminoBreadth, ...
    nonAminoBreadth, sigmaProtection, nPermutation);
writetable(scoreResults, fullfile(tableDir, ...
    'amino_acid_specificity_score_comparison.csv'));

jointResults = evaluate_joint_model(y, disease, aminoBreadth, ...
    nonAminoBreadth, sigmaProtection, nPermutation);
writetable(jointResults, fullfile(tableDir, ...
    'amino_acid_specificity_joint_model.csv'));

deletionResults = evaluate_amino_deletions(candidateCalls(:, aminoRows'), ...
    calls.Substrate(aminoRows), y, disease);
writetable(deletionResults, fullfile(tableDir, ...
    'amino_acid_leave_one_well_out.csv'));

predictionResults = make_prediction_table(traits.Strain, y, disease, ...
    aminoBreadth, nonAminoBreadth);
writetable(predictionResults, fullfile(tableDir, ...
    'amino_acid_specificity_out_of_fold_predictions.csv'));

figureFile = fullfile(figureDir, ...
    'amino_acid_specificity.png');
make_review_figure(figureFile, traits, aminoBreadth, nonAminoBreadth, ...
    deletionResults, scoreResults, predictionResults);

summaryFile = fullfile(tableDir, ...
    'amino_acid_specificity_summary.txt');
write_summary(summaryFile, scoreResults, jointResults, deletionResults, ...
    nPermutation);

fprintf('Wrote amino-acid specificity outputs to %s\n', reviewDir);

%% Local functions
function matrix = align_call_matrix(calls, strains)
matrix = false(numel(strains), height(calls));
for i = 1:numel(strains)
    variableName = strrep(strains(i), '-', '_');
    assert(ismember(variableName, string(calls.Properties.VariableNames)), ...
        'PM1 matrix is missing %s.', strains(i));
    matrix(i, :) = logical(calls.(variableName))';
end
end

function results = evaluate_scalar_scores(y, disease, amino, nonAmino, ...
    sigma, nPermutation)
Score = ["Amino-acid breadth"; "Non-amino-acid breadth"];
features = {amino, nonAmino};
nScore = numel(Score);
N_Wells = [16; 79];
Spearman_Rho = nan(nScore, 1);
Spearman_P_Value = nan(nScore, 1);
Partial_Rank_Rho = nan(nScore, 1);
Partial_Rank_P_Value = nan(nScore, 1);
Adjusted_OLS_Estimate = nan(nScore, 1);
Adjusted_OLS_P_Value = nan(nScore, 1);
Adjusted_GLS_Estimate = nan(nScore, 1);
Adjusted_GLS_LRT_P_Value = nan(nScore, 1);
Adjusted_Permutation_P_Value = nan(nScore, 1);
Unadjusted_LOSO_RMSE = nan(nScore, 1);
Adjusted_LOSO_RMSE = nan(nScore, 1);
Adjusted_CV_R2 = nan(nScore, 1);
n = numel(y);
diseaseGLS = fit_meta_gls(y, [ones(n, 1), disease], sigma);
for i = 1:nScore
    x = features{i};
    [Spearman_Rho(i), Spearman_P_Value(i)] = corr( ...
        x, y, 'Type', 'Spearman');
    [Partial_Rank_Rho(i), Partial_Rank_P_Value(i)] = ...
        partial_rank(x, y, disease);
    adjustedModel = fitlm([disease, x], y);
    Adjusted_OLS_Estimate(i) = adjustedModel.Coefficients.Estimate(3);
    Adjusted_OLS_P_Value(i) = adjustedModel.Coefficients.pValue(3);
    adjustedGLS = fit_meta_gls(y, [ones(n, 1), disease, x], sigma);
    Adjusted_GLS_Estimate(i) = adjustedGLS.Beta(3);
    Adjusted_GLS_LRT_P_Value(i) = likelihood_ratio_p( ...
        diseaseGLS, adjustedGLS, 1);
    Adjusted_Permutation_P_Value(i) = nested_permutation_test(y, ...
        [ones(n, 1), disease], [ones(n, 1), disease, x], ...
        nPermutation, 5140 + i);
    prediction = leave_one_out_predict(y, x);
    Unadjusted_LOSO_RMSE(i) = sqrt(mean((y - prediction).^2));
    adjustedPrediction = leave_one_out_predict(y, [disease, x]);
    Adjusted_LOSO_RMSE(i) = sqrt(mean((y - adjustedPrediction).^2));
    Adjusted_CV_R2(i) = cv_r2(y, adjustedPrediction);
end
results = table(Score, N_Wells, Spearman_Rho, Spearman_P_Value, ...
    Partial_Rank_Rho, Partial_Rank_P_Value, Adjusted_OLS_Estimate, ...
    Adjusted_OLS_P_Value, Adjusted_GLS_Estimate, ...
    Adjusted_GLS_LRT_P_Value, Adjusted_Permutation_P_Value, ...
    Unadjusted_LOSO_RMSE, Adjusted_LOSO_RMSE, Adjusted_CV_R2);
end

function results = evaluate_joint_model(y, disease, amino, nonAmino, ...
    sigma, nPermutation)
n = numel(y);
XFull = [ones(n, 1), disease, amino, nonAmino];
XNoAmino = [ones(n, 1), disease, nonAmino];
XNoNonAmino = [ones(n, 1), disease, amino];
fullOLS = fitlm([disease, amino, nonAmino], y);
fullGLS = fit_meta_gls(y, XFull, sigma);
noAminoGLS = fit_meta_gls(y, XNoAmino, sigma);
noNonAminoGLS = fit_meta_gls(y, XNoNonAmino, sigma);
Term = ["Amino-acid breadth | disease + non-amino breadth"; ...
    "Non-amino breadth | disease + amino-acid breadth"];
OLS_Estimate = [fullOLS.Coefficients.Estimate(3); ...
    fullOLS.Coefficients.Estimate(4)];
OLS_P_Value = [fullOLS.Coefficients.pValue(3); ...
    fullOLS.Coefficients.pValue(4)];
GLS_Estimate = [fullGLS.Beta(3); fullGLS.Beta(4)];
GLS_LRT_P_Value = [likelihood_ratio_p(noAminoGLS, fullGLS, 1); ...
    likelihood_ratio_p(noNonAminoGLS, fullGLS, 1)];
Permutation_P_Value = [nested_permutation_test(y, XNoAmino, XFull, ...
    nPermutation, 5271); nested_permutation_test(y, XNoNonAmino, ...
    XFull, nPermutation, 5272)];
standardizedModel = fitlm([zscore(disease), zscore(amino), ...
    zscore(nonAmino)], zscore(y));
Standardized_OLS_Estimate = [standardizedModel.Coefficients.Estimate(3); ...
    standardizedModel.Coefficients.Estimate(4)];
results = table(Term, OLS_Estimate, OLS_P_Value, GLS_Estimate, ...
    GLS_LRT_P_Value, Permutation_P_Value, Standardized_OLS_Estimate);
end

function results = evaluate_amino_deletions(aminoCalls, names, y, disease)
nWells = numel(names);
Removed_Well = names;
N_Positive_Calls = sum(aminoCalls, 1)';
Spearman_Rho = nan(nWells, 1);
Spearman_P_Value = nan(nWells, 1);
Partial_Rank_Rho = nan(nWells, 1);
Partial_Rank_P_Value = nan(nWells, 1);
Adjusted_OLS_Estimate = nan(nWells, 1);
Adjusted_OLS_P_Value = nan(nWells, 1);
Unadjusted_LOSO_RMSE = nan(nWells, 1);
Adjusted_LOSO_RMSE = nan(nWells, 1);
fullBreadth = sum(aminoCalls, 2);
for i = 1:nWells
    score = fullBreadth - aminoCalls(:, i);
    [Spearman_Rho(i), Spearman_P_Value(i)] = corr( ...
        score, y, 'Type', 'Spearman');
    [Partial_Rank_Rho(i), Partial_Rank_P_Value(i)] = ...
        partial_rank(score, y, disease);
    adjustedModel = fitlm([disease, score], y);
    Adjusted_OLS_Estimate(i) = adjustedModel.Coefficients.Estimate(3);
    Adjusted_OLS_P_Value(i) = adjustedModel.Coefficients.pValue(3);
    prediction = leave_one_out_predict(y, score);
    Unadjusted_LOSO_RMSE(i) = sqrt(mean((y - prediction).^2));
    adjustedPrediction = leave_one_out_predict(y, [disease, score]);
    Adjusted_LOSO_RMSE(i) = sqrt(mean((y - adjustedPrediction).^2));
end
results = table(Removed_Well, N_Positive_Calls, Spearman_Rho, ...
    Spearman_P_Value, Partial_Rank_Rho, Partial_Rank_P_Value, ...
    Adjusted_OLS_Estimate, Adjusted_OLS_P_Value, ...
    Unadjusted_LOSO_RMSE, Adjusted_LOSO_RMSE);
end

function results = make_prediction_table(strain, y, disease, amino, nonAmino)
Observed_Protection = y;
Disease_Only = leave_one_out_predict(y, disease);
Disease_Plus_Amino = leave_one_out_predict(y, [disease, amino]);
Disease_Plus_Non_Amino = leave_one_out_predict(y, [disease, nonAmino]);
Disease_Plus_Both = leave_one_out_predict(y, ...
    [disease, amino, nonAmino]);
results = table(strain, Observed_Protection, Disease_Only, ...
    Disease_Plus_Amino, Disease_Plus_Non_Amino, Disease_Plus_Both, ...
    'VariableNames', {'Strain', 'Observed_Protection', 'Disease_Only', ...
    'Disease_Plus_Amino', 'Disease_Plus_Non_Amino', ...
    'Disease_Plus_Both'});
end

function make_review_figure(file, traits, amino, nonAmino, deletion, ...
    scoreResults, predictions)
colors = struct('amino', [14, 128, 123] / 255, ...
    'other', [130, 142, 151] / 255, 'adjusted', [204, 82, 52] / 255, ...
    'focus', [26, 78, 115] / 255, 'grid', [221, 225, 228] / 255);
fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [1, 1, 7.6, 8.2]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
y = traits.Protection_Effect;
focus = ismember(traits.Strain, ["ST1-75", "ST1-6", "ST1-68"]);

axA = nexttile(layout, 1);
draw_scatter(axA, amino, y, focus, traits.Strain, colors.amino, ...
    scoreResults.Spearman_Rho(1));
xlabel(axA, 'Amino-acid positive-call breadth');
ylabel(axA, 'Protection score');
title(axA, 'Amino-acid breadth captures a strong class signal');
panel_label(axA, 'A');

axB = nexttile(layout, 2);
draw_scatter(axB, nonAmino, y, focus, traits.Strain, colors.other, ...
    scoreResults.Spearman_Rho(2));
xlabel(axB, 'Non-amino-acid positive-call breadth');
ylabel(axB, 'Protection score');
title(axB, 'The remaining 79 wells are less informative');
panel_label(axB, 'B');

axC = nexttile(layout, 3);
deletion = sortrows(deletion, 'Spearman_Rho', 'ascend');
positions = (1:height(deletion))';
plot(axC, deletion.Spearman_Rho, positions, 'o', ...
    'MarkerFaceColor', colors.amino, 'MarkerEdgeColor', 'w', ...
    'MarkerSize', 6.5, 'LineStyle', 'none');
hold(axC, 'on');
plot(axC, deletion.Partial_Rank_Rho, positions, 'o', ...
    'MarkerFaceColor', colors.adjusted, 'MarkerEdgeColor', 'w', ...
    'MarkerSize', 6.5, 'LineStyle', 'none');
set(axC, 'YTick', positions, 'YTickLabel', deletion.Removed_Well, ...
    'TickLabelInterpreter', 'none');
ylim(axC, [0.5, height(deletion) + 0.5]);
xline(axC, scoreResults.Spearman_Rho(1), '-', ...
    'Color', colors.amino, 'Alpha', 0.35);
xline(axC, scoreResults.Partial_Rank_Rho(1), '-', ...
    'Color', colors.adjusted, 'Alpha', 0.35);
xlabel(axC, 'Correlation after removing one amino-acid well');
title(axC, 'The signal survives every single-well deletion');
legend(axC, {'Raw Spearman', 'Disease-adjusted partial rank'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Box', 'off', 'FontSize', 7.2);
panel_label(axC, 'C');

axD = nexttile(layout, 4);
models = ["Disease only"; "+ amino-acid breadth"; ...
    "+ non-amino breadth"; "+ both breadth terms"];
predictionMatrix = [predictions.Disease_Only, ...
    predictions.Disease_Plus_Amino, ...
    predictions.Disease_Plus_Non_Amino, predictions.Disease_Plus_Both];
rmse = sqrt(mean((y - predictionMatrix).^2, 1))';
positions = (1:numel(models))';
pointColors = [colors.other; colors.amino; colors.other; colors.focus];
hold(axD, 'on');
for i = 1:numel(models)
    plot(axD, rmse(i), positions(i), 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', pointColors(i, :), 'MarkerEdgeColor', 'w');
end
set(axD, 'YTick', positions, 'YTickLabel', models, 'YDir', 'reverse');
ylim(axD, [0.5, numel(models) + 0.5]);
xline(axD, rmse(1), '--', 'Color', colors.grid);
xlabel(axD, {'Held-out-strain RMSE', ...
    'lower is better; dashed = disease'});
title(axD, 'Non-amino breadth adds no held-out gain');
panel_label(axD, 'D');

allAxes = [axA, axB, axC, axD];
set(allAxes, 'FontName', 'Arial', 'FontSize', 8.2, 'LineWidth', 0.8, ...
    'Box', 'off', 'TickDir', 'out');
grid(allAxes, 'on');
for ax = allAxes
    ax.GridColor = colors.grid;
    ax.GridAlpha = 0.45;
end
sgtitle(fig, 'Amino-acid breadth is distributed and distinct from residual breadth', ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
exportgraphics(fig, file, 'Resolution', 300);
close(fig);
end

function draw_scatter(ax, x, y, focus, names, color, rho)
scatter(ax, x(~focus), y(~focus), 34, [150, 158, 165] / 255, ...
    'filled', 'MarkerFaceAlpha', 0.75);
hold(ax, 'on');
scatter(ax, x(focus), y(focus), 52, color, 'filled', ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.7);
limits = [min(x), max(x)];
beta = pinv([ones(numel(x), 1), x]) * y;
plot(ax, limits, [ones(2, 1), limits'] * beta, ...
    'Color', color, 'LineWidth', 1.7);
label_focus(ax, x, y, names, focus);
text(ax, 0.04, 0.94, sprintf('Spearman \\rho = %.2f', rho), ...
    'Units', 'normalized', 'VerticalAlignment', 'top');
end

function label_focus(ax, x, y, names, focus)
rows = find(focus);
xLimits = [min(x), max(x)];
for i = 1:numel(rows)
    row = rows(i);
    if x(row) > xLimits(1) + 0.75 * diff(xLimits)
        label = names(row) + "  ";
        alignment = 'right';
    else
        label = "  " + names(row);
        alignment = 'left';
    end
    text(ax, x(row), y(row), label, 'FontName', 'Arial', ...
        'FontSize', 7.5, 'HorizontalAlignment', alignment, ...
        'VerticalAlignment', 'middle');
end
end

function panel_label(ax, label)
text(ax, -0.16, 1.02, label, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

function [rho, pValue] = partial_rank(x, y, disease)
n = numel(y);
design = [ones(n, 1), tiedrank(disease)];
residualizer = eye(n) - design * pinv(design);
xResidual = residualizer * tiedrank(x);
yResidual = residualizer * tiedrank(y);
rho = corr(xResidual, yResidual);
% One fitted disease covariate: partial-correlation t test has n - 3 df.
df = n - 3;
pValue = 2 * tcdf(-abs(rho) * sqrt(df / (1 - rho^2)), df);
end

function prediction = leave_one_out_predict(y, predictors)
n = numel(y);
prediction = nan(n, 1);
for i = 1:n
    train = true(n, 1);
    train(i) = false;
    XTrain = [ones(sum(train), 1), predictors(train, :)];
    beta = pinv(XTrain) * y(train);
    prediction(i) = [1, predictors(i, :)] * beta;
end
end

function value = cv_r2(y, prediction)
value = 1 - sum((y - prediction).^2) / sum((y - mean(y)).^2);
end

function sigma = recover_protection_covariance(strainOrder, expectedEffects)
protectionFile = project_data_file('processed', 'mouse', 'scores', ...
    'ProtectionScreen_CDI_mouse.csv');
raw = readtable(protectionFile, 'TextType', 'string');
raw.experiment = string(raw.experiment);
raw.tx = string(raw.tx);
raw.cdiffstrain = string(raw.cdiffstrain);
raw.mouse = string(raw.mouse);
screen = raw(raw.experiment ~= "ks65", :);
secondary = screen(endsWith(screen.cdiffstrain, ".vpi") | ...
    screen.cdiffstrain == "vpi", :);
[effects, model] = fit_score_model(secondary, "vpi");
[matched, rows] = ismember(strainOrder, effects.Strain);
assert(all(matched), 'Unable to recover every protection coefficient.');
assert(max(abs(effects.Effect(rows) - expectedEffects)) < 1e-8, ...
    'Refitted protection estimates do not match permanent scores.');
indices = effects.Coefficient_Index(rows);
sigma = model.CoefficientCovariance(indices, indices);
sigma = (sigma + sigma') / 2;
end

function [effects, model] = fit_score_model(raw, reference)
data = raw;
data.relweight(ismissing(data.relweight)) = 0;
data.exp_id = categorical(data.experiment + "_" + ...
    data.cdiffstrain + "_" + data.mouse);
data.cdiffstrain = categorical(data.cdiffstrain);
data.cdiffstrain = reordercats(data.cdiffstrain, ...
    [reference; setdiff(categories(data.cdiffstrain), reference)]);
model = fitlme(data, ...
    'relweight ~ cdiffstrain + (1|day) + (1|exp_id)');
coef = model.Coefficients;
names = string(coef.Name);
keep = startsWith(names, "cdiffstrain_st1.");
condition = erase(names(keep), "cdiffstrain_");
condition = erase(condition, ".vpi");
Strain = upper(replace(condition, ".", "-"));
Effect = coef.Estimate(keep);
Coefficient_Index = find(keep);
effects = table(Strain, Effect, Coefficient_Index);
end

function model = fit_meta_gls(y, X, sigma)
n = numel(y);
p = size(X, 2);
upper = max(var(y) * 20, 1);
objective = @(tau2) gls_objective(y, X, sigma, tau2);
[tau2, nll] = fminbnd(objective, 0, upper, ...
    optimset('Display', 'off', 'TolX', 1e-10));
V = sigma + tau2 * eye(n) + 1e-10 * eye(n);
W = pinv(V);
beta = pinv(X' * W * X) * (X' * W * y);
covBeta = pinv(X' * W * X);
se = sqrt(max(diag(covBeta), 0));
tStat = beta ./ se;
pValue = 2 * (1 - tcdf(abs(tStat), max(n - p, 1)));
k = p + 1;
aic = 2 * nll + 2 * k;
aicc = aic + 2 * k * (k + 1) / max(n - k - 1, 1);
model = struct('Beta', beta, 'SE', se, 'PValue', pValue, ...
    'CovBeta', covBeta, 'Tau2', tau2, 'LogLikelihood', -nll, ...
    'AICc', aicc);
end

function nll = gls_objective(y, X, sigma, tau2)
n = numel(y);
V = sigma + tau2 * eye(n) + 1e-10 * eye(n);
[R, flag] = chol(V);
if flag ~= 0
    nll = realmax / 100;
    return;
end
W = pinv(V);
beta = pinv(X' * W * X) * (X' * W * y);
residual = y - X * beta;
logdet = 2 * sum(log(diag(R)));
nll = 0.5 * (n * log(2 * pi) + logdet + residual' * W * residual);
end

function p = likelihood_ratio_p(reduced, full, df)
statistic = max(0, 2 * (full.LogLikelihood - reduced.LogLikelihood));
p = 1 - chi2cdf(statistic, df);
end

function pValue = nested_permutation_test(y, XReduced, XFull, ...
    nPermutation, seed)
reducedFitted = XReduced * (pinv(XReduced) * y);
reducedResidual = y - reducedFitted;
observed = nested_rss_gain(y, XReduced, XFull);
exceed = 0;
rng(seed, 'twister');
for i = 1:nPermutation
    permuted = reducedFitted + reducedResidual(randperm(numel(y)));
    statistic = nested_rss_gain(permuted, XReduced, XFull);
    exceed = exceed + (statistic >= observed - 1e-12);
end
pValue = (exceed + 1) / (nPermutation + 1);
end

function value = nested_rss_gain(y, XReduced, XFull)
reducedResidual = y - XReduced * (pinv(XReduced) * y);
fullResidual = y - XFull * (pinv(XFull) * y);
value = sum(reducedResidual.^2) - sum(fullResidual.^2);
end

function write_summary(file, score, joint, deletion, nPermutation)
fid = fopen(file, 'w');
assert(fid > 0, 'Unable to open amino-acid specificity summary.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Amino-acid breadth specificity analysis\n');
fprintf(fid, 'Prespecified permutations: %d.\n\n', nPermutation);
for i = 1:height(score)
    fprintf(fid, ['%s (%d wells): raw rho %.3f (p %.4g); ' ...
        'partial rho %.3f (p %.4g); adjusted permutation p %.4g; ' ...
        'adjusted LOSO RMSE %.3f; CV R2 %.3f.\n'], ...
        score.Score(i), score.N_Wells(i), score.Spearman_Rho(i), ...
        score.Spearman_P_Value(i), score.Partial_Rank_Rho(i), ...
        score.Partial_Rank_P_Value(i), ...
        score.Adjusted_Permutation_P_Value(i), ...
        score.Adjusted_LOSO_RMSE(i), score.Adjusted_CV_R2(i));
end
fprintf(fid, '\nJoint disease + amino + non-amino model:\n');
for i = 1:height(joint)
    fprintf(fid, ['%s: OLS estimate %.3f (p %.4g); GLS estimate %.3f; ' ...
        'GLS LRT p %.4g; permutation p %.4g; standardized beta %.3f.\n'], ...
        joint.Term(i), joint.OLS_Estimate(i), joint.OLS_P_Value(i), ...
        joint.GLS_Estimate(i), joint.GLS_LRT_P_Value(i), ...
        joint.Permutation_P_Value(i), ...
        joint.Standardized_OLS_Estimate(i));
end
fprintf(fid, ['\nLeave-one-amino-acid-well-out raw rho range %.3f to ' ...
    '%.3f; partial-rank range %.3f to %.3f; adjusted LOSO RMSE ' ...
    'range %.3f to %.3f.\n'], min(deletion.Spearman_Rho), ...
    max(deletion.Spearman_Rho), min(deletion.Partial_Rank_Rho), ...
    max(deletion.Partial_Rank_Rho), min(deletion.Adjusted_LOSO_RMSE), ...
    max(deletion.Adjusted_LOSO_RMSE));
end
