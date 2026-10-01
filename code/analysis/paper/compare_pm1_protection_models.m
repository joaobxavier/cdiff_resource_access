%% Compare scalar PM1 protection models
% Report model comparisons, out-of-fold predictions and influence analyses.

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

nPermutation = 5000;
rng(270829, 'twister');

traitFile = fullfile(repoRoot, 'results', 'tables', 'main', ...
    'figure4_phylogeny_aligned_traits_and_breadth.csv');
callFile = fullfile(repoRoot, 'results', 'tables', 'main', ...
    'figure4_pm1_call_matrix_plotted_order.csv');
traits = readtable(traitFile, 'TextType', 'string');
calls = readtable(callFile, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
assert(height(traits) == 21 && height(calls) == 95, ...
    'Expected 21 ST1 strains and 95 PM1 substrates.');

[candidateCalls, vpiCalls] = align_call_matrix(calls, traits.Strain);
assert(sum(candidateCalls, 'all') + sum(vpiCalls) == 700, ...
    'The audited 700-call PM1 matrix changed unexpectedly.');

y = traits.Protection_Effect;
disease = traits.Mono_Disease_Effect;
total = traits.Total_Breadth;
private = traits.ST1_Private;
vpiPrivate = traits.VPI_Private;
shared = traits.Shared_With_VPI;
ratio = private ./ vpiPrivate;
logRatio = log1p(ratio);
privateShare = private ./ (private + vpiPrivate);
persistence = traits.Model_ST1_Persistence_Fraction;
aminoRows = calls.Chemical_Class == "amino acid";
aminoTotal = sum(candidateCalls(:, aminoRows'), 2);
aminoPrivate = sum(candidateCalls(:, aminoRows' & ~vpiCalls), 2);
assert(all(vpiPrivate > 0) && all(isfinite(ratio)), ...
    'The private-breadth ratio is undefined for one or more strains.');

sigmaProtection = recover_protection_covariance(traits.Strain, y);

specifications = build_model_specifications(disease, total, private, ...
    ratio, logRatio, privateShare, aminoTotal, aminoPrivate, ...
    persistence, shared);
[modelResults, predictions] = evaluate_models(specifications, y, ...
    sigmaProtection, nPermutation);
writetable(modelResults, fullfile(tableDir, ...
    'simple_pm1_model_comparison.csv'));
writetable([table(traits.Strain, y, 'VariableNames', ...
    {'Strain', 'Observed_Protection'}), predictions], ...
    fullfile(tableDir, 'simple_pm1_model_out_of_fold_predictions.csv'));

influenceResults = evaluate_influence(specifications, traits.Strain, y);
writetable(influenceResults, fullfile(tableDir, ...
    'simple_pm1_model_focal_sensitivity.csv'));

classResults = evaluate_univariate_classes(candidateCalls, ...
    calls.Chemical_Class, y);
writetable(classResults, fullfile(tableDir, ...
    'simple_pm1_univariate_class_comparison.csv'));

figureFile = fullfile(figureDir, ...
    'pm1_model_comparison.png');
make_model_decision_figure(figureFile, traits, aminoTotal, ...
    modelResults);

summaryFile = fullfile(tableDir, 'simple_pm1_model_summary.txt');
write_summary(summaryFile, modelResults, influenceResults, classResults, ...
    nPermutation);

fprintf('Wrote simple PM1 model comparison outputs to %s\n', reviewDir);

%% Local functions
function [matrix, vpi] = align_call_matrix(calls, strains)
matrix = false(numel(strains), height(calls));
for i = 1:numel(strains)
    variableName = strrep(strains(i), '-', '_');
    assert(ismember(variableName, string(calls.Properties.VariableNames)), ...
        'PM1 matrix is missing %s.', strains(i));
    matrix(i, :) = logical(calls.(variableName))';
end
vpi = logical(calls.VPI10463)';
end

function specifications = build_model_specifications(disease, total, ...
    private, ratio, logRatio, privateShare, aminoTotal, aminoPrivate, ...
    persistence, shared)
Model = ["Intercept only"; "Disease only"; ...
    "1. Total breadth"; "2. Total breadth + disease"; ...
    "3. Candidate-private breadth"; ...
    "Candidate-private breadth + disease"; ...
    "4. Private / VPI-private ratio"; ...
    "Private / VPI-private ratio + disease"; ...
    "5. Amino-acid breadth"; "6. Amino-acid breadth + disease"; ...
    "Amino-acid private breadth"; ...
    "Amino-acid private breadth + disease"; ...
    "log(1 + private/VPI-private)"; ...
    "log(1 + private/VPI-private) + disease"; ...
    "Private share of differentiated calls"; ...
    "Private share + disease"; ...
    "Model persistence fraction"; ...
    "Model persistence fraction + disease"; ...
    "Private/shared partition + disease"];
Predictors = {zeros(numel(disease), 0), disease, total, ...
    [total, disease], private, [private, disease], ratio, ...
    [ratio, disease], aminoTotal, [aminoTotal, disease], ...
    aminoPrivate, [aminoPrivate, disease], logRatio, ...
    [logRatio, disease], privateShare, [privateShare, disease], ...
    persistence, [persistence, disease], [private, shared, disease]};
ResourceColumns = {[], [], 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
    1, 1, 1, 1, [1, 2]};
DiseaseColumn = [NaN; 1; NaN; 2; NaN; 2; NaN; 2; NaN; 2; ...
    NaN; 2; NaN; 2; NaN; 2; NaN; 2; 3];
User_Specified = [false; false; true; true; true; false; true; false; ...
    true; true; false; false; false; false; false; false; false; ...
    false; false];
Narrative_Group = ["reference"; "reference"; "total"; "total"; ...
    "private"; "private"; "ratio"; "ratio"; "amino acid"; ...
    "amino acid"; "amino acid private"; "amino acid private"; ...
    "ratio sensitivity"; "ratio sensitivity"; ...
    "ratio sensitivity"; "ratio sensitivity"; "model derived"; ...
    "model derived"; "partition"];
specifications = table(Model, Predictors', ResourceColumns', ...
    DiseaseColumn, User_Specified, Narrative_Group, ...
    'VariableNames', {'Model', 'Predictors', 'Resource_Columns', ...
    'Disease_Column', 'User_Specified', 'Narrative_Group'});
end

function [results, predictions] = evaluate_models(specifications, y, ...
    sigma, nPermutation)
nModels = height(specifications);
n = numel(y);
N_Predictors = nan(nModels, 1);
OLS_R2 = nan(nModels, 1);
OLS_Adjusted_R2 = nan(nModels, 1);
OLS_AICc = nan(nModels, 1);
GLS_AICc = nan(nModels, 1);
Delta_GLS_AICc_vs_Disease = nan(nModels, 1);
Resource_OLS_P_Value = nan(nModels, 1);
Resource_GLS_LRT_P_Value = nan(nModels, 1);
Resource_Permutation_P_Value = nan(nModels, 1);
LOSO_RMSE = nan(nModels, 1);
LOSO_MAE = nan(nModels, 1);
Cross_Validated_R2 = nan(nModels, 1);
LOSO_Spearman_Rho = nan(nModels, 1);
predictions = table();

diseaseIndex = find(specifications.Model == "Disease only", 1);
diseasePredictor = specifications.Predictors{diseaseIndex};
XBaseDisease = [ones(n, 1), diseasePredictor];
glsDisease = fit_meta_gls(y, XBaseDisease, sigma);

for i = 1:nModels
    predictors = specifications.Predictors{i};
    X = [ones(n, 1), predictors];
    ols = fitlm(predictors, y);
    gls = fit_meta_gls(y, X, sigma);
    prediction = leave_one_out_predict(y, predictors);
    variableName = matlab.lang.makeValidName(specifications.Model(i));
    predictions.(variableName) = prediction;

    N_Predictors(i) = size(predictors, 2);
    OLS_R2(i) = ols.Rsquared.Ordinary;
    OLS_Adjusted_R2(i) = ols.Rsquared.Adjusted;
    OLS_AICc(i) = ordinary_aicc(y, X);
    GLS_AICc(i) = gls.AICc;
    Delta_GLS_AICc_vs_Disease(i) = gls.AICc - glsDisease.AICc;
    LOSO_RMSE(i) = sqrt(mean((y - prediction).^2));
    LOSO_MAE(i) = mean(abs(y - prediction));
    Cross_Validated_R2(i) = 1 - sum((y - prediction).^2) / ...
        sum((y - mean(y)).^2);
    if std(prediction) > 0
        LOSO_Spearman_Rho(i) = corr(y, prediction, 'Type', 'Spearman');
    end

    resourceColumns = specifications.Resource_Columns{i};
    if isempty(resourceColumns)
        continue;
    end
    coefficientColumns = resourceColumns + 1;
    contrast = zeros(numel(resourceColumns), size(X, 2));
    for j = 1:numel(resourceColumns)
        contrast(j, coefficientColumns(j)) = 1;
    end
    Resource_OLS_P_Value(i) = coefTest(ols, contrast);
    diseaseColumn = specifications.Disease_Column(i);
    if isnan(diseaseColumn)
        XReduced = ones(n, 1);
    else
        XReduced = [ones(n, 1), predictors(:, diseaseColumn)];
    end
    reducedGLS = fit_meta_gls(y, XReduced, sigma);
    df = rank(X) - rank(XReduced);
    Resource_GLS_LRT_P_Value(i) = likelihood_ratio_p( ...
        reducedGLS, gls, df);
    Resource_Permutation_P_Value(i) = nested_permutation_test( ...
        y, XReduced, X, nPermutation, 3100 + i);
end

diseaseRMSE = LOSO_RMSE(diseaseIndex);
Delta_LOSO_RMSE_vs_Disease = LOSO_RMSE - diseaseRMSE;
Resource_BH_Q_Value_All_Models = bh_adjust(Resource_Permutation_P_Value);
specified = specifications.User_Specified & ...
    ~isnan(Resource_Permutation_P_Value);
Resource_BH_Q_Value_User_Specified = nan(nModels, 1);
Resource_BH_Q_Value_User_Specified(specified) = ...
    bh_adjust(Resource_Permutation_P_Value(specified));

results = [specifications(:, {'Model', 'User_Specified', ...
    'Narrative_Group'}), table(N_Predictors, OLS_R2, OLS_Adjusted_R2, ...
    OLS_AICc, GLS_AICc, Delta_GLS_AICc_vs_Disease, ...
    Resource_OLS_P_Value, Resource_GLS_LRT_P_Value, ...
    Resource_Permutation_P_Value, Resource_BH_Q_Value_All_Models, ...
    Resource_BH_Q_Value_User_Specified, LOSO_RMSE, ...
    Delta_LOSO_RMSE_vs_Disease, LOSO_MAE, Cross_Validated_R2, ...
    LOSO_Spearman_Rho)];
end

function results = evaluate_influence(specifications, strains, y)
mainModels = ["1. Total breadth"; "2. Total breadth + disease"; ...
    "3. Candidate-private breadth"; ...
    "Candidate-private breadth + disease"; ...
    "4. Private / VPI-private ratio"; ...
    "Private / VPI-private ratio + disease"; ...
    "5. Amino-acid breadth"; "6. Amino-acid breadth + disease"];
keepSets = {true(numel(y), 1), strains ~= "ST1-49", ...
    strains ~= "ST1-68", strains ~= "ST1-49" & strains ~= "ST1-68"};
sensitivityNames = ["All strains"; "Without ST1-49"; ...
    "Without ST1-68"; "Without ST1-49 and ST1-68"];
Model = strings(0, 1);
Sensitivity = strings(0, 1);
N_Strains = zeros(0, 1);
Feature_Estimate = zeros(0, 1);
Feature_P_Value = zeros(0, 1);
OLS_R2 = zeros(0, 1);
Spearman_Rho_Feature_Protection = zeros(0, 1);
for m = 1:numel(mainModels)
    row = find(specifications.Model == mainModels(m), 1);
    predictors = specifications.Predictors{row};
    feature = predictors(:, specifications.Resource_Columns{row}(1));
    for i = 1:numel(keepSets)
        keep = keepSets{i};
        model = fitlm(predictors(keep, :), y(keep));
        coefficientRow = specifications.Resource_Columns{row}(1) + 1;
        Model(end + 1, 1) = mainModels(m); %#ok<AGROW>
        Sensitivity(end + 1, 1) = sensitivityNames(i); %#ok<AGROW>
        N_Strains(end + 1, 1) = sum(keep); %#ok<AGROW>
        Feature_Estimate(end + 1, 1) = ...
            model.Coefficients.Estimate(coefficientRow); %#ok<AGROW>
        Feature_P_Value(end + 1, 1) = ...
            model.Coefficients.pValue(coefficientRow); %#ok<AGROW>
        OLS_R2(end + 1, 1) = model.Rsquared.Ordinary; %#ok<AGROW>
        Spearman_Rho_Feature_Protection(end + 1, 1) = corr( ...
            feature(keep), y(keep), 'Type', 'Spearman'); %#ok<AGROW>
    end
end
results = table(Model, Sensitivity, N_Strains, Feature_Estimate, ...
    Feature_P_Value, OLS_R2, Spearman_Rho_Feature_Protection);
end

function results = evaluate_univariate_classes(matrix, substrateClass, y)
classOrder = ["amino acid"; "carbohydrate"; "carboxylic acid"; ...
    "fatty acid"; "amine"; "alcohol"; "amide"; "ester"];
nClass = numel(classOrder);
Class = title_case(classOrder);
N_Substrates = nan(nClass, 1);
Spearman_Rho = nan(nClass, 1);
Spearman_P_Value = nan(nClass, 1);
OLS_R2 = nan(nClass, 1);
OLS_P_Value = nan(nClass, 1);
LOSO_RMSE = nan(nClass, 1);
for i = 1:nClass
    wells = substrateClass == classOrder(i);
    feature = sum(matrix(:, wells'), 2);
    N_Substrates(i) = sum(wells);
    [Spearman_Rho(i), Spearman_P_Value(i)] = corr( ...
        feature, y, 'Type', 'Spearman');
    model = fitlm(feature, y);
    OLS_R2(i) = model.Rsquared.Ordinary;
    OLS_P_Value(i) = model.Coefficients.pValue(2);
    prediction = leave_one_out_predict(y, feature);
    LOSO_RMSE(i) = sqrt(mean((y - prediction).^2));
end
Spearman_BH_Q_Value = bh_adjust(Spearman_P_Value);
results = table(Class, N_Substrates, Spearman_Rho, Spearman_P_Value, ...
    Spearman_BH_Q_Value, OLS_R2, OLS_P_Value, LOSO_RMSE);
results = sortrows(results, 'Spearman_P_Value', 'ascend');
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
missingWeight = ismissing(data.relweight);
assert(all(data.death(missingWeight) == 1), ...
    'Only recorded terminal-event rows may receive score zeros.');
data.relweight(missingWeight) = 0;
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

function aicc = ordinary_aicc(y, X)
n = numel(y);
residual = y - X * (pinv(X) * y);
rss = sum(residual.^2);
k = size(X, 2) + 1;
aic = n * (log(2 * pi) + 1 + log(rss / n)) + 2 * k;
aicc = aic + 2 * k * (k + 1) / max(n - k - 1, 1);
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

function q = bh_adjust(p)
q = nan(size(p));
valid = ~isnan(p);
values = p(valid);
[sorted, order] = sort(values);
m = numel(sorted);
adjusted = sorted .* m ./ (1:m)';
adjusted = flipud(cummin(flipud(adjusted)));
adjusted = min(adjusted, 1);
restored = nan(m, 1);
restored(order) = adjusted;
q(valid) = restored;
end

function text = title_case(text)
text = string(text);
for i = 1:numel(text)
    if strlength(text(i)) > 0
        text(i) = upper(extractBefore(text(i), 2)) + ...
            extractAfter(text(i), 1);
    end
end
end

function make_model_decision_figure(file, traits, aminoTotal, results)
colors = struct('amino', [19, 126, 125] / 255, ...
    'disease', [75, 85, 99] / 255, 'private', [45, 91, 132] / 255, ...
    'total', [143, 149, 155] / 255, 'focus', [200, 79, 50] / 255, ...
    'grid', [220, 224, 228] / 255);
figureModels = ["Disease only"; "1. Total breadth"; ...
    "2. Total breadth + disease"; "3. Candidate-private breadth"; ...
    "Candidate-private breadth + disease"; ...
    "4. Private / VPI-private ratio"; ...
    "5. Amino-acid breadth"; "6. Amino-acid breadth + disease"];
[matched, rows] = ismember(figureModels, results.Model);
assert(all(matched), 'A comparison model is unavailable for plotting.');
plotData = results(rows, :);
labels = ["Disease"; "Total"; "Total + disease"; "Private"; ...
    "Private + disease"; "Private/VPI-private"; "Amino acid"; ...
    "Amino acid + disease"];
barColors = repmat(colors.total, numel(labels), 1);
barColors(1, :) = colors.disease;
barColors(4:5, :) = repmat(colors.private, 2, 1);
barColors(7:8, :) = repmat(colors.amino, 2, 1);

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [1, 1, 7.2, 7.6]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

axA = nexttile(layout, 1);
yPosition = (1:numel(labels))';
hold(axA, 'on');
for i = 1:numel(labels)
    plot(axA, plotData.Delta_GLS_AICc_vs_Disease(i), yPosition(i), ...
        'o', 'MarkerSize', 7, 'MarkerFaceColor', barColors(i, :), ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.7);
end
xline(axA, 0, '--', 'Disease reference', 'Color', colors.disease, ...
    'LabelVerticalAlignment', 'bottom');
set(axA, 'YTick', yPosition, 'YTickLabel', labels, 'YDir', 'reverse');
xlabel(axA, '\DeltaAICc versus disease only');
title(axA, 'AICc favors amino acid + disease');
panel_label(axA, 'A');

axB = nexttile(layout, 2);
hold(axB, 'on');
for i = 1:numel(labels)
    plot(axB, plotData.LOSO_RMSE(i), yPosition(i), 'o', ...
        'MarkerSize', 7, 'MarkerFaceColor', barColors(i, :), ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.7);
end
diseaseRMSE = plotData.LOSO_RMSE(1);
xline(axB, diseaseRMSE, '--', 'Disease reference', ...
    'Color', colors.disease, 'LabelVerticalAlignment', 'bottom');
set(axB, 'YTick', yPosition, 'YTickLabel', [], 'YDir', 'reverse');
xlabel(axB, 'Leave-one-strain-out RMSE');
title(axB, 'Amino acid is best resource-only');
panel_label(axB, 'B');

axC = nexttile(layout, 3);
focus = ismember(traits.Strain, ["ST1-75", "ST1-6", "ST1-68"]);
scatter(axC, aminoTotal(~focus), traits.Protection_Effect(~focus), ...
    34, colors.total, 'filled', 'MarkerFaceAlpha', 0.74);
hold(axC, 'on');
scatter(axC, aminoTotal(focus), traits.Protection_Effect(focus), ...
    52, colors.amino, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.7);
add_fit_line(axC, aminoTotal, traits.Protection_Effect, colors.amino);
label_focus(axC, aminoTotal, traits.Protection_Effect, ...
    traits.Strain, focus);
[rho, pValue] = corr(aminoTotal, traits.Protection_Effect, ...
    'Type', 'Spearman');
text(axC, 0.04, 0.94, sprintf('Spearman \\rho = %.2f\np = %.4f', ...
    rho, pValue), 'Units', 'normalized', 'VerticalAlignment', 'top');
xlabel(axC, 'Amino-acid positive-call breadth');
ylabel(axC, 'Protection score');
title(axC, 'Strong standalone association');
panel_label(axC, 'C');

axD = nexttile(layout, 4);
[xResidual, yResidual] = added_variable(aminoTotal, ...
    traits.Protection_Effect, traits.Mono_Disease_Effect);
scatter(axD, xResidual(~focus), yResidual(~focus), 34, colors.total, ...
    'filled', 'MarkerFaceAlpha', 0.74);
hold(axD, 'on');
scatter(axD, xResidual(focus), yResidual(focus), 52, colors.focus, ...
    'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.7);
add_fit_line(axD, xResidual, yResidual, colors.focus);
label_focus(axD, xResidual, yResidual, traits.Strain, focus);
xline(axD, 0, ':', 'Color', colors.grid);
yline(axD, 0, ':', 'Color', colors.grid);
adjustedRow = results(results.Model == ...
    "6. Amino-acid breadth + disease", :);
text(axD, 0.04, 0.94, sprintf(['incremental GLS p = %.4f\n' ...
    'LOSO RMSE = %.2f'], adjustedRow.Resource_GLS_LRT_P_Value, ...
    adjustedRow.LOSO_RMSE), 'Units', 'normalized', ...
    'VerticalAlignment', 'top');
xlabel(axD, 'Amino-acid breadth | disease');
ylabel(axD, 'Protection | disease');
title(axD, 'Adjustment is ST1-68-sensitive');
panel_label(axD, 'D');

allAxes = [axA, axB, axC, axD];
set(allAxes, 'FontName', 'Arial', 'FontSize', 8.5, 'LineWidth', 0.8, ...
    'Box', 'off', 'TickDir', 'out');
grid(allAxes, 'on');
for ax = allAxes
    ax.GridColor = colors.grid;
    ax.GridAlpha = 0.45;
end
sgtitle(fig, 'Simple PM1 model comparison for the manuscript', ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
exportgraphics(fig, file, 'Resolution', 300);
close(fig);
end

function [xResidual, yResidual] = added_variable(x, y, covariate)
design = [ones(numel(y), 1), covariate];
xResidual = x - design * (pinv(design) * x);
yResidual = y - design * (pinv(design) * y);
end

function add_fit_line(ax, x, y, color)
limits = [min(x), max(x)];
beta = pinv([ones(numel(x), 1), x]) * y;
plot(ax, limits, [ones(2, 1), limits'] * beta, ...
    'Color', color, 'LineWidth', 1.7);
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
text(ax, -0.13, 1.02, label, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

function write_summary(file, results, influence, classes, nPermutation)
fid = fopen(file, 'w');
assert(fid > 0, 'Unable to open model summary output.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Scalar PM1 protection model comparison\n');
fprintf(fid, 'Residual permutations per resource test: %d\n\n', ...
    nPermutation);
userRows = results.User_Specified | results.Model == "Disease only" | ...
    results.Model == "Candidate-private breadth + disease" | ...
    results.Model == "Private / VPI-private ratio + disease";
selected = results(userRows, :);
for i = 1:height(selected)
    fprintf(fid, ['%s: GLS AICc %.3f; delta vs disease %.3f; ' ...
        'resource GLS p %.4g; permutation p %.4g; LOSO RMSE %.3f; ' ...
        'CV R2 %.3f.\n'], selected.Model(i), selected.GLS_AICc(i), ...
        selected.Delta_GLS_AICc_vs_Disease(i), ...
        selected.Resource_GLS_LRT_P_Value(i), ...
        selected.Resource_Permutation_P_Value(i), ...
        selected.LOSO_RMSE(i), selected.Cross_Validated_R2(i));
end
fprintf(fid, '\nAmino-acid class comparison:\n');
amino = classes(classes.Class == "Amino acid", :);
fprintf(fid, ['Spearman rho %.3f; p %.4g; BH q %.4g; ' ...
    'LOSO RMSE %.3f.\n'], amino.Spearman_Rho, amino.Spearman_P_Value, ...
    amino.Spearman_BH_Q_Value, amino.LOSO_RMSE);
fprintf(fid, '\nAmino-acid focal sensitivity:\n');
rows = contains(influence.Model, "Amino-acid");
subset = influence(rows, :);
for i = 1:height(subset)
    fprintf(fid, '%s | %s: beta %.3f; p %.4g; R2 %.3f; rho %.3f.\n', ...
        subset.Model(i), subset.Sensitivity(i), ...
        subset.Feature_Estimate(i), subset.Feature_P_Value(i), ...
        subset.OLS_R2(i), subset.Spearman_Rho_Feature_Protection(i));
end
end
