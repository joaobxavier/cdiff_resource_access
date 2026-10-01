%% Matched-batch intracellular GC-MS exploration
% Blank-corrects each bacterial sample against Media_1, Media_2, or Media_3
% from the same confirmed experimental batch. The analysis asks whether the
% five strains occupy distinct intracellular metabolic states. It does not
% interpret pool sizes as fluxes or mechanisms.

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

rng(20260821, 'twister');
nPermGlobal = 9999;
nPermLeaveOneOut = 1999;
nPermPls = 999;
strainOrder = ["ST1-75", "ST1-68", "ST1-6", "VPI10463", "R20291"];
strainColors = [
    36, 104, 162;
    60, 157, 86;
    142, 90, 169;
    217, 71, 63;
    228, 154, 47] ./ 255;

%% Load and verify the balanced matched-batch design
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
assert(numel(unique(sampleID)) == numel(sampleID), ...
    'Sample identifiers must be unique.');
assert(all(isfinite(rawSignal), 'all'), 'All source values must be finite.');
assert(all(rawSignal >= 0, 'all'), 'Source values must be nonnegative.');

bacterialRow = group ~= "Media";
bacterialIndex = find(bacterialRow);
bacterialSample = sampleID(bacterialRow);
bacterialGroup = group(bacterialRow);
bacterialBatch = batch(bacterialRow);
blankCorrected = nan(sum(bacterialRow), numel(metabolite));
matchedMediaSample = strings(sum(bacterialRow), 1);
for i = 1:numel(bacterialIndex)
    sourceRow = bacterialIndex(i);
    mediaRow = find(group == "Media" & batch == batch(sourceRow), 1);
    assert(~isempty(mediaRow), 'Missing matched media for batch %d.', batch(sourceRow));
    blankCorrected(i, :) = rawSignal(sourceRow, :) - rawSignal(mediaRow, :);
    matchedMediaSample(i) = sampleID(mediaRow);
end

%% Define the analyte set before multivariate analysis
% Retain metabolites above their matched blank in at least 12 of 15 bacterial
% samples. This excludes consistently below-blank features while allowing up
% to three nonpositive observations to be left-censored in log/CLR analyses.
nPositive = sum(blankCorrected > 0, 1);
include = nPositive >= 12;
includedMetabolite = metabolite(include);
Xblank = blankCorrected(:, include);
assert(numel(includedMetabolite) == 13, ...
    'The >=12/15 diagnostic rule should retain 13 metabolites.');

minimumPositive = nan(1, sum(include));
for j = 1:sum(include)
    positive = Xblank(Xblank(:, j) > 0, j);
    minimumPositive(j) = min(positive);
end

Xhalf = Xblank;
Xtenth = Xblank;
for j = 1:size(Xblank, 2)
    Xhalf(Xhalf(:, j) <= 0, j) = 0.5 * minimumPositive(j);
    Xtenth(Xtenth(:, j) <= 0, j) = 0.1 * minimumPositive(j);
end
XlogHalf = log2(Xhalf);
XlogTenth = log2(Xtenth);
XclrHalf = XlogHalf - mean(XlogHalf, 2);
XclrTenth = XlogTenth - mean(XlogTenth, 2);
XblankZ = standardizeColumns(Xblank);

%% Export source-level and preprocessing audits
audit = table(metabolite', sum(blankCorrected > 0, 1)', ...
    sum(blankCorrected <= 0, 1)', include', ...
    'VariableNames', {'Metabolite', 'NAboveMatchedBlank', ...
    'NAtOrBelowMatchedBlank', 'IncludedInMultivariateAnalysis'});
writetable(audit, fullfile(tableDir, 'metabolite_inclusion_audit.csv'));

correctedTable = table(bacterialSample, bacterialGroup, bacterialBatch, ...
    matchedMediaSample, 'VariableNames', {'Sample', 'Strain', 'Batch', ...
    'MatchedMediaSample'});
correctedValues = array2table(blankCorrected, ...
    'VariableNames', matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(metabolite)));
correctedTable = [correctedTable, correctedValues];
writetable(correctedTable, fullfile(tableDir, 'blank_corrected_intracellular_values.csv'));

%% PCA on absolute blank-corrected pools and relative composition
[coeffRaw, scoreRaw, ~, ~, explainedRaw] = pca(XblankZ);
[coeffClr, scoreClr, ~, ~, explainedClr] = pca(XclrHalf);
[coeffClrTenth, scoreClrTenth, ~, ~, explainedClrTenth] = pca(XclrTenth);

pcaScores = table(bacterialSample, bacterialGroup, bacterialBatch, ...
    scoreRaw(:, 1), scoreRaw(:, 2), scoreClr(:, 1), scoreClr(:, 2), ...
    scoreClrTenth(:, 1), scoreClrTenth(:, 2), ...
    'VariableNames', {'Sample', 'Strain', 'Batch', 'RawCorrectedPC1', ...
    'RawCorrectedPC2', 'ClrHalfMinPC1', 'ClrHalfMinPC2', ...
    'ClrTenthMinPC1', 'ClrTenthMinPC2'});
writetable(pcaScores, fullfile(tableDir, 'pca_scores.csv'));

pcaLoadings = table(includedMetabolite', coeffRaw(:, 1), coeffRaw(:, 2), ...
    coeffClr(:, 1), coeffClr(:, 2), coeffClrTenth(:, 1), ...
    coeffClrTenth(:, 2), 'VariableNames', {'Metabolite', ...
    'RawCorrectedPC1', 'RawCorrectedPC2', 'ClrHalfMinPC1', ...
    'ClrHalfMinPC2', 'ClrTenthMinPC1', 'ClrTenthMinPC2'});
writetable(pcaLoadings, fullfile(tableDir, 'pca_loadings.csv'));

varianceRecords = [
    makeVarianceRows("Raw blank-corrected, feature standardized", explainedRaw);
    makeVarianceRows("CLR, half-minimum censoring", explainedClr);
    makeVarianceRows("CLR, tenth-minimum censoring", explainedClrTenth)];
writetable(varianceRecords, fullfile(tableDir, 'pca_explained_variance.csv'));

%% Batch-aware global strain tests
[rawF, rawP, rawR2Strain, rawR2Batch] = balancedStrainPermutationTest( ...
    XblankZ, bacterialGroup, bacterialBatch, nPermGlobal);
[clrF, clrP, clrR2Strain, clrR2Batch] = balancedStrainPermutationTest( ...
    XclrHalf, bacterialGroup, bacterialBatch, nPermGlobal);
[clrTenthF, clrTenthP, clrTenthR2Strain, clrTenthR2Batch] = ...
    balancedStrainPermutationTest(XclrTenth, bacterialGroup, ...
    bacterialBatch, nPermGlobal);

globalTests = table([
    "Raw blank-corrected, feature standardized";
    "CLR, half-minimum censoring";
    "CLR, tenth-minimum censoring"], ...
    [rawF; clrF; clrTenthF], [rawP; clrP; clrTenthP], ...
    [rawR2Strain; clrR2Strain; clrTenthR2Strain], ...
    [rawR2Batch; clrR2Batch; clrTenthR2Batch], ...
    repmat(nPermGlobal, 3, 1), ...
    'VariableNames', {'Representation', 'StrainPseudoF', ...
    'BatchConstrainedPermutationP', 'StrainR2', 'BatchR2', ...
    'NPermutations'});
writetable(globalTests, fullfile(tableDir, 'global_strain_tests.csv'));

% Test whether the CLR strain result depends on any single retained metabolite.
leaveOneOut = table(includedMetabolite', nan(sum(include), 1), ...
    nan(sum(include), 1), nan(sum(include), 1), nan(sum(include), 1), ...
    repmat(nPermLeaveOneOut, sum(include), 1), ...
    'VariableNames', {'OmittedMetabolite', 'StrainPseudoF', ...
    'BatchConstrainedPermutationP', 'StrainR2', 'BatchR2', ...
    'NPermutations'});
for j = 1:sum(include)
    retained = true(1, sum(include));
    retained(j) = false;
    [leaveOneOut.StrainPseudoF(j), ...
        leaveOneOut.BatchConstrainedPermutationP(j), ...
        leaveOneOut.StrainR2(j), leaveOneOut.BatchR2(j)] = ...
        balancedStrainPermutationTest(XclrHalf(:, retained), ...
        bacterialGroup, bacterialBatch, nPermLeaveOneOut);
end
writetable(leaveOneOut, fullfile(tableDir, ...
    'clr_leave_one_metabolite_out_sensitivity.csv'));

%% Matched ST1-75 versus ST1-68 metabolite contrasts
focusedHalf = matchedMetaboliteContrasts(XlogHalf, Xblank, ...
    bacterialGroup, bacterialBatch, includedMetabolite, "half_minimum");
focusedTenth = matchedMetaboliteContrasts(XlogTenth, Xblank, ...
    bacterialGroup, bacterialBatch, includedMetabolite, "tenth_minimum");
focused = [focusedHalf; focusedTenth];
writetable(focused, fullfile(tableDir, 'st175_st168_matched_metabolite_effects.csv'));

%% Strain-centroid distances
distanceRaw = strainCentroidDistances(XblankZ, bacterialGroup, strainOrder, ...
    "Raw blank-corrected, feature standardized");
distanceClr = strainCentroidDistances(XclrHalf, bacterialGroup, strainOrder, ...
    "CLR, half-minimum censoring");
writetable([distanceRaw; distanceClr], ...
    fullfile(tableDir, 'strain_centroid_distances.csv'));

%% Exploratory batch-held-out PLS-DA diagnostic
% Component selection and scaling occur inside the folds. The feature screen
% and censoring floors above use all samples, so this is not an independent,
% unbiased classification-performance estimate. It is distinct from the
% descriptive full-data projection shown in Figure 3A.
classLabel = nan(numel(bacterialGroup), 1);
for g = 1:numel(strainOrder)
    classLabel(bacterialGroup == strainOrder(g)) = g;
end
assert(all(isfinite(classLabel)), 'Every bacterial sample requires a class label.');

[plsAccuracy, plsPrediction, selectedComponents] = nestedBatchPlsDa( ...
    XclrHalf, classLabel, bacterialBatch, numel(strainOrder));
permutedAccuracy = nan(nPermPls, 1);
for p = 1:nPermPls
    permutedLabel = permuteLabelsWithinBatch(classLabel, bacterialBatch);
    permutedAccuracy(p) = nestedBatchPlsDa(XclrHalf, permutedLabel, ...
        bacterialBatch, numel(strainOrder));
end
plsP = (1 + sum(permutedAccuracy >= plsAccuracy - 1e-12)) / (nPermPls + 1);

[plsTenthAccuracy, ~, ~] = nestedBatchPlsDa( ...
    XclrTenth, classLabel, bacterialBatch, numel(strainOrder));
permutedTenthAccuracy = nan(nPermPls, 1);
for p = 1:nPermPls
    permutedLabel = permuteLabelsWithinBatch(classLabel, bacterialBatch);
    permutedTenthAccuracy(p) = nestedBatchPlsDa(XclrTenth, permutedLabel, ...
        bacterialBatch, numel(strainOrder));
end
plsTenthP = (1 + sum(permutedTenthAccuracy >= plsTenthAccuracy - 1e-12)) / ...
    (nPermPls + 1);

plsSummary = table([
    "CLR, half-minimum censoring";
    "CLR, tenth-minimum censoring"], ...
    [plsAccuracy; plsTenthAccuracy], [plsP; plsTenthP], ...
    repmat(nPermPls, 2, 1), ...
    'VariableNames', {'Representation', 'LeaveOneBatchOutAccuracy', ...
    'WithinBatchPermutationP', 'NPermutations'});
writetable(plsSummary, fullfile(tableDir, 'plsda_validation_summary.csv'));

predictionTable = table(bacterialSample, bacterialGroup, bacterialBatch, ...
    strainOrder(plsPrediction)', selectedComponents, ...
    'VariableNames', {'Sample', 'ObservedStrain', 'Batch', ...
    'PredictedStrain', 'ComponentsSelectedWithoutOuterBatch'});
writetable(predictionTable, fullfile(tableDir, 'plsda_leave_one_batch_out_predictions.csv'));

confusion = confusionmat(classLabel, plsPrediction, 'Order', 1:numel(strainOrder));
confusionTable = array2table(confusion, ...
    'VariableNames', matlab.lang.makeValidName(strainOrder), ...
    'RowNames', cellstr(strainOrder));
writetable(confusionTable, fullfile(tableDir, 'plsda_confusion_matrix.csv'), ...
    'WriteRowNames', true);

permutationTable = table((1:nPermPls)', permutedAccuracy, ...
    permutedTenthAccuracy, 'VariableNames', {'Permutation', ...
    'ClrHalfMinAccuracy', 'ClrTenthMinAccuracy'});
writetable(permutationTable, fullfile(tableDir, 'plsda_permutation_distribution.csv'));

%% Figures
makePcaFigure(scoreRaw, scoreClr, explainedRaw, explainedClr, ...
    bacterialGroup, bacterialBatch, strainOrder, strainColors, ...
    rawP, rawR2Strain, rawR2Batch, clrP, clrR2Strain, clrR2Batch, figureDir);
makeStateFigure(XclrHalf, bacterialGroup, strainOrder, strainColors, ...
    includedMetabolite, focusedHalf, figureDir);
makePlsValidationFigure(permutedAccuracy, plsAccuracy, confusion, ...
    strainOrder, plsP, figureDir);

%% Human-readable summary
summaryFile = fullfile(tableDir, 'analysis_summary.txt');
fid = fopen(summaryFile, 'w');
assert(fid >= 0, 'Unable to open analysis summary for writing.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Matched-batch intracellular GC-MS analysis\n');
fprintf(fid, '=========================================\n');
fprintf(fid, 'Source: %s\n', sourceFile);
fprintf(fid, 'Matched blank rule: bacterial sample suffix i minus Media_i.\n');
fprintf(fid, 'Included metabolites: %d of %d (above matched blank in >=12/15 bacterial samples).\n', ...
    sum(include), numel(metabolite));
fprintf(fid, 'Global raw corrected strain effect: pseudo-F %.4f, p %.6f, strain R2 %.4f, batch R2 %.4f.\n', ...
    rawF, rawP, rawR2Strain, rawR2Batch);
fprintf(fid, 'Global CLR strain effect: pseudo-F %.4f, p %.6f, strain R2 %.4f, batch R2 %.4f.\n', ...
    clrF, clrP, clrR2Strain, clrR2Batch);
fprintf(fid, 'CLR tenth-min sensitivity: pseudo-F %.4f, p %.6f, strain R2 %.4f, batch R2 %.4f.\n', ...
    clrTenthF, clrTenthP, clrTenthR2Strain, clrTenthR2Batch);
fprintf(fid, ['CLR leave-one-metabolite-out sensitivity: strain R2 range ' ...
    '%.4f-%.4f; permutation p range %.6f-%.6f.\n'], ...
    min(leaveOneOut.StrainR2), max(leaveOneOut.StrainR2), ...
    min(leaveOneOut.BatchConstrainedPermutationP), ...
    max(leaveOneOut.BatchConstrainedPermutationP));
fprintf(fid, 'Nested leave-one-batch-out PLS-DA accuracy: %.4f; within-batch permutation p %.6f.\n', ...
    plsAccuracy, plsP);
fprintf(fid, 'Tenth-min PLS-DA sensitivity accuracy: %.4f; permutation p %.6f.\n', ...
    plsTenthAccuracy, plsTenthP);
fprintf(fid, 'PC1+PC2 variance, raw corrected: %.2f%%; CLR: %.2f%%.\n', ...
    sum(explainedRaw(1:2)), sum(explainedClr(1:2)));
fprintf(fid, 'No mechanistic or flux interpretation is made from intracellular pool sizes.\n');

disp(fileread(summaryFile));

%% Local functions
function standardized = standardizeColumns(values)
means = mean(values, 1);
scales = std(values, 0, 1);
scales(scales == 0) = 1;
standardized = (values - means) ./ scales;
end

function rows = makeVarianceRows(representation, explained)
n = min(5, numel(explained));
rows = table(repmat(string(representation), n, 1), (1:n)', explained(1:n), ...
    'VariableNames', {'Representation', 'PrincipalComponent', ...
    'ExplainedVariancePercent'});
end

function [pseudoF, pValue, strainR2, batchR2] = ...
        balancedStrainPermutationTest(values, strain, batch, nPerm)
grand = mean(values, 1);
totalSS = sum((values - grand).^2, 'all');
uniqueStrain = unique(strain, 'stable');
uniqueBatch = unique(batch, 'stable');

strainSS = 0;
for g = uniqueStrain'
    mask = strain == g;
    strainSS = strainSS + sum(mask) * sum((mean(values(mask, :), 1) - grand).^2);
end
batchSS = 0;
for b = uniqueBatch'
    mask = batch == b;
    batchSS = batchSS + sum(mask) * sum((mean(values(mask, :), 1) - grand).^2);
end
residualSS = max(0, totalSS - strainSS - batchSS);
strainDf = numel(uniqueStrain) - 1;
residualDf = size(values, 1) - numel(uniqueStrain) - numel(uniqueBatch) + 1;
pseudoF = (strainSS / strainDf) / (residualSS / residualDf);
strainR2 = strainSS / totalSS;
batchR2 = batchSS / totalSS;

permutedF = nan(nPerm, 1);
for p = 1:nPerm
    permutedStrain = strain;
    for b = uniqueBatch'
        rows = find(batch == b);
        permutedStrain(rows) = strain(rows(randperm(numel(rows))));
    end
    permutedSS = 0;
    for g = uniqueStrain'
        mask = permutedStrain == g;
        permutedSS = permutedSS + sum(mask) * ...
            sum((mean(values(mask, :), 1) - grand).^2);
    end
    permutedResidual = max(eps, totalSS - permutedSS - batchSS);
    permutedF(p) = (permutedSS / strainDf) / ...
        (permutedResidual / residualDf);
end
pValue = (1 + sum(permutedF >= pseudoF - 1e-12)) / (nPerm + 1);
end

function result = matchedMetaboliteContrasts(logValues, rawCorrected, ...
        strain, batch, metabolite, floorMethod)
rows75 = find(strain == "ST1-75");
rows68 = find(strain == "ST1-68");
[~, order75] = sort(batch(rows75));
[~, order68] = sort(batch(rows68));
rows75 = rows75(order75);
rows68 = rows68(order68);
assert(isequal(batch(rows75), batch(rows68)), ...
    'Focused strain rows must be matched by batch.');

nMetabolite = numel(metabolite);
meanDifference = nan(nMetabolite, 1);
sdDifference = nan(nMetabolite, 1);
pairedP = nan(nMetabolite, 1);
signFlipP = nan(nMetabolite, 1);
n75Nonpositive = nan(nMetabolite, 1);
n68Nonpositive = nan(nMetabolite, 1);
for j = 1:nMetabolite
    difference = logValues(rows75, j) - logValues(rows68, j);
    meanDifference(j) = mean(difference);
    sdDifference(j) = std(difference, 0);
    if any(abs(difference - difference(1)) > 1e-12) || abs(difference(1)) > 1e-12
        [~, pairedP(j)] = ttest(difference, 0);
    end
    signFlipP(j) = exactSignFlipP(difference);
    n75Nonpositive(j) = sum(rawCorrected(rows75, j) <= 0);
    n68Nonpositive(j) = sum(rawCorrected(rows68, j) <= 0);
end
qValue = bhAdjust(pairedP);
result = table(repmat(string(floorMethod), nMetabolite, 1), metabolite', ...
    meanDifference, 2.^meanDifference, sdDifference, pairedP, qValue, ...
    signFlipP, n75Nonpositive, n68Nonpositive, ...
    'VariableNames', {'CensoringMethod', 'Metabolite', ...
    'MeanLog2Difference_ST175MinusST168', ...
    'GeometricMeanRatio_ST175OverST168', 'SDLog2Difference', ...
    'PairedTP', 'BenjaminiHochbergQ', 'ExactSignFlipP', ...
    'N_ST175_AtOrBelowBlank', 'N_ST168_AtOrBelowBlank'});
end

function pValue = exactSignFlipP(difference)
n = numel(difference);
observed = abs(mean(difference));
allMean = nan(2^n, 1);
for k = 0:(2^n - 1)
    signValue = 2 * bitget(k, 1:n) - 1;
    allMean(k + 1) = abs(mean(difference .* signValue'));
end
pValue = mean(allMean >= observed - 1e-12);
end

function adjusted = bhAdjust(pValue)
adjusted = nan(size(pValue));
finite = find(isfinite(pValue));
if isempty(finite)
    return;
end
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

function result = strainCentroidDistances(values, strain, strainOrder, representation)
recordA = strings(0, 1);
recordB = strings(0, 1);
distance = zeros(0, 1);
correlation = zeros(0, 1);
for a = 1:(numel(strainOrder) - 1)
    centroidA = mean(values(strain == strainOrder(a), :), 1);
    for b = (a + 1):numel(strainOrder)
        centroidB = mean(values(strain == strainOrder(b), :), 1);
        recordA(end + 1, 1) = strainOrder(a); %#ok<AGROW>
        recordB(end + 1, 1) = strainOrder(b); %#ok<AGROW>
        distance(end + 1, 1) = norm(centroidA - centroidB); %#ok<AGROW>
        correlation(end + 1, 1) = corr(centroidA', centroidB', ...
            'Type', 'Spearman'); %#ok<AGROW>
    end
end
result = table(repmat(string(representation), numel(distance), 1), ...
    recordA, recordB, distance, correlation, 'VariableNames', ...
    {'Representation', 'StrainA', 'StrainB', 'EuclideanDistance', ...
    'SpearmanCorrelation'});
end

function [accuracy, prediction, selectedComponents] = nestedBatchPlsDa( ...
        values, label, batch, nClass)
uniqueBatch = unique(batch, 'stable');
prediction = nan(size(label));
selectedComponents = nan(size(label));
for outerBatch = uniqueBatch'
    outerTest = batch == outerBatch;
    outerTrain = ~outerTest;
    trainingBatch = unique(batch(outerTrain), 'stable');
    maximumComponent = min([3, sum(outerTrain) - 1, size(values, 2)]);
    componentAccuracy = nan(maximumComponent, 1);
    for component = 1:maximumComponent
        innerPrediction = nan(sum(outerTrain), 1);
        trainingRows = find(outerTrain);
        for innerBatch = trainingBatch'
            innerTestGlobal = outerTrain & batch == innerBatch;
            innerTrainGlobal = outerTrain & batch ~= innerBatch;
            predicted = fitPredictPls(values(innerTrainGlobal, :), ...
                label(innerTrainGlobal), values(innerTestGlobal, :), ...
                component, nClass);
            [~, location] = ismember(find(innerTestGlobal), trainingRows);
            innerPrediction(location) = predicted;
        end
        componentAccuracy(component) = mean(innerPrediction == label(outerTrain));
    end
    best = find(componentAccuracy == max(componentAccuracy), 1, 'first');
    prediction(outerTest) = fitPredictPls(values(outerTrain, :), ...
        label(outerTrain), values(outerTest, :), best, nClass);
    selectedComponents(outerTest) = best;
end
accuracy = mean(prediction == label);
end

function prediction = fitPredictPls(trainX, trainLabel, testX, component, nClass)
trainMean = mean(trainX, 1);
trainScale = std(trainX, 0, 1);
trainScale(trainScale == 0) = 1;
trainZ = (trainX - trainMean) ./ trainScale;
testZ = (testX - trainMean) ./ trainScale;
component = min([component, size(trainZ, 1) - 1, rank(trainZ), size(trainZ, 2)]);
component = max(component, 1);
response = zeros(size(trainZ, 1), nClass);
response(sub2ind(size(response), (1:size(response, 1))', trainLabel)) = 1;
[~, ~, ~, ~, beta] = plsregress(trainZ, response, component);
score = [ones(size(testZ, 1), 1), testZ] * beta;
[~, prediction] = max(score, [], 2);
end

function permuted = permuteLabelsWithinBatch(label, batch)
permuted = label;
for b = unique(batch, 'stable')'
    rows = find(batch == b);
    permuted(rows) = label(rows(randperm(numel(rows))));
end
end

function makePcaFigure(scoreRaw, scoreClr, explainedRaw, explainedClr, ...
        strain, batch, strainOrder, colors, rawP, rawR2Strain, rawR2Batch, ...
        clrP, clrR2Strain, clrR2Batch, figureDir)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1420, 610]);
layout = tiledlayout(figureHandle, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
ax1 = nexttile(layout);
plotPcaPanel(ax1, scoreRaw, strain, batch, strainOrder, colors);
title(ax1, sprintf(['A  Blank-corrected intracellular pools\n' ...
    'strain R^2=%.2f, p=%.4f; batch R^2=%.2f'], ...
    rawR2Strain, rawP, rawR2Batch), 'FontWeight', 'bold');
xlabel(ax1, sprintf('PC1 (%.1f%%)', explainedRaw(1)));
ylabel(ax1, sprintf('PC2 (%.1f%%)', explainedRaw(2)));

ax2 = nexttile(layout);
plotPcaPanel(ax2, scoreClr, strain, batch, strainOrder, colors);
title(ax2, sprintf(['B  Relative intracellular pool composition\n' ...
    'strain R^2=%.2f, p=%.4f; batch R^2=%.2f'], ...
    clrR2Strain, clrP, clrR2Batch), 'FontWeight', 'bold');
xlabel(ax2, sprintf('PC1 (%.1f%%)', explainedClr(1)));
ylabel(ax2, sprintf('PC2 (%.1f%%)', explainedClr(2)));

legendHandle = gobjects(numel(strainOrder), 1);
for g = 1:numel(strainOrder)
    legendHandle(g) = plot(ax2, nan, nan, 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', colors(g, :), 'MarkerEdgeColor', 'w');
end
legend(ax2, legendHandle, strainOrder, 'Location', 'southoutside', ...
    'NumColumns', 5, 'Box', 'off');
sgtitle(layout, ['Matched-blank intracellular GC-MS supports ' ...
    'strain-dependent metabolic states'], 'FontSize', 16, 'FontWeight', 'bold');
exportgraphics(figureHandle, fullfile(figureDir, ...
    'intracellular_blank_corrected_pca.png'), 'Resolution', 300);
close(figureHandle);
end

function plotPcaPanel(ax, score, strain, batch, strainOrder, colors)
hold(ax, 'on');
for g = 1:numel(strainOrder)
    rows = find(strain == strainOrder(g));
    [~, order] = sort(batch(rows));
    rows = rows(order);
    plot(ax, score(rows, 1), score(rows, 2), '-', 'Color', ...
        [colors(g, :), 0.28], 'LineWidth', 1.0);
    scatter(ax, score(rows, 1), score(rows, 2), 92, colors(g, :), ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.9);
    for i = 1:numel(rows)
        text(ax, score(rows(i), 1), score(rows(i), 2), ...
            string(batch(rows(i))), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', 'w', ...
            'FontSize', 8, 'FontWeight', 'bold');
    end
end
xline(ax, 0, 'Color', [0.85, 0.85, 0.85]);
yline(ax, 0, 'Color', [0.85, 0.85, 0.85]);
axis(ax, 'equal');
box(ax, 'off');
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1);
end

function makeStateFigure(clrValues, strain, strainOrder, colors, ...
        metabolite, focused, figureDir)
groupMean = nan(numel(strainOrder), numel(metabolite));
for g = 1:numel(strainOrder)
    groupMean(g, :) = mean(clrValues(strain == strainOrder(g), :), 1);
end
featureVariation = var(groupMean, 0, 1);
[~, featureOrder] = sort(featureVariation, 'descend');
heat = standardizeColumns(groupMean)';
heat = heat(featureOrder, :);

primary = focused(focused.CensoringMethod == "half_minimum", :);
[~, effectOrder] = sort(primary.MeanLog2Difference_ST175MinusST168);
primary = primary(effectOrder, :);
tCritical = tinv(0.975, 2);
halfWidth = tCritical .* primary.SDLog2Difference ./ sqrt(3);

figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1450, 760]);
layout = tiledlayout(figureHandle, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
ax1 = nexttile(layout);
imagesc(ax1, heat, [-2, 2]);
colormap(ax1, blueRedMap(256));
colorbarHandle = colorbar(ax1, 'Location', 'eastoutside');
colorbarHandle.Label.String = 'Z score across strain means';
set(ax1, 'XTick', 1:numel(strainOrder), 'XTickLabel', strainOrder, ...
    'XTickLabelRotation', 28, 'YTick', 1:numel(metabolite), ...
    'YTickLabel', metabolite(featureOrder), 'FontName', 'Arial', ...
    'FontSize', 10);
title(ax1, 'A  Relative intracellular pool states', ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left');
xlabel(ax1, 'Strain');

ax2 = nexttile(layout);
hold(ax2, 'on');
xline(ax2, 0, '--', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 1.0);
for i = 1:height(primary)
    color = [0.38, 0.44, 0.50];
    if primary.BenjaminiHochbergQ(i) < 0.05
        color = colors(1, :);
    end
    errorbar(ax2, primary.MeanLog2Difference_ST175MinusST168(i), i, ...
        halfWidth(i), halfWidth(i), 'horizontal', 'o', ...
        'Color', color, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.3, 'CapSize', 4, 'MarkerSize', 6);
end
set(ax2, 'YTick', 1:height(primary), 'YTickLabel', primary.Metabolite, ...
    'YDir', 'reverse', 'FontName', 'Arial', 'FontSize', 10);
xlabel(ax2, 'Matched mean log_2 ratio, ST1-75 / ST1-68');
title(ax2, 'B  Focused intracellular pool differences', ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left');
box(ax2, 'off');
sgtitle(layout, ['Blank-corrected intracellular states differ across strains ' ...
    'without resolving to a single metabolite'], 'FontSize', 15, ...
    'FontWeight', 'bold');
exportgraphics(figureHandle, fullfile(figureDir, ...
    'intracellular_state_heatmap_and_focused_effects.png'), 'Resolution', 300);
close(figureHandle);
end

function makePlsValidationFigure(permutedAccuracy, observedAccuracy, ...
        confusion, strainOrder, pValue, figureDir)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1200, 500]);
layout = tiledlayout(figureHandle, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
ax1 = nexttile(layout);
histogram(ax1, permutedAccuracy, 'BinEdges', (-0.5:1:15.5) ./ 15, ...
    'FaceColor', [0.65, 0.70, 0.75], 'EdgeColor', 'w');
xline(ax1, observedAccuracy, '-', 'Color', [0.13, 0.38, 0.63], ...
    'LineWidth', 2.5, 'Label', sprintf('observed = %.0f%%', ...
    observedAccuracy * 100), 'LabelVerticalAlignment', 'middle');
xlabel(ax1, 'Nested leave-one-batch-out accuracy');
ylabel(ax1, 'Within-batch label permutations');
xlim(ax1, [0, 1]);
title(ax1, sprintf('A  PLS-DA validation, permutation p=%.3f', pValue), ...
    'FontWeight', 'bold');
box(ax1, 'off');

ax2 = nexttile(layout);
imagesc(ax2, confusion);
colormap(ax2, parula(256));
colorbar(ax2);
set(ax2, 'XTick', 1:numel(strainOrder), 'XTickLabel', strainOrder, ...
    'XTickLabelRotation', 30, 'YTick', 1:numel(strainOrder), ...
    'YTickLabel', strainOrder, 'FontName', 'Arial', 'FontSize', 10);
xlabel(ax2, 'Predicted strain');
ylabel(ax2, 'Observed strain');
title(ax2, 'B  Held-out-batch confusion matrix', 'FontWeight', 'bold');
for row = 1:size(confusion, 1)
    for column = 1:size(confusion, 2)
        text(ax2, column, row, string(confusion(row, column)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', 'w', 'FontWeight', 'bold');
    end
end
sgtitle(layout, 'Supervised strain separation requires held-out-batch validation', ...
    'FontSize', 15, 'FontWeight', 'bold');
exportgraphics(figureHandle, fullfile(figureDir, ...
    'intracellular_plsda_validation.png'), 'Resolution', 300);
close(figureHandle);
end

function map = blueRedMap(n)
blue = [0.08, 0.28, 0.52];
white = [0.98, 0.98, 0.97];
red = [0.65, 0.05, 0.15];
half = floor(n / 2);
first = [linspace(blue(1), white(1), half)', ...
    linspace(blue(2), white(2), half)', ...
    linspace(blue(3), white(3), half)'];
second = [linspace(white(1), red(1), n - half)', ...
    linspace(white(2), red(2), n - half)', ...
    linspace(white(3), red(3), n - half)'];
map = [first; second];
end
