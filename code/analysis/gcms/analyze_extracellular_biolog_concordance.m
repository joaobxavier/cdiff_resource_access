%% Extracellular CoRe structure and concordance with audited PM1 calls
% Compares the five-strain extracellular GC-MS footprint with Vishwas's
% authoritative 700-call PM1 matrix. R20291 is retained in all GC-MS analyses
% but is absent from the PM1 panel and is therefore excluded only from direct
% cross-assay comparisons.

clear;

rng(20260821, 'twister');

scriptPath = mfilename('fullpath');
localRoot = fileparts(fileparts(fileparts(scriptPath)));
repoRoot = fileparts(localRoot);
gcmsFile = package_input('secretome.xlsx');
pm1File = fullfile(package_root(),'results','tables','biolog','Biolog_growth_matrix.csv');
moaFile = package_input('moas.xlsx');
tableDir = fullfile(repoRoot, 'results', 'tables', 'gcms', ...
    'gcms_biolog_concordance');
figureDir = fullfile(repoRoot, 'results', 'diagnostic_figures', 'gcms', ...
    'gcms_biolog_concordance');
if ~isfolder(tableDir); mkdir(tableDir); end
if ~isfolder(figureDir); mkdir(figureDir); end

strainOrder = ["ST1-75", "ST1-68", "ST1-6", "VPI10463", "R20291"];
crossStrains = ["ST1-75", "ST1-68", "ST1-6", "VPI10463"];
pm1Columns = ["ST1_75", "ST1_68", "ST1_6", "VPI"];
strainColors = [
    36, 104, 162;
    60, 157, 86;
    142, 90, 169;
    217, 71, 63;
    228, 154, 47] ./ 255;

%% Load and verify the extracellular matched-batch design
source = readtable(gcmsFile, 'VariableNamingRule', 'preserve');
sampleID = string(source{:, 1});
metabolite = string(source.Properties.VariableNames(2:end));
metabolite = regexprep(metabolite, '_\d+Results$', '');
rawSignal = double(source{:, 2:end});
group = extractBefore(sampleID, '_');
group(group == "VPI") = "VPI10463";
batch = str2double(extractAfter(sampleID, '_'));

expectedGroups = ["Media", strainOrder];
for g = expectedGroups
    assert(sum(group == g) == 3, 'Expected three samples for %s.', g);
    assert(isequal(sort(batch(group == g))', 1:3), ...
        'Expected batch suffixes 1, 2, and 3 for %s.', g);
end
assert(all(rawSignal >= 0, 'all') && all(isfinite(rawSignal), 'all'), ...
    'GC-MS signals must be finite and nonnegative.');

bacterialRow = group ~= "Media";
bacterialSample = sampleID(bacterialRow);
bacterialGroup = group(bacterialRow);
bacterialBatch = batch(bacterialRow);

% Thirteen metabolites are positive in at least two of three samples in every
% group, including media; this defines the continuous extracellular subset.
continuous = true(1, numel(metabolite));
for g = expectedGroups
    continuous = continuous & sum(rawSignal(group == g, :) > 0, 1) >= 2;
end
assert(sum(continuous) == 13, ...
    'The audited continuous extracellular subset should contain 13 metabolites.');

% Half-minimum replacement is used for log-scale displays, tests and distance
% diagnostics. Raw media-range calls below remain imputation-free.
workingSignal = rawSignal;
for j = 1:numel(metabolite)
    positive = rawSignal(rawSignal(:, j) > 0, j);
    floorValue = 0.5 * min(positive);
    workingSignal(workingSignal(:, j) <= 0, j) = floorValue;
end

matchedLog2Ratio = nan(sum(bacterialRow), numel(metabolite));
matchedMediaSignal = nan(sum(bacterialRow), numel(metabolite));
rawBacterialSignal = rawSignal(bacterialRow, :);
for i = 1:sum(bacterialRow)
    sourceRow = find(bacterialRow);
    sourceRow = sourceRow(i);
    mediaRow = find(group == "Media" & batch == batch(sourceRow), 1);
    matchedMediaSignal(i, :) = rawSignal(mediaRow, :);
    matchedLog2Ratio(i, :) = log2(workingSignal(sourceRow, :)) - ...
        log2(workingSignal(mediaRow, :));
end

groupMeanLog2Ratio = nan(numel(strainOrder), numel(metabolite));
rangeState = nan(numel(strainOrder), numel(metabolite));
mediaMin = min(rawSignal(group == "Media", :), [], 1);
mediaMax = max(rawSignal(group == "Media", :), [], 1);
for g = 1:numel(strainOrder)
    bacterialRows = bacterialGroup == strainOrder(g);
    sourceRows = group == strainOrder(g);
    groupMeanLog2Ratio(g, :) = mean(matchedLog2Ratio(bacterialRows, :), 1);
    strainMin = min(rawSignal(sourceRows, :), [], 1);
    strainMax = max(rawSignal(sourceRows, :), [], 1);
    rangeState(g, strainMax < mediaMin) = -1;
    rangeState(g, strainMin > mediaMax) = 1;
    rangeState(g, isnan(rangeState(g, :))) = 0;
end

% Export the complete source-to-figure mapping in long form. Raw signals and
% matched media values remain untouched; the log ratio uses the documented
% feature-specific half-minimum replacement only when a source value is zero.
longRows = table;
for i = 1:numel(bacterialSample)
    strainIndex = find(strainOrder == bacterialGroup(i), 1);
    for j = 1:numel(metabolite)
        newRow = table(bacterialSample(i), bacterialGroup(i), ...
            bacterialBatch(i), metabolite(j), rawBacterialSignal(i, j), ...
            matchedMediaSignal(i, j), matchedLog2Ratio(i, j), ...
            continuous(j), rangeState(strainIndex, j), ...
            'VariableNames', {'Sample', 'Strain', 'Batch', 'Metabolite', ...
            'RawProcessedSignal', 'MatchedMediaRawProcessedSignal', ...
            'MatchedLog2ResidualRatio', 'ContinuousFeature', ...
            'StrainMediaRangeState'});
        longRows = [longRows; newRow]; %#ok<AGROW>
    end
end
writetable(longRows, fullfile(tableDir, ...
    'extracellular_five_strain_long_form.csv'));

%% Load the authoritative audited PM1 matrix
pm1 = readtable(pm1File, 'VariableNamingRule', 'preserve', 'TextType', 'string');
pm1 = pm1(pm1.Metabolites ~= "Negative Control", :);
allCalls = double(table2array(pm1(:, 2:end)));
assert(height(pm1) == 95 && sum(allCalls, 'all') == 700, ...
    'PM1 input is not the authoritative audited 700-call matrix.');

pm1Call = nan(height(pm1), numel(crossStrains));
for s = 1:numel(crossStrains)
    pm1Call(:, s) = double(pm1{:, pm1Columns(s)});
end
pm1Breadth = sum(pm1Call, 1)';
assert(isequal(pm1Breadth', [82, 24, 76, 23]), ...
    'Focused PM1 breadth no longer matches the audited matrix.');

% The source repository supplies chemical-class annotations for the PM1 wells.
% Record amino-acid-class breadth because the extracellular GC-MS panel is
% amino-acid centered. This is a class-level comparison, not evidence that an
% individual binary PM1 call predicts endpoint depletion in complex medium.
moa = readtable(moaFile, 'VariableNamingRule', 'preserve', ...
    'TextType', 'string');
[isAnnotated, moaRow] = ismember(pm1.Metabolites, moa.Chemical);
assert(all(isAnnotated), 'Every audited PM1 substrate should have a class annotation.');
pm1Class = moa.MoA(moaRow);
aminoAcidRow = contains(lower(pm1Class), 'amino acid');
assert(sum(aminoAcidRow) == 16, ...
    'Expected 16 amino-acid-class substrates in the PM1 annotation.');
pm1AminoAcidBreadth = sum(pm1Call(aminoAcidRow, :), 1)';
aminoAcidBreadthTable = table(crossStrains', pm1AminoAcidBreadth, ...
    pm1AminoAcidBreadth ./ sum(aminoAcidRow), ...
    repmat(sum(aminoAcidRow), numel(crossStrains), 1), ...
    'VariableNames', {'Strain', 'Pm1AminoAcidPositiveCallBreadth', ...
    'Pm1AminoAcidPositiveCallFraction', ...
    'NAnnotatedPm1AminoAcidSubstrates'});
writetable(aminoAcidBreadthTable, fullfile(tableDir, ...
    'pm1_amino_acid_class_breadth.csv'));

%% Explicit direct-metabolite mapping
% The GC-MS labels do not encode stereochemistry. The L-form PM1 well is used
% as the primary mapping because these are conventional amino-acid labels; an
% any-stereoisomer sensitivity is retained where a D-form PM1 well exists.
mapGcms = ["Lactate"; "Alanine"; "Aspartate"; "Glutamate"; ...
    "Proline"; "Serine"; "Threonine"];
mapPm1L = ["L-Lactic acid"; "L-Alanine"; "L-Aspartic acid"; ...
    "L-Glutamic acid"; "L-Proline"; "L-Serine"; "L-Threonine"];
mapPm1D = [""; "D-Alanine"; "D-Aspartic acid"; ""; ""; ...
    "D-Serine"; "D-Threonine"];
mappingNote = repmat("GC-MS stereochemistry not encoded; L-form PM1 is primary", ...
    numel(mapGcms), 1);
mapping = table(mapGcms, mapPm1L, mapPm1D, mappingNote, ...
    'VariableNames', {'GcmsMetabolite', 'PrimaryPm1Substrate', ...
    'AlternatePm1Stereoisomer', 'MappingNote'});
writetable(mapping, fullfile(tableDir, 'gcms_pm1_metabolite_mapping.csv'));

nMapped = numel(mapGcms);
primaryMappedCall = nan(nMapped, numel(crossStrains));
anyStereoCall = nan(nMapped, numel(crossStrains));
mappedLog2Ratio = nan(nMapped, numel(crossStrains));
mappedRangeState = nan(nMapped, numel(crossStrains));
for m = 1:nMapped
    pm1Row = find(pm1.Metabolites == mapPm1L(m), 1);
    gcmsColumn = find(metabolite == mapGcms(m), 1);
    assert(~isempty(pm1Row) && ~isempty(gcmsColumn), ...
        'Mapped substrate is absent from a source table.');
    primaryMappedCall(m, :) = pm1Call(pm1Row, :);
    anyStereoCall(m, :) = primaryMappedCall(m, :);
    if strlength(mapPm1D(m)) > 0
        alternateRow = find(pm1.Metabolites == mapPm1D(m), 1);
        anyStereoCall(m, :) = max(anyStereoCall(m, :), ...
            pm1Call(alternateRow, :));
    end
    for s = 1:numel(crossStrains)
        gcmsGroup = find(strainOrder == crossStrains(s), 1);
        mappedLog2Ratio(m, s) = groupMeanLog2Ratio(gcmsGroup, gcmsColumn);
        mappedRangeState(m, s) = rangeState(gcmsGroup, gcmsColumn);
    end
end

crossRows = table;
for m = 1:nMapped
    for s = 1:numel(crossStrains)
        newRow = table(mapGcms(m), mapPm1L(m), crossStrains(s), ...
            primaryMappedCall(m, s), anyStereoCall(m, s), ...
            mappedLog2Ratio(m, s), mappedRangeState(m, s), ...
            'VariableNames', {'GcmsMetabolite', 'Pm1Substrate', 'Strain', ...
            'PrimaryLFormCall', 'AnyStereoisomerCall', ...
            'MeanMatchedLog2ResidualRatio', 'RawMediaRangeState'});
        crossRows = [crossRows; newRow]; %#ok<AGROW>
    end
end
writetable(crossRows, fullfile(tableDir, ...
    'gcms_pm1_direct_overlap_comparison.csv'));

%% Cross-assay concordance tests
% Test whether PM1-positive strain/substrate pairs show deeper extracellular
% depletion in the matched GC-MS assay. Calls are permuted among strains within
% each substrate, preserving the number of positive calls per substrate.
[primaryEffect, primaryP, primaryInformative] = ...
    callDepletionConcordance(primaryMappedCall, -mappedLog2Ratio, 9999);
[stereoEffect, stereoP, stereoInformative] = ...
    callDepletionConcordance(anyStereoCall, -mappedLog2Ratio, 9999);
concordanceTests = table([
    "Primary L-form mapping";
    "Any-stereoisomer sensitivity"], ...
    [primaryEffect; stereoEffect], [primaryP; stereoP], ...
    [primaryInformative; stereoInformative], repmat(9999, 2, 1), ...
    'VariableNames', {'Mapping', ...
    'MeanDepletionDifference_CallPositiveMinusNegative', ...
    'WithinSubstratePermutationP', 'NInformativeSubstrates', ...
    'NPermutations'});
writetable(concordanceTests, fullfile(tableDir, ...
    'gcms_pm1_direct_concordance_tests.csv'));

% Compare whole-profile distances across the four shared strains. GC-MS uses
% feature-standardized mean matched log2 residual ratios for the 13 continuous
% metabolites; PM1 uses Jaccard distance across all 95 calls.
crossGroupIndex = arrayfun(@(x) find(strainOrder == x, 1), crossStrains);
gcmsProfile = groupMeanLog2Ratio(crossGroupIndex, continuous);
gcmsProfile = (gcmsProfile - mean(gcmsProfile, 1)) ./ std(gcmsProfile, 0, 1);
gcmsDistance = squareform(pdist(gcmsProfile, 'euclidean'));
pm1Distance = nan(numel(crossStrains));
for a = 1:numel(crossStrains)
    for b = 1:numel(crossStrains)
        unionCount = sum(pm1Call(:, a) | pm1Call(:, b));
        intersectionCount = sum(pm1Call(:, a) & pm1Call(:, b));
        pm1Distance(a, b) = 1 - intersectionCount / unionCount;
    end
end
[mantelRho, mantelP] = exactMantel(gcmsDistance, pm1Distance);

nPairs = nchoosek(numel(crossStrains), 2);
pairA = strings(nPairs, 1);
pairB = strings(nPairs, 1);
gcmsPairDistance = zeros(nPairs, 1);
pm1PairDistance = zeros(nPairs, 1);
pairIndex = 0;
for a = 1:(numel(crossStrains) - 1)
    for b = (a + 1):numel(crossStrains)
        pairIndex = pairIndex + 1;
        pairA(pairIndex) = crossStrains(a);
        pairB(pairIndex) = crossStrains(b);
        gcmsPairDistance(pairIndex) = gcmsDistance(a, b);
        pm1PairDistance(pairIndex) = pm1Distance(a, b);
    end
end
distanceTable = table(pairA, pairB, gcmsPairDistance, pm1PairDistance, ...
    repmat(mantelRho, numel(pairA), 1), repmat(mantelP, numel(pairA), 1), ...
    'VariableNames', {'StrainA', 'StrainB', ...
    'GcmsStandardizedProfileDistance', 'Pm1JaccardDistance', ...
    'ExactMantelSpearmanRho', 'ExactMantelP'});
writetable(distanceTable, fullfile(tableDir, ...
    'gcms_pm1_profile_distance_concordance.csv'));

%% Breadth and ST1-6 metabolic differences
gcmsDepletionBreadth = nan(numel(crossStrains), 1);
for s = 1:numel(crossStrains)
    g = find(strainOrder == crossStrains(s), 1);
    gcmsDepletionBreadth(s) = sum(rangeState(g, :) == -1);
end
breadthTable = table(crossStrains', pm1Breadth, ...
    pm1Breadth ./ 95, gcmsDepletionBreadth, gcmsDepletionBreadth ./ 18, ...
    'VariableNames', {'Strain', 'Pm1PositiveCallBreadth', ...
    'Pm1PositiveCallFraction', 'GcmsBelowMediaRangeCount', ...
    'GcmsBelowMediaRangeFraction'});
writetable(breadthTable, fullfile(tableDir, ...
    'gcms_pm1_breadth_comparison.csv'));

% Batch-constrained global strain tests for each continuous extracellular
% metabolite test strain effects; the heat map is a diagnostic output only.
% Figure 3B instead displays the independently fitted extracellular PLS-DA.
globalF = nan(sum(continuous), 1);
globalP = nan(sum(continuous), 1);
globalR2 = nan(sum(continuous), 1);
continuousColumns = find(continuous);
for j = 1:sum(continuous)
    [globalF(j), globalP(j), globalR2(j)] = balancedGroupPermutationTest( ...
        matchedLog2Ratio(:, continuousColumns(j)), ...
        bacterialGroup, bacterialBatch, 9999);
end
globalQ = bhAdjust(globalP);
globalTests = table(metabolite(continuous)', globalF, globalP, globalQ, ...
    globalR2, repmat(9999, sum(continuous), 1), ...
    'VariableNames', {'Metabolite', 'StrainPseudoF', ...
    'BatchConstrainedPermutationP', 'BenjaminiHochbergQ', 'StrainR2', ...
    'NPermutations'});
writetable(globalTests, fullfile(tableDir, ...
    'extracellular_matched_global_strain_tests.csv'));

% Pairwise ST1-75/ST1-68 analysis of extracellular remodeling.
% The same-suffix comparison is the confirmed matched-batch design. Because
% matched media cancel within a batch, differences in matched log2 ratios equal
% differences in the corresponding half-minimum log2 processed signals.
rows75 = find(bacterialGroup == "ST1-75");
rows68 = find(bacterialGroup == "ST1-68");
[~, order75] = sort(bacterialBatch(rows75));
[~, order68] = sort(bacterialBatch(rows68));
rows75 = rows75(order75);
rows68 = rows68(order68);
assert(isequal(bacterialBatch(rows75), bacterialBatch(rows68)), ...
    'ST1-75 and ST1-68 batch suffixes do not align.');
focusedDifference = matchedLog2Ratio(rows75, :) - ...
    matchedLog2Ratio(rows68, :);
focusedP = nan(numel(metabolite), 1);
focusedSignFlipP = nan(numel(metabolite), 1);
for j = 1:numel(metabolite)
    difference = focusedDifference(:, j);
    if ~all(abs(difference) < 1e-12)
        [~, focusedP(j)] = ttest(difference, 0);
        focusedSignFlipP(j) = exactSignFlipP(difference);
    end
end
focusedQ = bhAdjust(focusedP);
focusedTests = table(metabolite', ...
    repmat(numel(rows75), numel(metabolite), 1), ...
    mean(focusedDifference, 1)', 2 .^ mean(focusedDifference, 1)', ...
    sum(focusedDifference < 0, 1)', focusedP, focusedQ, ...
    focusedSignFlipP, repmat(sum(isfinite(focusedP)), ...
    numel(metabolite), 1), ...
    'VariableNames', {'Metabolite', 'NMatchedBatches', ...
    'MeanLog2ResidualDifference_ST1_75_Minus_ST1_68', ...
    'GeometricResidualSignalRatio_ST1_75_Over_ST1_68', ...
    'NMatchedBatches_ST1_75_Lower', 'PairedTP', ...
    'BenjaminiHochbergQAcrossNonconstantMetabolites', ...
    'ExactTwoSidedSignFlipP', 'NNonconstantMetabolitesInFdrFamily'});
writetable(focusedTests, fullfile(tableDir, ...
    'extracellular_st175_st168_matched_metabolite_tests.csv'));

focusedPairRows = table;
for j = 1:numel(metabolite)
    for i = 1:numel(rows75)
        newRow = table(metabolite(j), bacterialBatch(rows75(i)), ...
            bacterialSample(rows75(i)), bacterialSample(rows68(i)), ...
            rawBacterialSignal(rows75(i), j), ...
            rawBacterialSignal(rows68(i), j), ...
            matchedLog2Ratio(rows75(i), j), ...
            matchedLog2Ratio(rows68(i), j), focusedDifference(i, j), ...
            'VariableNames', {'Metabolite', 'Batch', 'ST1_75_Sample', ...
            'ST1_68_Sample', 'ST1_75_RawProcessedSignal', ...
            'ST1_68_RawProcessedSignal', 'ST1_75_MatchedLog2Ratio', ...
            'ST1_68_MatchedLog2Ratio', ...
            'Log2ResidualDifference_ST1_75_Minus_ST1_68'});
        focusedPairRows = [focusedPairRows; newRow]; %#ok<AGROW>
    end
end
writetable(focusedPairRows, fullfile(tableDir, ...
    'extracellular_st175_st168_matched_pairs.csv'));

% Batch-paired ST1-6 contrasts across the 13 continuous metabolites.
focusedComparators = ["ST1-75", "ST1-68", "VPI10463", "R20291"];
st16Contrasts = table;
for comparator = focusedComparators
    rows6 = bacterialGroup == "ST1-6";
    rowsOther = bacterialGroup == comparator;
    [~, order6] = sort(bacterialBatch(rows6));
    [~, orderOther] = sort(bacterialBatch(rowsOther));
    values6 = matchedLog2Ratio(rows6, continuous);
    valuesOther = matchedLog2Ratio(rowsOther, continuous);
    values6 = values6(order6, :);
    valuesOther = valuesOther(orderOther, :);
    difference = values6 - valuesOther;
    pValue = nan(sum(continuous), 1);
    for j = 1:sum(continuous)
        [~, pValue(j)] = ttest(difference(:, j), 0);
    end
    qValue = bhAdjust(pValue);
    newRows = table(repmat("ST1-6", sum(continuous), 1), ...
        repmat(comparator, sum(continuous), 1), metabolite(continuous)', ...
        mean(difference, 1)', pValue, qValue, ...
        'VariableNames', {'StrainA', 'StrainB', 'Metabolite', ...
        'MeanMatchedLog2ResidualDifference_AminusB', 'PairedTP', ...
        'BenjaminiHochbergQ'});
    st16Contrasts = [st16Contrasts; newRows]; %#ok<AGROW>
end
writetable(st16Contrasts, fullfile(tableDir, ...
    'st16_extracellular_matched_contrasts.csv'));

%% Figures
makeFiveStrainFigure(groupMeanLog2Ratio(:, continuous), ...
    metabolite(continuous), strainOrder, strainColors, ...
    matchedLog2Ratio, bacterialGroup, bacterialBatch, metabolite, figureDir);
makeCrossAssayFigure(primaryMappedCall, mappedLog2Ratio, mapGcms, ...
    crossStrains, breadthTable, figureDir);

%% Human-readable summary
summaryFile = fullfile(tableDir, 'analysis_summary.txt');
fid = fopen(summaryFile, 'w');
assert(fid >= 0, 'Unable to write analysis summary.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Extracellular GC-MS and PM1 concordance analysis\n');
fprintf(fid, '================================================\n');
fprintf(fid, 'GC-MS source: %s\n', gcmsFile);
fprintf(fid, 'PM1 source: %s\n', pm1File);
fprintf(fid, 'PM1 breadth ST1-75/ST1-68/ST1-6/VPI10463: %s.\n', ...
    strjoin(string(pm1Breadth'), ', '));
fprintf(fid, ['PM1 amino-acid-class breadth in the same order: %s ' ...
    '(of %d annotated substrates).\n'], ...
    strjoin(string(pm1AminoAcidBreadth'), ', '), sum(aminoAcidRow));
fprintf(fid, 'GC-MS below-media counts in the same order: %s.\n', ...
    strjoin(string(gcmsDepletionBreadth'), ', '));
fprintf(fid, ['Direct overlap, primary mapping: mean within-substrate ' ...
    'depletion difference %.4f, permutation p %.4f (%d informative).\n'], ...
    primaryEffect, primaryP, primaryInformative);
fprintf(fid, ['Any-stereoisomer sensitivity: mean within-substrate ' ...
    'depletion difference %.4f, permutation p %.4f (%d informative).\n'], ...
    stereoEffect, stereoP, stereoInformative);
fprintf(fid, ['Whole-profile distance concordance: exact Mantel Spearman ' ...
    'rho %.4f, p %.4f (four shared strains).\n'], mantelRho, mantelP);
fprintf(fid, 'Global extracellular metabolites at BH q<0.05: %s.\n', ...
    strjoin(globalTests.Metabolite(globalTests.BenjaminiHochbergQ < 0.05), ', '));
prolineFocused = focusedTests(focusedTests.Metabolite == "Proline", :);
fprintf(fid, ['Focused ST1-75/ST1-68 proline geometric residual-signal ' ...
    'ratio %.4f; paired p %.5f; BH q %.5f across %d nonconstant ' ...
    'metabolites; exact sign-flip p %.4f.\n'], ...
    prolineFocused.GeometricResidualSignalRatio_ST1_75_Over_ST1_68, ...
    prolineFocused.PairedTP, ...
    prolineFocused.BenjaminiHochbergQAcrossNonconstantMetabolites, ...
    prolineFocused.NNonconstantMetabolitesInFdrFamily, ...
    prolineFocused.ExactTwoSidedSignFlipP);
fprintf(fid, ['R20291 is absent from the PM1 matrix and is excluded only ' ...
    'from direct cross-assay comparisons.\n']);
fprintf(fid, ['PM1 calls represent substrate-response capacity; GC-MS ' ...
    'represents net extracellular remodeling in one medium.\n']);
disp(fileread(summaryFile));

%% Local functions
function [effect, pValue, nInformative] = ...
        callDepletionConcordance(callMatrix, depletion, nPerm)
informative = false(size(callMatrix, 1), 1);
for m = 1:size(callMatrix, 1)
    informative(m) = numel(unique(callMatrix(m, :))) > 1 && ...
        std(depletion(m, :)) > 1e-12;
end
nInformative = sum(informative);
effect = meanWithinSubstrateEffect(callMatrix(informative, :), ...
    depletion(informative, :));
nullEffect = nan(nPerm, 1);
for p = 1:nPerm
    permuted = callMatrix(informative, :);
    for m = 1:size(permuted, 1)
        permuted(m, :) = permuted(m, randperm(size(permuted, 2)));
    end
    nullEffect(p) = meanWithinSubstrateEffect(permuted, ...
        depletion(informative, :));
end
pValue = (1 + sum(abs(nullEffect) >= abs(effect) - 1e-12)) / (nPerm + 1);
end

function value = meanWithinSubstrateEffect(callMatrix, depletion)
perSubstrate = nan(size(callMatrix, 1), 1);
for m = 1:size(callMatrix, 1)
    perSubstrate(m) = mean(depletion(m, callMatrix(m, :) == 1)) - ...
        mean(depletion(m, callMatrix(m, :) == 0));
end
value = mean(perSubstrate);
end

function [rho, pValue] = exactMantel(distanceA, distanceB)
upper = triu(true(size(distanceA)), 1);
observedA = distanceA(upper);
observedB = distanceB(upper);
rho = corr(observedA, observedB, 'Type', 'Spearman');
allPermutation = perms(1:size(distanceA, 1));
nullRho = nan(size(allPermutation, 1), 1);
for p = 1:size(allPermutation, 1)
    permuted = distanceB(allPermutation(p, :), allPermutation(p, :));
    nullRho(p) = corr(observedA, permuted(upper), 'Type', 'Spearman');
end
pValue = mean(abs(nullRho) >= abs(rho) - 1e-12);
end

function adjusted = bhAdjust(pValue)
adjusted = nan(size(pValue));
finiteIndex = find(isfinite(pValue));
[sortedP, order] = sort(pValue(finiteIndex));
m = numel(sortedP);
if m == 0
    return;
end
sortedQ = sortedP .* m ./ (1:m)';
for i = (m - 1):-1:1
    sortedQ(i) = min(sortedQ(i), sortedQ(i + 1));
end
sortedQ = min(sortedQ, 1);
restored = nan(m, 1);
restored(order) = sortedQ;
adjusted(finiteIndex) = restored;
end

function pValue = exactSignFlipP(difference)
n = numel(difference);
signPattern = dec2bin(0:(2 ^ n - 1)) - '0';
signPattern(signPattern == 0) = -1;
nullMean = abs(mean(signPattern .* difference', 2));
observed = abs(mean(difference));
pValue = mean(nullMean >= observed - 1e-12);
end

function [pseudoF, pValue, groupR2] = balancedGroupPermutationTest( ...
        values, group, batch, nPerm)
grand = mean(values);
totalSS = sum((values - grand) .^ 2);
uniqueGroup = unique(group, 'stable');
uniqueBatch = unique(batch, 'stable');
groupSS = 0;
for g = uniqueGroup'
    rows = group == g;
    groupSS = groupSS + sum(rows) * (mean(values(rows)) - grand) ^ 2;
end
batchSS = 0;
for b = uniqueBatch'
    rows = batch == b;
    batchSS = batchSS + sum(rows) * (mean(values(rows)) - grand) ^ 2;
end
residualSS = max(eps, totalSS - groupSS - batchSS);
groupDf = numel(uniqueGroup) - 1;
residualDf = numel(values) - numel(uniqueGroup) - numel(uniqueBatch) + 1;
pseudoF = (groupSS / groupDf) / (residualSS / residualDf);
groupR2 = groupSS / totalSS;
nullF = nan(nPerm, 1);
for p = 1:nPerm
    permutedGroup = group;
    for b = uniqueBatch'
        rows = find(batch == b);
        permutedGroup(rows) = group(rows(randperm(numel(rows))));
    end
    permutedSS = 0;
    for g = uniqueGroup'
        rows = permutedGroup == g;
        permutedSS = permutedSS + sum(rows) * ...
            (mean(values(rows)) - grand) ^ 2;
    end
    permutedResidual = max(eps, totalSS - permutedSS - batchSS);
    nullF(p) = (permutedSS / groupDf) / ...
        (permutedResidual / residualDf);
end
pValue = (1 + sum(nullF >= pseudoF - 1e-12)) / (nPerm + 1);
end

function makeFiveStrainFigure(groupMean, metabolite, strainOrder, colors, ...
        matchedLog2Ratio, bacterialGroup, bacterialBatch, allMetabolite, ...
        figureDir)
[~, order] = sort(var(groupMean, 0, 1), 'descend');
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1450, 760]);
layout = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(layout, [2, 2]);
colorLimits = [-9, 3];
imagesc(ax1, groupMean(:, order)', colorLimits);
colormap(ax1, zeroAnchoredBlueWhiteRed(colorLimits, 256));
cb = colorbar(ax1);
cb.Label.String = 'Mean matched log_2 residual signal';
set(ax1, 'XTick', 1:numel(strainOrder), 'XTickLabel', strainOrder, ...
    'XTickLabelRotation', 30, 'YTick', 1:numel(metabolite), ...
    'YTickLabel', metabolite(order), 'FontName', 'Arial', 'FontSize', 10);
title(ax1, 'A  Five-strain extracellular CoRe profiles', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');

plotFocusedResidual(nexttile(layout), "Proline", matchedLog2Ratio, ...
    bacterialGroup, bacterialBatch, allMetabolite, strainOrder, colors);
title('B  Proline', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
plotFocusedResidual(nexttile(layout), "Glutamate", matchedLog2Ratio, ...
    bacterialGroup, bacterialBatch, allMetabolite, strainOrder, colors);
title('C  Glutamate', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
sgtitle(layout, ['Protective strains can occupy different extracellular ' ...
    'metabolic programs'], 'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold');
exportgraphics(figureHandle, fullfile(figureDir, ...
    'extracellular_five_strain_programs.png'), 'Resolution', 300);
close(figureHandle);
end

function plotFocusedResidual(ax, target, values, group, batch, ...
        metabolite, strainOrder, colors)
hold(ax, 'on');
column = find(metabolite == target, 1);
jitter = [-0.10; 0; 0.10];
if target == "Proline"
    rows75 = find(group == "ST1-75");
    rows68 = find(group == "ST1-68");
    [~, order75] = sort(batch(rows75));
    [~, order68] = sort(batch(rows68));
    rows75 = rows75(order75);
    rows68 = rows68(order68);
    for i = 1:numel(rows75)
        plot(ax, [1 + jitter(i), 2 + jitter(i)], ...
            [values(rows75(i), column), values(rows68(i), column)], '-', ...
            'Color', [0.72, 0.74, 0.76], 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end
end
for g = 1:numel(strainOrder)
    rows = find(group == strainOrder(g));
    [~, order] = sort(batch(rows));
    rows = rows(order);
    x = g + jitter;
    scatter(ax, x, values(rows, column), 54, colors(g, :), 'filled', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
    plot(ax, [g - 0.16, g + 0.16], ...
        repmat(mean(values(rows, column)), 1, 2), '-', ...
        'Color', colors(g, :), 'LineWidth', 2.2);
end
yline(ax, 0, '--', 'Color', [0.55, 0.55, 0.55]);
set(ax, 'XTick', 1:numel(strainOrder), 'XTickLabel', strainOrder, ...
    'XTickLabelRotation', 35, 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax, 'Matched log_2 residual / media');
box(ax, 'off');
end

function makeCrossAssayFigure(callMatrix, log2Ratio, metabolite, ...
        strainOrder, breadth, figureDir)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1500, 650]);
layout = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(layout);
imagesc(ax1, callMatrix, [0, 1]);
colormap(ax1, [0.93, 0.93, 0.91; 0.16, 0.42, 0.66]);
set(ax1, 'XTick', 1:numel(strainOrder), 'XTickLabel', strainOrder, ...
    'XTickLabelRotation', 30, 'YTick', 1:numel(metabolite), ...
    'YTickLabel', metabolite, 'FontName', 'Arial', 'FontSize', 10);
title(ax1, 'A  Audited PM1 calls', 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold');
for m = 1:size(callMatrix, 1)
    for s = 1:size(callMatrix, 2)
        text(ax1, s, m, string(callMatrix(m, s)), ...
            'HorizontalAlignment', 'center', 'Color', ...
            callMatrix(m, s) * [1, 1, 1] + ...
            (1 - callMatrix(m, s)) * [0.25, 0.25, 0.25], ...
            'FontWeight', 'bold');
    end
end

ax2 = nexttile(layout);
colorLimits = [-9, 3];
imagesc(ax2, log2Ratio, colorLimits);
colormap(ax2, zeroAnchoredBlueWhiteRed(colorLimits, 256));
cb = colorbar(ax2);
cb.Label.String = 'Matched log_2 residual signal';
set(ax2, 'XTick', 1:numel(strainOrder), 'XTickLabel', strainOrder, ...
    'XTickLabelRotation', 30, 'YTick', 1:numel(metabolite), ...
    'YTickLabel', metabolite, 'FontName', 'Arial', 'FontSize', 10);
title(ax2, 'B  Extracellular CoRe outcome', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');

ax3 = nexttile(layout);
barHandle = bar(ax3, [breadth.Pm1PositiveCallFraction, ...
    breadth.GcmsBelowMediaRangeFraction], 'grouped');
barHandle(1).FaceColor = [0.16, 0.42, 0.66];
barHandle(2).FaceColor = [0.62, 0.66, 0.69];
ylim(ax3, [0, 1]);
set(ax3, 'XTick', 1:height(breadth), 'XTickLabel', breadth.Strain, ...
    'XTickLabelRotation', 30, 'FontName', 'Arial', 'FontSize', 10);
ylabel(ax3, 'Fraction of assay panel');
legend(ax3, {'PM1 positive calls (95)', 'GC-MS below media (18)'}, ...
    'Location', 'northoutside', 'Box', 'off');
title(ax3, 'C  Breadth depends on assay context', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
box(ax3, 'off');

sgtitle(layout, ['BIOLOG resource access and extracellular CoRe describe ' ...
    'complementary metabolic layers'], 'FontName', 'Arial', ...
    'FontSize', 16, 'FontWeight', 'bold');
exportgraphics(figureHandle, fullfile(figureDir, ...
    'gcms_pm1_cross_assay_comparison.png'), 'Resolution', 300);
close(figureHandle);
end

function map = zeroAnchoredBlueWhiteRed(limits, n)
blue = [0.08, 0.30, 0.58];
white = [0.98, 0.98, 0.97];
red = [0.70, 0.09, 0.12];
assert(limits(1) < 0 && limits(2) > 0, ...
    'Color limits must span zero.');
zeroIndex = 1 + round((0 - limits(1)) / diff(limits) * (n - 1));
first = [linspace(blue(1), white(1), zeroIndex)', ...
    linspace(blue(2), white(2), zeroIndex)', ...
    linspace(blue(3), white(3), zeroIndex)'];
secondLength = n - zeroIndex + 1;
second = [linspace(white(1), red(1), secondLength)', ...
    linspace(white(2), red(2), secondLength)', ...
    linspace(white(3), red(3), secondLength)'];
map = [first; second(2:end, :)];
end
