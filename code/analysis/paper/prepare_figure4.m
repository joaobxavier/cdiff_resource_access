%% Prepare strain-aligned inputs for Figure 4
% Combine the recombination-filtered phylogeny, mouse disease/protection
% scores, PM1 calls, chemical classes and consumer-resource model outputs.

clear;

scriptDir = fileparts(mfilename('fullpath'));
localRoot = fileparts(fileparts(scriptDir));
repoRoot = fileparts(localRoot);
addpath(genpath(fullfile(localRoot, 'analysis')));

figureDir = fullfile(repoRoot, 'results', 'figures');
tableDir = fullfile(repoRoot, 'results', 'tables', 'main');
if ~isfolder(figureDir); mkdir(figureDir); end
if ~isfolder(tableDir); mkdir(tableDir); end
set(0, 'DefaultFigureVisible', 'off');

matrixFile = fullfile(package_root(),'results','tables','biolog','Biolog_growth_matrix.csv');
classFile = package_input('moas.xlsx');
scoreFile = fullfile(repoRoot, 'results', 'tables', ...
    'mouse_weight_survival_reanalysis', ...
    'primary_terminal_event_zero_scores.csv');
treeFile = fullfile(package_root(),'results','phylogeny','tree_tests','masked_new22.analysis.nwk');

requiredFiles = string({matrixFile, classFile, scoreFile, treeFile});
assert(all(isfile(requiredFiles)), ...
    'One or more authoritative Figure 4 inputs are unavailable.');

pm1 = readtable(matrixFile, 'VariableNamingRule', 'preserve', ...
    'TextType', 'string');
classes = readtable(classFile, 'VariableNamingRule', 'preserve', ...
    'TextType', 'string');
scores = readtable(scoreFile, 'TextType', 'string');
pm1 = pm1(pm1.Metabolites ~= "Negative Control", :);

assert(height(pm1) == 95 && width(pm1) == 23, ...
    'Expected 95 substrates, 21 ST1 strains, and VPI.');
allCalls = double(table2array(pm1(:, 2:end)));
assert(all(ismember(allCalls(:), [0; 1])) && sum(allCalls, 'all') == 700, ...
    'The authoritative PM1 matrix must contain 700 binary positive calls.');
assert(height(scores) == 21 && numel(unique(scores.Strain)) == 21, ...
    'Expected 21 unique permanent mouse-score rows.');

treeData = extract_top_tree(treeFile);
treeOrder = treeData.TipOrder;

[matrixStrains, st1Columns, vpiColumn] = normalize_matrix_strains(pm1);
assert(numel(st1Columns) == 21 && ~isempty(vpiColumn), ...
    'Expected exactly 21 ST1 PM1 columns and one VPI column.');
[matchedMatrix, treeColumns] = ismember(treeOrder, matrixStrains);
assert(all(matchedMatrix) && numel(unique(treeColumns)) == 21, ...
    'Tree tips and PM1 candidate columns do not match one-to-one.');
[matchedScores, scoreRows] = ismember(treeOrder, scores.Strain);
assert(all(matchedScores) && numel(unique(scoreRows)) == 21, ...
    'Tree tips and permanent mouse scores do not match one-to-one.');

st1Calls = double(table2array(pm1(:, 1 + treeColumns)));
vpiCalls = double(table2array(pm1(:, 1 + vpiColumn)));
assert(sum(st1Calls, 'all') == 677 && sum(vpiCalls) == 23, ...
    'Candidate and VPI call totals changed unexpectedly.');

classOrder = ["amino acid", "carbohydrate", "carboxylic acid", ...
    "fatty acid", "amine", "alcohol", "amide", "ester"]';
[substrateClass, metadataRows] = align_substrate_classes(pm1, classes);
[rowOrder, classStarts, classEnds, classCounts] = ...
    order_substrates_by_class(substrateClass, metadataRows, classOrder);
assert(isequal(classCounts, [16; 38; 32; 3; 2; 2; 1; 1]), ...
    'The authoritative PM1 chemical-class counts changed unexpectedly.');

[niche, lambdaGrid] = compute_niche_metrics(pm1, scores);
association = compute_association_stats(niche);
chemical = compute_chemical_class_stats(pm1, classes, niche);
visualTrend = compute_visual_trend(niche);
[phylogeneticSignal, focalDistance] = ...
    compute_phylogenetic_context(treeFile, niche);

nicheTree = align_rows(niche, treeOrder);

focusOrder = ["ST1-75", "ST1-6", "ST1-68"]';
focused = align_rows(niche, focusOrder);
assert(isequal([focused.Shared_With_VPI, focused.ST1_Private, ...
    focused.VPI_Private], [22, 60, 1; 21, 55, 2; 16, 8, 7]), ...
    'Focused shared/private PM1 counts changed unexpectedly.');
assert(max(abs(focused.Model_ST1_Persistence_Fraction - ...
    [0.998003992015968; 0.926147704590818; 0.293413173652695])) ...
    < 1e-12, 'Focused model persistence fractions changed unexpectedly.');

primary = association(association.Analysis == "primary_total_breadth", :);
bootstrap = association(association.Analysis == "bootstrap_total_breadth", :);
assert(abs(primary.Spearman_Rho - 0.4569387937135) < 1e-12 && ...
    abs(primary.P_Value - 0.0373039598709594) < 1e-12 && ...
    abs(bootstrap.CI_Lower - (-0.048405583389638)) < 1e-12 && ...
    abs(bootstrap.CI_Upper - 0.802375781985214) < 1e-12, ...
    'The verified panel-wide association changed unexpectedly.');

amino = chemical(chemical.Category == "amino acid", :);
assert(height(amino) == 1 && ...
    abs(amino.Spearman_Rho - 0.702371354675232) < 1e-12 && ...
    abs(amino.Protection_Partial_Rank_Rho - 0.645063265979801) < 1e-12, ...
    'The verified amino-acid-class analysis changed unexpectedly.');
assert(abs(visualTrend.Slope(1) - 0.121446956516204) < 1e-12 && ...
    abs(visualTrend.Intercept(1) - 4.35419487636448) < 1e-12, ...
    'The descriptive Theil-Sen visual guide changed unexpectedly.');


function [matrixStrains, st1Columns, vpiColumn] = ...
    normalize_matrix_strains(pm1)
rawNames = string(pm1.Properties.VariableNames(2:end));
matrixStrains = replace(rawNames, "_", "-");
matrixStrains(rawNames == "VPI") = "VPI10463";
st1Columns = find(startsWith(matrixStrains, "ST1-"));
vpiColumn = find(matrixStrains == "VPI10463", 1);
end

function [substrateClass, metadataRows] = ...
    align_substrate_classes(pm1, classes)
[matched, metadataRows] = ismember(pm1.Metabolites, classes.Chemical);
assert(all(matched) && numel(unique(metadataRows)) == height(pm1), ...
    'Every PM1 substrate must match one unique class annotation.');
substrateClass = regexprep(classes.MoA(metadataRows), '^C-Source, ', '');
assert(~any(ismissing(substrateClass) | substrateClass == ""), ...
    'A PM1 substrate lacks a chemical-class annotation.');
end

function [rowOrder, classStarts, classEnds, classCounts] = ...
    order_substrates_by_class(substrateClass, metadataRows, classOrder)
rowOrder = zeros(0, 1);
classStarts = zeros(numel(classOrder), 1);
classEnds = zeros(numel(classOrder), 1);
classCounts = zeros(numel(classOrder), 1);
for i = 1:numel(classOrder)
    rows = find(substrateClass == classOrder(i));
    [~, orderWithinClass] = sort(metadataRows(rows), 'ascend');
    rows = rows(orderWithinClass);
    classStarts(i) = numel(rowOrder) + 1;
    rowOrder = [rowOrder; rows]; %#ok<AGROW>
    classEnds(i) = numel(rowOrder);
    classCounts(i) = numel(rows);
end
assert(numel(rowOrder) == 95 && numel(unique(rowOrder)) == 95, ...
    'Chemical-class ordering must retain every PM1 substrate exactly once.');
end

function ordered = align_rows(tableIn, strainOrder)
[matched, row] = ismember(strainOrder, tableIn.Strain);
assert(all(matched) && numel(unique(row)) == numel(strainOrder), ...
    'Unable to align one or more strain-resolved table rows.');
ordered = tableIn(row, :);
end

function treeData = extract_top_tree(treeFile)
treeData = tree_coordinates(treeFile);
end

function [niche, lambdaGrid] = compute_niche_metrics(pm1, scores)
rawNames = string(pm1.Properties.VariableNames(2:end));
vpiColumn = find(rawNames == "VPI", 1);
st1Columns = find(startsWith(rawNames, "ST1_"));
assert(~isempty(vpiColumn) && numel(st1Columns) == 21, ...
    'Expected VPI and 21 ST1 columns in the PM1 matrix.');

data = double(table2array(pm1(:, 2:end)));
vpi = data(:, vpiColumn) == 1;
lambdaGrid = linspace(0, max(sum(data, 1)), 501);
n = numel(st1Columns);

Strain = replace(rawNames(st1Columns)', "_", "-");
Total_Breadth = zeros(n, 1);
Shared_With_VPI = zeros(n, 1);
ST1_Private = zeros(n, 1);
VPI_Private = zeros(n, 1);
Private_Advantage = zeros(n, 1);
Model_ST1_Persistence_Fraction = zeros(n, 1);
Model_ST1_Exclusion_Fraction = zeros(n, 1);
Model_Coexistence_Fraction = zeros(n, 1);
Model_VPI_Exclusion_Fraction = zeros(n, 1);

for i = 1:n
    st1 = data(:, st1Columns(i)) == 1;
    shared = sum(st1 & vpi);
    st1Private = sum(st1 & ~vpi);
    vpiPrivate = sum(~st1 & vpi);
    Total_Breadth(i) = sum(st1);
    Shared_With_VPI(i) = shared;
    ST1_Private(i) = st1Private;
    VPI_Private(i) = vpiPrivate;
    Private_Advantage(i) = st1Private - vpiPrivate;
    [Model_ST1_Persistence_Fraction(i), ...
        Model_ST1_Exclusion_Fraction(i), ...
        Model_Coexistence_Fraction(i), ...
        Model_VPI_Exclusion_Fraction(i)] = ...
        model_fractions(shared, st1Private, vpiPrivate, lambdaGrid);
end

niche = table(Strain, Total_Breadth, Shared_With_VPI, ST1_Private, ...
    VPI_Private, Private_Advantage, Model_ST1_Persistence_Fraction, ...
    Model_ST1_Exclusion_Fraction, Model_Coexistence_Fraction, ...
    Model_VPI_Exclusion_Fraction);
[matched, scoreRows] = ismember(niche.Strain, scores.Strain);
assert(all(matched), 'A PM1 ST1 strain is missing from the score table.');
niche.Protection_Effect = scores.Protection_Effect(scoreRows);
niche.Protection_CI_Lower = scores.Protection_CI_Lower(scoreRows);
niche.Protection_CI_Upper = scores.Protection_CI_Upper(scoreRows);
niche.Mono_Disease_Effect = scores.Mono_Disease_Effect(scoreRows);
niche.Mono_Disease_CI_Lower = scores.Mono_Disease_CI_Lower(scoreRows);
niche.Mono_Disease_CI_Upper = scores.Mono_Disease_CI_Upper(scoreRows);

% Preserve the verified protection-descending order for the seeded bootstrap.
niche = sortrows(niche, 'Protection_Effect', 'descend');
end

function [persistence, st1Exclusion, coexistence, vpiExclusion] = ...
    model_fractions(shared, st1Private, vpiPrivate, lambdaGrid)
outcome = strings(size(lambdaGrid));
for k = 1:numel(lambdaGrid)
    lambda = lambdaGrid(k);
    st1Feasible = (shared + st1Private) > lambda;
    vpiFeasible = (shared + vpiPrivate) > lambda;
    if st1Feasible && ~vpiFeasible
        outcome(k) = "ST1 excludes VPI";
    elseif vpiFeasible && ~st1Feasible
        outcome(k) = "VPI excludes ST1";
    elseif ~st1Feasible && ~vpiFeasible
        outcome(k) = "neither feasible";
    else
        st1Invasion = st1Private - ...
            lambda * vpiPrivate / (shared + vpiPrivate);
        vpiInvasion = vpiPrivate - ...
            lambda * st1Private / (shared + st1Private);
        if st1Invasion > 0 && vpiInvasion > 0
            outcome(k) = "coexistence";
        elseif st1Invasion > 0 && vpiInvasion <= 0
            outcome(k) = "ST1 excludes VPI";
        elseif st1Invasion <= 0 && vpiInvasion > 0
            outcome(k) = "VPI excludes ST1";
        else
            outcome(k) = "priority";
        end
    end
end
st1Exclusion = mean(outcome == "ST1 excludes VPI");
coexistence = mean(outcome == "coexistence");
vpiExclusion = mean(outcome == "VPI excludes ST1");
persistence = st1Exclusion + coexistence;
end

function association = compute_association_stats(niche)
x = niche.Total_Breadth;
y = niche.Protection_Effect;
disease = niche.Mono_Disease_Effect;
[rhoTotal, pTotal] = corr(x, y, 'Type', 'Spearman');

rng(1);
bootstrapRho = nan(5000, 1);
for b = 1:numel(bootstrapRho)
    rows = randi(height(niche), height(niche), 1);
    if std(x(rows)) > 0 && std(y(rows)) > 0
        bootstrapRho(b) = corr(x(rows), y(rows), 'Type', 'Spearman');
    end
end
bootstrapCI = prctile(bootstrapRho(~isnan(bootstrapRho)), [2.5, 97.5]);

[rhoPrivate, pPrivate] = corr(niche.ST1_Private, y, 'Type', 'Spearman');
[rhoPersistence, pPersistence] = corr( ...
    niche.Model_ST1_Persistence_Fraction, y, 'Type', 'Spearman');
xRank = tiedrank(x);
yRank = tiedrank(y);
diseaseRank = tiedrank(disease);
covariates = [ones(height(niche), 1), diseaseRank];
xResidual = xRank - covariates * (covariates \ xRank);
yResidual = yRank - covariates * (covariates \ yRank);
[rhoPartial, pPartial] = corr(xResidual, yResidual, 'Type', 'Pearson');

Analysis = ["primary_total_breadth"; "bootstrap_total_breadth"; ...
    "st1_private_breadth"; "model_persistence_fraction"; ...
    "partial_rank_adjusted_for_mono_disease"];
N_Strains = repmat(height(niche), numel(Analysis), 1);
Spearman_Rho = [rhoTotal; rhoTotal; rhoPrivate; rhoPersistence; rhoPartial];
P_Value = [pTotal; pTotal; pPrivate; pPersistence; pPartial];
CI_Lower = [NaN; bootstrapCI(1); NaN; NaN; NaN];
CI_Upper = [NaN; bootstrapCI(2); NaN; NaN; NaN];
Notes = ["Permanent primary protection score"; ...
    "Percentile bootstrap, 5000 resamples, seed 1"; ...
    "ST1 assay-private PM1 calls"; ...
    "Derived from the same binary PM1 calls"; ...
    "Pearson correlation of rank residuals after disease-score rank"];
association = table(Analysis, N_Strains, Spearman_Rho, P_Value, ...
    CI_Lower, CI_Upper, Notes);
end

function visualTrend = compute_visual_trend(niche)
% A descriptive Theil-Sen line guides the eye without becoming the primary
% inferential model. The ribbon resamples strains and therefore represents
% uncertainty in the visual trend across the measured panel.
x = double(niche.Total_Breadth);
y = double(niche.Protection_Effect);
xGrid = linspace(min(x), max(x), 121)';
[slope, intercept] = theil_sen_fit(x, y);
fitValue = intercept + slope * xGrid;

nBootstrap = 5000;
rng(4, 'twister');
bootstrapFits = nan(nBootstrap, numel(xGrid));
for b = 1:nBootstrap
    rows = randi(numel(x), numel(x), 1);
    [bootstrapSlope, bootstrapIntercept] = ...
        theil_sen_fit(x(rows), y(rows));
    bootstrapFits(b, :) = bootstrapIntercept + ...
        bootstrapSlope * xGrid';
end
assert(all(isfinite(bootstrapFits), 'all'), ...
    'The visual-trend bootstrap produced a non-finite estimate.');
ciLower = prctile(bootstrapFits, 2.5, 1)';
ciUpper = prctile(bootstrapFits, 97.5, 1)';
Method = repmat("Theil-Sen visual guide; 5000 strain bootstraps; seed 4", ...
    numel(xGrid), 1);
N_Strains = repmat(numel(x), numel(xGrid), 1);
Slope = repmat(slope, numel(xGrid), 1);
Intercept = repmat(intercept, numel(xGrid), 1);
visualTrend = table(xGrid, fitValue, ciLower, ciUpper, Slope, ...
    Intercept, N_Strains, Method, 'VariableNames', ...
    {'Total_Breadth', 'Trend_Fit', 'CI_Lower', 'CI_Upper', ...
    'Slope', 'Intercept', 'N_Strains', 'Method'});
end

function [slope, intercept] = theil_sen_fit(x, y)
x = x(:);
y = y(:);
slopes = zeros(numel(x) * (numel(x) - 1) / 2, 1);
nSlopes = 0;
for i = 1:(numel(x) - 1)
    deltaX = x((i + 1):end) - x(i);
    deltaY = y((i + 1):end) - y(i);
    valid = deltaX ~= 0;
    nValid = sum(valid);
    if nValid > 0
        slopes((nSlopes + 1):(nSlopes + nValid)) = ...
            deltaY(valid) ./ deltaX(valid);
        nSlopes = nSlopes + nValid;
    end
end
assert(nSlopes > 0, ...
    'The Theil-Sen visual trend requires at least two distinct breadths.');
slope = median(slopes(1:nSlopes));
intercept = median(y - slope * x);
end

function chemical = compute_chemical_class_stats(pm1, classes, niche)
classes.Category = regexprep(classes.MoA, '^C-Source, ', '');
[matched, metadataRows] = ismember(pm1.Metabolites, classes.Chemical);
assert(all(matched), 'Unable to align PM1 substrates and class metadata.');
categoryBySubstrate = classes.Category(metadataRows);

rawNames = string(pm1.Properties.VariableNames(2:end));
st1Columns = find(startsWith(rawNames, "ST1_"));
st1Names = replace(rawNames(st1Columns)', "_", "-");
[matchedStrains, order] = ismember(niche.Strain, st1Names);
assert(all(matchedStrains), 'Unable to align PM1 and mouse strains.');

categories = unique(categoryBySubstrate);
n = numel(categories);
Category = categories;
Total_Substrates = zeros(n, 1);
Spearman_Rho = nan(n, 1);
P_Value = nan(n, 1);
Protection_Partial_Rank_Rho = nan(n, 1);
Protection_Partial_Rank_P_Value = nan(n, 1);
N_Strains = repmat(height(niche), n, 1);

for j = 1:n
    rows = categoryBySubstrate == categories(j);
    counts = sum(double(table2array( ...
        pm1(rows, 1 + st1Columns(order)))), 1)';
    Total_Substrates(j) = sum(rows);
    if std(counts) > 0
        [Spearman_Rho(j), P_Value(j)] = corr(counts, ...
            niche.Protection_Effect, 'Type', 'Spearman');
        countRank = tiedrank(counts);
        protectionRank = tiedrank(niche.Protection_Effect);
        breadthRank = tiedrank(niche.Total_Breadth);
        covariates = [ones(height(niche), 1), breadthRank];
        countResidual = countRank - covariates * ...
            (covariates \ countRank);
        protectionResidual = protectionRank - covariates * ...
            (covariates \ protectionRank);
        if std(countResidual) > 0 && std(protectionResidual) > 0
            [Protection_Partial_Rank_Rho(j), ...
                Protection_Partial_Rank_P_Value(j)] = corr( ...
                countResidual, protectionResidual, 'Type', 'Pearson');
        end
    end
end
Q_Value = bh_fdr(P_Value);
Protection_Partial_Rank_Q_Value = ...
    bh_fdr(Protection_Partial_Rank_P_Value);
chemical = table(Category, Total_Substrates, Spearman_Rho, P_Value, ...
    Q_Value, Protection_Partial_Rank_Rho, ...
    Protection_Partial_Rank_P_Value, ...
    Protection_Partial_Rank_Q_Value, N_Strains);
chemical = sortrows(chemical, 'Spearman_Rho', 'descend', ...
    'MissingPlacement', 'last');
end

function q = bh_fdr(p)
q = nan(size(p));
valid = ~isnan(p);
pv = p(valid);
[sorted, order] = sort(pv, 'ascend');
m = numel(sorted);
adjusted = sorted .* m ./ (1:m)';
adjusted = flipud(cummin(flipud(adjusted)));
adjusted = min(adjusted, 1);
restored = nan(m, 1);
restored(order) = adjusted;
q(valid) = restored;
end

function [signal, focalDistance] = ...
    compute_phylogenetic_context(treeFile, niche)
tree = phytreeread(treeFile);
leafNames = string(get(tree, 'LeafNames'));
distanceMatrix = squareform(pdist(tree, 'Nodes', 'leaves'));
isST1 = startsWith(leafNames, "ST1.");
leafNames = upper(replace(leafNames(isST1), ".", "-"));
distanceMatrix = distanceMatrix(isST1, isST1);
[leafNames, canonicalOrder] = sort(leafNames);
distanceMatrix = distanceMatrix(canonicalOrder, canonicalOrder);
[matched, nicheRows] = ismember(leafNames, niche.Strain);
assert(all(matched) && numel(leafNames) == 21, ...
    'The tree does not match all 21 candidate strains.');

breadth = niche.Total_Breadth(nicheRows);
protection = niche.Protection_Effect(nicheRows);
disease = niche.Mono_Disease_Effect(nicheRows);
traits = {breadth, protection, disease};
Trait = ["Total PM1 call breadth"; "Co-infection protection effect"; ...
    "Mono-colonization disease effect"];
upperTriangle = triu(true(size(distanceMatrix)), 1);
phylogeneticDistance = distanceMatrix(upperTriangle);

Mantel_Rho = nan(3, 1);
Permutation_P_Value = nan(3, 1);
N_Strains = repmat(numel(leafNames), 3, 1);
N_Permutations = repmat(10000, 3, 1);
Median_Patristic_Distance = repmat(median(phylogeneticDistance), 3, 1);
Maximum_Patristic_Distance = repmat(max(phylogeneticDistance), 3, 1);
rng(1);
for i = 1:numel(traits)
    values = traits{i};
    traitDifference = abs(values - values');
    Mantel_Rho(i) = corr(phylogeneticDistance, ...
        traitDifference(upperTriangle), 'Type', 'Spearman', ...
        'Rows', 'complete');
    nullRho = nan(N_Permutations(i), 1);
    for j = 1:N_Permutations(i)
        permuted = values(randperm(numel(values)));
        permutedDifference = abs(permuted - permuted');
        nullRho(j) = corr(phylogeneticDistance, ...
            permutedDifference(upperTriangle), 'Type', 'Spearman', ...
            'Rows', 'complete');
    end
    Permutation_P_Value(i) = (1 + sum(abs(nullRho) >= ...
        abs(Mantel_Rho(i)))) / (N_Permutations(i) + 1);
end
Notes = repmat("Two-sided Mantel-style label-permutation test; seed 1", ...
    3, 1);
signal = table(Trait, Mantel_Rho, Permutation_P_Value, N_Strains, ...
    N_Permutations, Median_Patristic_Distance, ...
    Maximum_Patristic_Distance, Notes);

idx75 = find(leafNames == "ST1-75", 1);
idx68 = find(leafNames == "ST1-68", 1);
Strain_1 = "ST1-75";
Strain_2 = "ST1-68";
Patristic_Distance_Substitutions_Per_Site = ...
    distanceMatrix(idx75, idx68);
focalDistance = table(Strain_1, Strain_2, ...
    Patristic_Distance_Substitutions_Per_Site);

end
