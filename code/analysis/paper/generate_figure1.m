%% Generate Figure 1: protection, enrichment, and quantification
% Combines observed weight trajectories, a composite event-free Kaplan-Meier
% display, relative qPCR composition, a worked display of the
% terminal-event-zero score convention, and model-derived protection effects.
% Death is the display label for the source death/humane-removal endpoint.
% The 1:5 estimate remains separate from the panel-wide
% mouse-screen model.

clear;

scriptDir = fileparts(mfilename('fullpath'));
localRoot = fileparts(fileparts(scriptDir));
repoRoot = fileparts(localRoot);
addpath(genpath(fullfile(localRoot, 'analysis')));

figureDir = fullfile(repoRoot, 'results', 'figures');
tableDir = fullfile(repoRoot, 'results', 'tables', 'main');
if ~isfolder(figureDir); mkdir(figureDir); end
if ~isfolder(tableDir); mkdir(tableDir); end

directFile = project_data_file('processed', 'mouse', 'scores', ...
    'ProtectionScreen_CDI_mouse.csv');
weightFile = project_data_file('raw', 'qpcr', 'weights.xlsx');
qpcrFile = fullfile(package_root(),'results','tables','main','qpcr_fractions.csv');
screenScoreFile = fullfile(repoRoot, 'results', 'tables', ...
    'mouse_weight_survival_reanalysis', ...
    'coinfection_protection_scores.csv');
survivalFile = fullfile(repoRoot, 'results', 'tables', ...
    'mouse_weight_survival_reanalysis', ...
    'coinfection_survival_comparisons.csv');
designSource = fullfile(package_root(),'README.md');
assert(all(isfile(string({directFile, weightFile, qpcrFile, ...
    screenScoreFile, survivalFile, designSource}))), ...
    'One or more Figure 1 sources are unavailable.');

%% Approximately equal-input challenge used in the screen
direct = readtable(directFile, 'TextType', 'string');
direct.experiment = string(direct.experiment);
direct.tx = string(direct.tx);
direct.cdiffstrain = string(direct.cdiffstrain);
direct.mouse = string(direct.mouse);
direct.animal_id = direct.experiment + "|" + direct.tx + "|" + ...
    direct.cdiffstrain + "|" + direct.mouse;
direct = direct(direct.experiment ~= "ks65", :);
assert(all(direct.death(ismissing(direct.relweight)) == 1), ...
    'Screen weights are missing outside terminal-event rows.');

sharedExperiments = intersect( ...
    unique(direct.experiment(direct.cdiffstrain == "vpi")), ...
    unique(direct.experiment(direct.cdiffstrain == "st1.75.vpi")));
directFocused = direct(ismember(direct.experiment, sharedExperiments) & ...
    ismember(direct.cdiffstrain, ["vpi", "st1.75.vpi"]), :);
assert(numel(sharedExperiments) == 7, ...
    'Expected seven contemporaneous ST1-75/VPI10463 experiments.');

screenScores = readtable(screenScoreFile, 'TextType', 'string');
screenSt175 = screenScores(screenScores.Strain == "ST1-75", :);
assert(height(screenSt175) == 1 && screenSt175.N_Animals == 21 && ...
    screenSt175.N_Experiments == 7, ...
    'The locked ST1-75 screen protection estimate changed unexpectedly.');
survival = readtable(survivalFile, 'TextType', 'string');
survival = survival(survival.Strain == "ST1-75", :);
assert(height(survival) == 1 && survival.Comparator_Events == 0 && ...
    survival.Reference_Events == 8, ...
    'The permanent ST1-75 survival comparison changed unexpectedly.');

%% Focused nominal 1:5 experiment
weights = readtable(weightFile, 'VariableNamingRule', 'preserve');
qpcr = readtable(qpcrFile, 'VariableNamingRule', 'preserve');
weights = renamevars(weights, "Var1", "Day");
weights = renamevars(weights, "ST1-75", "ST1_75");
qpcr = renamevars(qpcr, "ST1-75", "ST1_75");
weights = weights(~all(ismissing(weights), 2), :);
assert(isequal(weights.Day, [0; 2; 3]), ...
    'Unexpected focused weight time points.');

mixedQpcr = qpcr(ismember(qpcr.Mouse, [3, 4, 5]), :);
validPair = ~ismissing(mixedQpcr.ST1_75) & ~ismissing(mixedQpcr.VPI);
assert(all(abs(mixedQpcr.ST1_75(validPair) + ...
    mixedQpcr.VPI(validPair) - 100) < 1e-6), ...
    'Relative qPCR pairs do not sum to 100%.');

% Apply the same model form to the KS65 mixed versus VPI trajectories.
% This estimate is descriptive because the 1:5 cohort has only one VPI
% control animal. It remains separate from the panel-wide model.
focusedScore = buildFocusedScoreTable(weights);
focusedModel = fitlme(focusedScore, ...
    'Score ~ Condition + (1|DayGroup) + (1|Animal)', ...
    'FitMethod', 'REML');
focusedCoefficients = focusedModel.Coefficients;
focusedEffectRow = startsWith(string(focusedCoefficients.Name), ...
    "Condition_");
assert(sum(focusedEffectRow) == 1, ...
    'Expected one focused mixed-versus-VPI coefficient.');
focusedEffect = focusedCoefficients.Estimate(focusedEffectRow);
focusedLower = focusedCoefficients.Lower(focusedEffectRow);
focusedUpper = focusedCoefficients.Upper(focusedEffectRow);

Estimate_Type = ["Screen estimate"; "Focused descriptive estimate"];
Comparison = ["Approximately equal-input ST1-75+VPI10463 vs VPI10463"; ...
    "Nominal 1:5 mixture vs focused VPI10463 control"];
Effect = [screenSt175.Effect; focusedEffect];
CI_Lower = [screenSt175.CI_Lower; focusedLower];
CI_Upper = [screenSt175.CI_Upper; focusedUpper];
N_Mixed = [screenSt175.N_Animals; 3];
N_VPI = [24; 1];
Display_Inference = ["Estimate with 95% CI"; ...
    "Model-derived 95% CI shown; repeated weights with focused VPI n=1"];
effectTable = table(Estimate_Type, Comparison, Effect, CI_Lower, CI_Upper, ...
    N_Mixed, N_VPI, Display_Inference);
qpcrSummary = summarizeMixedQpcr(mixedQpcr);
weightLong = buildFocusedWeightLongTable(weights);
kaplanMeierDisplay = buildKaplanMeierDisplay(directFocused);
mixKm = kaplanMeierDisplay( ...
    kaplanMeierDisplay.Condition == "ST1-75 + VPI10463", :);
vpiKm = kaplanMeierDisplay( ...
    kaplanMeierDisplay.Condition == "VPI10463", :);
assert(sum(mixKm.N_Events) == 0 && sum(vpiKm.N_Events) == 8 && ...
    mixKm.Event_Free_Probability(end) == 1 && ...
    abs(vpiKm.Event_Free_Probability(end) - 16 / 24) < 1e-12, ...
    'The pooled Kaplan-Meier display changed unexpectedly.');
writetable(effectTable, fullfile(tableDir, ...
    'figure1_protection_effects.csv'));
writetable(directFocused(:, {'experiment', 'day', 'mouse', 'tx', ...
    'cdiffstrain', 'relweight', 'death', 'animal_id'}), ...
    fullfile(tableDir, 'figure1_equal_input_direct_challenge.csv'));
writetable(weightLong, fullfile(tableDir, ...
    'figure1_focused_weight_observations.csv'));
writetable(mixedQpcr, fullfile(tableDir, ...
    'figure1_focused_qpcr_observations.csv'));
writetable(qpcrSummary, fullfile(tableDir, ...
    'figure1_focused_qpcr_summary.csv'));
writetable(survival, fullfile(tableDir, ...
    'figure1_st175_survival_comparison.csv'));
writetable(kaplanMeierDisplay, fullfile(tableDir, ...
    'figure1_kaplan_meier_display.csv'));

%% Figure layout
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.2, 0.2, 16.5, 7.65]);

axA = axes(fig, 'Position', [0.045, 0.525, 0.32, 0.42]);
plotDirectProtection(axA, directFocused, sharedExperiments);
panelLabel(axA, 'A');

axB = axes(fig, 'Position', [0.405, 0.525, 0.18, 0.42]);
plotKaplanMeier(axB, kaplanMeierDisplay, survival);
panelLabel(axB, 'B');

axC = axes(fig, 'Position', [0.625, 0.525, 0.16, 0.42]);
plotFocusedWeights(axC, weights);
panelLabel(axC, 'C');

axD = axes(fig, 'Position', [0.825, 0.525, 0.15, 0.42]);
plotRelativeQpcr(axD, mixedQpcr);
panelLabel(axD, 'D');
identity = label_ks65_mice(axC, axD, fullfile(package_root(), ...
    'input', 'qpcr', 'Analysis', 'KS65.xlsx'));
writetable(identity, fullfile(tableDir, 'figure1_ks65_mouse_identity.csv'));

% Panel E: worked scoring convention and model definition.
annotation(fig, 'textbox', [0.013, 0.405, 0.59, 0.045], ...
    'String', 'E   Protection-score calculation', ...
    'EdgeColor', 'none', 'FontName', 'Arial', 'FontSize', 13, ...
    'FontWeight', 'bold', 'Margin', 0);
plotScoringWorkflow(fig, directFocused);

axF = axes(fig, 'Position', [0.75, 0.06, 0.22, 0.34]);
plotProtectionEffects(axF, screenSt175, focusedEffect, focusedLower, ...
    focusedUpper);
annotation(fig, 'textbox', [0.715, 0.405, 0.025, 0.045], ...
    'String', 'F', 'EdgeColor', 'none', 'FontName', 'Arial', ...
    'FontSize', 15, 'FontWeight', 'bold', 'Margin', 0);

outputPng = fullfile(figureDir, ...
    'figure1_st175_protection_and_competitive_enrichment.png');
export_vector_asset(fig, outputPng);
exportgraphics(fig, outputPng, 'Resolution', 300);
close(fig);



Input = ["Panel co-infection mouse table"; ...
    "Permanent protection-score model"; ...
    "Permanent stratified survival comparison"; ...
    "Focused relative-weight workbook"; ...
    "Focused relative-qPCR workbook"; ...
    "Approximate 1:5 design-target decision"; ...
    "Figure 1 MATLAB workflow"];
Path = string({directFile; screenScoreFile; survivalFile; weightFile; ...
    qpcrFile; designSource; [mfilename('fullpath'), '.m']});
Use = ["Approximately equal-input ST1-75+VPI10463 protection phenotype"; ...
    "Locked equal-input protection estimate and confidence interval"; ...
    "Kaplan-Meier display, terminal events, and stratified comparison"; ...
    "All available focused weight trajectories and descriptive model"; ...
    "Observed relative strain fractions in mixed-infection mice"; ...
    "Design reference only; not a measured baseline"; ...
    "Figure assembly and vector export"];
writetable(table(Input, Path, Use), fullfile(tableDir, ...
    'figure1_source_provenance.csv'));

summaryFile = fullfile(tableDir, 'figure1_analysis_summary.txt');
fid = fopen(summaryFile, 'w');
assert(fid >= 0, 'Unable to create Figure 1 summary.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Figure 1: ST1-75 protection, enrichment, and quantification\n');
fprintf(fid, ['Approximately equal-input comparison: 21 ST1-75+VPI10463 ' ...
    'mice and 24 contemporaneous VPI10463 mice across 7 experiments.\n']);
fprintf(fid, 'Terminal events: 0/21 ST1-75+VPI10463 and 8/24 VPI10463.\n');
fprintf(fid, 'Stratified log-rank p %.12g.\n', ...
    survival.P_Value);
fprintf(fid, ['Locked approximately equal-input protection effect: %.6f ' ...
    '(95%% CI %.6f to %.6f).\n'], screenSt175.Effect, ...
    screenSt175.CI_Lower, screenSt175.CI_Upper);
fprintf(fid, ['Focused nominal 1:5 descriptive effect: %.6f ' ...
    '(model 95%% CI %.6f to %.6f; repeated weights with VPI n=1).\n'], ...
    focusedEffect, focusedLower, focusedUpper);
for i = 1:height(qpcrSummary)
    fprintf(fid, ['Day %.3g: mean ST1-75 relative-qPCR fraction %.6f%%, ' ...
        'n=%d.\n'], qpcrSummary.Day(i), ...
        qpcrSummary.Mean_ST1_75_Fraction_Percent(i), ...
        qpcrSummary.N_Mice_With_Observed_Pair(i));
end
fprintf(fid, 'PNG: %s\n', outputPng);
disp(fileread(summaryFile));

function summary = summarizeMixedQpcr(qpcr)
days = sort(unique(qpcr.Day));
Day = days;
N_Mice_With_Observed_Pair = zeros(numel(days), 1);
Mean_ST1_75_Fraction_Percent = nan(numel(days), 1);
Minimum_ST1_75_Fraction_Percent = nan(numel(days), 1);
Maximum_ST1_75_Fraction_Percent = nan(numel(days), 1);
for i = 1:numel(days)
    values = qpcr.ST1_75(qpcr.Day == days(i));
    values = values(~ismissing(values));
    N_Mice_With_Observed_Pair(i) = numel(values);
    Mean_ST1_75_Fraction_Percent(i) = mean(values);
    Minimum_ST1_75_Fraction_Percent(i) = min(values);
    Maximum_ST1_75_Fraction_Percent(i) = max(values);
end
summary = table(Day, N_Mice_With_Observed_Pair, ...
    Mean_ST1_75_Fraction_Percent, Minimum_ST1_75_Fraction_Percent, ...
    Maximum_ST1_75_Fraction_Percent);
end

function long = buildFocusedWeightLongTable(weights)
conditions = ["ST1-75 alone", "VPI10463 alone", ...
    "Nominal 1:5 mixed trajectory 1", ...
    "Nominal 1:5 mixed trajectory 2", ...
    "Nominal 1:5 mixed trajectory 3"];
variables = ["ST1_75", "VPI", "mix", "mix_1", "mix_2"];
Day = repmat(weights.Day, numel(variables), 1);
Condition = strings(numel(Day), 1);
Relative_Weight_Percent = nan(numel(Day), 1);
row = 0;
for j = 1:numel(variables)
    values = weights.(variables(j));
    for i = 1:height(weights)
        row = row + 1;
        Condition(row) = conditions(j);
        Relative_Weight_Percent(row) = values(i);
    end
end
long = table(Day, Condition, Relative_Weight_Percent);
end

function focused = buildFocusedScoreTable(weights)
mixed = [weights.mix, weights.mix_1, weights.mix_2];
Day = [weights.Day; repmat(weights.Day, 3, 1)];
Condition = [repmat("VPI10463", height(weights), 1); ...
    repmat("Nominal 1:5 mixture", 3 * height(weights), 1)];
Animal = [repmat("VPI_control_1", height(weights), 1); ...
    repelem(["Mixed_trajectory_1"; "Mixed_trajectory_2"; ...
    "Mixed_trajectory_3"], height(weights))];
Score = [weights.VPI; mixed(:)];
DayGroup = categorical(string(Day));
Condition = categorical(Condition, ...
    ["VPI10463", "Nominal 1:5 mixture"]);
Animal = categorical(Animal);
focused = table(Score, Condition, DayGroup, Animal, Day);
end

function displayTable = buildKaplanMeierDisplay(raw)
conditions = ["st1.75.vpi", "vpi"];
labels = ["ST1-75 + VPI10463", "VPI10463"];
allRows = cell(numel(conditions), 1);
for g = 1:numel(conditions)
    group = raw(raw.cdiffstrain == conditions(g), :);
    ids = unique(group.animal_id);
    followUpDay = nan(numel(ids), 1);
    event = false(numel(ids), 1);
    for i = 1:numel(ids)
        one = sortrows(group(group.animal_id == ids(i), :), 'day');
        terminalRow = find(one.death == 1, 1, 'first');
        if isempty(terminalRow)
            followUpDay(i) = max(one.day);
        else
            followUpDay(i) = one.day(terminalRow);
            event(i) = true;
        end
    end

    days = sort(unique(group.day));
    nAtRisk = zeros(numel(days), 1);
    nEvents = zeros(numel(days), 1);
    probability = ones(numel(days), 1);
    currentProbability = 1;
    for d = 1:numel(days)
        nAtRisk(d) = sum(followUpDay >= days(d));
        nEvents(d) = sum(event & followUpDay == days(d));
        if nEvents(d) > 0
            currentProbability = currentProbability * ...
                (1 - nEvents(d) / nAtRisk(d));
        end
        probability(d) = currentProbability;
    end
    condition = repmat(labels(g), numel(days), 1);
    allRows{g} = table(condition, days, nAtRisk, nEvents, probability, ...
        'VariableNames', {'Condition', 'Day', 'N_At_Risk', 'N_Events', ...
        'Event_Free_Probability'});
end
displayTable = vertcat(allRows{:});
end

function plotDirectProtection(ax, raw, sharedExperiments)
conditions = ["vpi", "st1.75.vpi"];
colors = [0.31, 0.31, 0.31; 0.08, 0.39, 0.72];
labels = ["VPI10463", "ST1-75 + VPI10463 (1:1 co-infection)"];
days = sort(unique(raw.day));
handles = gobjects(2, 1);
hold(ax, 'on');
for g = 1:2
    group = raw(raw.cdiffstrain == conditions(g), :);
    ids = unique(group.animal_id);
    pale = 0.82 * [1, 1, 1] + 0.18 * colors(g, :);
    for i = 1:numel(ids)
        one = sortrows(group(group.animal_id == ids(i), :), 'day');
        observed = ~ismissing(one.relweight);
        plot(ax, one.day(observed), one.relweight(observed), '-', ...
            'Color', pale, 'LineWidth', 0.75, 'HandleVisibility', 'off');
        if any(one.death == 1)
            lastObserved = find(observed, 1, 'last');
            plot(ax, one.day(lastObserved), one.relweight(lastObserved), ...
                'x', 'Color', [0, 0, 0], 'MarkerSize', 7, ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end

    experimentMean = nan(numel(sharedExperiments), numel(days));
    for e = 1:numel(sharedExperiments)
        oneExperiment = group(group.experiment == sharedExperiments(e), :);
        for d = 1:numel(days)
            values = oneExperiment.relweight(oneExperiment.day == days(d));
            experimentMean(e, d) = mean(values, 'omitnan');
        end
    end
    displayedMean = mean(experimentMean, 1, 'omitnan');
    handles(g) = plot(ax, days, displayedMean, '-', ...
        'Color', colors(g, :), 'LineWidth', 2.8, ...
        'DisplayName', labels(g));
end

eventHandle = plot(ax, nan, nan, 'x', ...
    'Color', [0, 0, 0], 'MarkerSize', 7, 'LineWidth', 1.5, ...
    'DisplayName', 'Death');
yline(ax, 100, ':', 'Color', [0.55, 0.55, 0.55], ...
    'HandleVisibility', 'off');
xlabel(ax, 'Day after simultaneous challenge');
ylabel(ax, 'Body weight (% of baseline)');
title(ax, '1:1 co-infection dampens disease', 'FontWeight', 'bold');
xlim(ax, [min(days), max(days)]);
ylim(ax, [70, 116]);
legend(ax, [handles; eventHandle], 'Location', 'southwest', ...
    'Box', 'off', 'FontSize', 7.7);
formatAxes(ax);
end

function plotKaplanMeier(ax, displayTable, survival)
conditions = ["ST1-75 + VPI10463", "VPI10463"];
colors = [0.08, 0.39, 0.72; 0.31, 0.31, 0.31];
hold(ax, 'on');
for g = 1:numel(conditions)
    one = displayTable(displayTable.Condition == conditions(g), :);
    stairs(ax, one.Day, one.Event_Free_Probability, ...
        'Color', colors(g, :), 'LineWidth', 2.4, ...
        'HandleVisibility', 'off');
end
text(ax, 7.16, 1.00, 'Mix', 'Color', colors(1, :), ...
    'FontSize', 7.1, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
text(ax, 7.16, 0.67, 'VPI', 'Color', colors(2, :), ...
    'FontSize', 7.1, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');

riskDays = [0, 2, 4, 6, 7];
riskY = [-0.08, -0.18];
for g = 1:numel(conditions)
    one = displayTable(displayTable.Condition == conditions(g), :);
    for d = 1:numel(riskDays)
        row = one.Day == riskDays(d);
        assert(sum(row) == 1, 'A requested at-risk day is unavailable.');
        text(ax, riskDays(d), riskY(g), sprintf('%d', one.N_At_Risk(row)), ...
            'HorizontalAlignment', 'center', 'FontSize', 7.1, ...
            'Color', colors(g, :));
    end
end
text(ax, 3.5, -0.01, 'Number at risk', 'FontSize', 7.1, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(ax, -1.20, riskY(1), 'Mix', 'FontSize', 7.1, ...
    'Color', colors(1, :), 'HorizontalAlignment', 'left');
text(ax, -1.20, riskY(2), 'VPI', 'FontSize', 7.1, ...
    'Color', colors(2, :), 'HorizontalAlignment', 'left');
text(ax, 0.98, 0.46, sprintf(['Stratified log-rank\n' ...
    'p = %.3g'], survival.P_Value), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontSize', 7.0, ...
    'Color', [0.22, 0.22, 0.22]);
xlabel(ax, 'Day after simultaneous challenge');
ylabel(ax, 'Survival probability');
title(ax, 'Co-infection increases survival', 'FontWeight', 'bold');
xlim(ax, [-1.25, 7.75]);
ylim(ax, [-0.24, 1.05]);
xticks(ax, riskDays);
yticks(ax, [0, 0.25, 0.5, 0.75, 1]);
formatAxes(ax);
end

function plotFocusedWeights(ax, weights)
hold(ax, 'on');
blue = [0.08, 0.39, 0.72];
red = [0.78, 0.20, 0.16];
mixedColor = [0.12, 0.58, 0.58];
plot(ax, weights.Day, weights.ST1_75, '-o', 'Color', blue, ...
    'LineWidth', 2.0, 'MarkerFaceColor', blue, ...
    'DisplayName', 'ST1-75 alone (n=1)');
plot(ax, weights.Day, weights.VPI, '-o', 'Color', red, ...
    'LineWidth', 2.0, 'MarkerFaceColor', red, ...
    'DisplayName', 'VPI10463 alone (n=1)');
mixed = [weights.mix, weights.mix_1, weights.mix_2];
for i = 1:size(mixed, 2)
    visibility = 'off';
    displayName = '';
    if i == 1
        visibility = 'on';
        displayName = 'Nominal 1:5 mixture (n=3)';
    end
    plot(ax, weights.Day, mixed(:, i), '-o', 'Color', mixedColor, ...
        'LineWidth', 1.5, 'MarkerFaceColor', 'w', ...
        'DisplayName', displayName, 'HandleVisibility', visibility);
end
yline(ax, 100, ':', 'Color', [0.55, 0.55, 0.55], ...
    'HandleVisibility', 'off');
xlabel(ax, 'Day after inoculation');
ylabel(ax, 'Body weight (% of baseline)');
title(ax, '1:5 co-infection dampens disease', 'FontWeight', 'bold');
xlim(ax, [0, 3]);
ylim(ax, [70, 112]);
xticks(ax, [0, 1, 2, 3]);
legend(ax, 'Location', 'southwest', 'Box', 'off', 'FontSize', 6.5);
formatAxes(ax);
end

function plotRelativeQpcr(ax, qpcr)
hold(ax, 'on');
mice = [3, 4, 5];
colors = [0.08, 0.34, 0.70; 0.15, 0.55, 0.75; 0.38, 0.72, 0.78];
for i = 1:numel(mice)
    one = sortrows(qpcr(qpcr.Mouse == mice(i), :), 'Day');
    % Inoculum reference guide only; never append it to measured observations.
    first = find(~ismissing(one.ST1_75), 1);
    plot(ax, [0, one.Day(first)], [100/6, one.ST1_75(first)], '--', ...
        'Color', [0.55, 0.25, 0.10], 'LineWidth', 0.6, ...
        'HandleVisibility', 'off');
    plot(ax, one.Day, one.ST1_75, '-o', 'Color', colors(i, :), ...
        'LineWidth', 1.7, 'MarkerFaceColor', colors(i, :), ...
        'DisplayName', sprintf('Mouse %d', mice(i)));
end

plot(ax, 0, 100 / 6, 'd', 'Color', [0.55, 0.25, 0.10], ...
    'MarkerFaceColor', 'w', 'MarkerSize', 7, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
text(ax, 0.08, 100 / 6, '1:5 inoculum', ...
    'Color', [0.45, 0.20, 0.08], 'FontSize', 7.1, ...
    'VerticalAlignment', 'middle');
text(ax, 3.0, 101.5, 'All >99%', 'HorizontalAlignment', 'right', ...
    'FontSize', 8.2, 'FontWeight', 'bold', 'Color', [0.04, 0.22, 0.48]);
xlabel(ax, 'Day after inoculation');
ylabel(ax, 'ST1-75 fraction of qPCR signal (%)');
title(ax, 'ST1-75 increases in frequency', 'FontWeight', 'bold');
xlim(ax, [-0.05, 3.12]);
ylim(ax, [0, 106]);
xticks(ax, [0, 0.5, 1, 3]);
xticklabels(ax, {'0', '0.5', '1', '3'});
legend(ax, 'Location', 'southeast', 'Box', 'off', 'FontSize', 6.8);
formatAxes(ax);
end

function plotScoringWorkflow(fig, directFocused)
terminalIds = unique(directFocused.animal_id( ...
    directFocused.cdiffstrain == "vpi" & directFocused.death == 1));
assert(~isempty(terminalIds), 'No terminal trajectory available for example.');
one = sortrows(directFocused(directFocused.animal_id == terminalIds(1), :), ...
    'day');
observed = ~ismissing(one.relweight);

axExample = axes(fig, 'Position', [0.045, 0.07, 0.145, 0.29]);
hold(axExample, 'on');
plot(axExample, one.day(observed), one.relweight(observed), '-o', ...
    'Color', [0.32, 0.32, 0.32], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.32, 0.32, 0.32]);
lastObserved = find(observed, 1, 'last');
plot(axExample, one.day(lastObserved), one.relweight(lastObserved), 'x', ...
    'Color', [0.78, 0.12, 0.10], 'MarkerSize', 8, 'LineWidth', 1.6);
yline(axExample, 100, ':', 'Color', [0.65, 0.65, 0.65]);
xlabel(axExample, 'Day');
ylabel(axExample, 'Body weight (% of baseline)');
title(axExample, '1. Observed weights', 'FontWeight', 'bold', ...
    'FontSize', 9.5);
ylim(axExample, [70, 104]);
formatAxes(axExample);

annotation(fig, 'arrow', [0.215, 0.245], [0.30, 0.30], ...
    'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.2);
annotation(fig, 'textbox', [0.272, 0.32, 0.22, 0.06], ...
    'String', '2. Composite disease score', 'EdgeColor', 'none', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 11, 'Margin', 0);
annotation(fig, 'textbox', [0.272, 0.19, 0.21, 0.13], ...
    'String', { ...
    'Observed weights retained', ...
    'Terminal event without a weight:', ...
    'score = 0 (not a measured weight)'}, ...
    'EdgeColor', 'none', 'BackgroundColor', 'none', ...
    'FontName', 'Arial', 'FontSize', 11, 'Margin', 0, ...
    'FitBoxToText', 'off');

annotation(fig, 'arrow', [0.445, 0.475], [0.30, 0.30], ...
    'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.2);
annotation(fig, 'textbox', [0.50, 0.32, 0.20, 0.06], ...
    'String', '3. Estimate protection', 'EdgeColor', 'none', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 11, 'Margin', 0);
annotation(fig, 'textbox', [0.50, 0.17, 0.20, 0.15], ...
    'String', { ...
    'score ~ condition', ...
    '       + (1|day)', ...
    '       + (1|animal)', ...
    '', ...
    'Positive coefficient = protection'}, ...
    'EdgeColor', 'none', 'BackgroundColor', 'none', ...
    'FontName', 'Arial', 'FontSize', 11, 'Margin', 0, ...
    'FitBoxToText', 'off');
end

function plotProtectionEffects(ax, screenSt175, focusedEffect, ...
    focusedLower, focusedUpper)
hold(ax, 'on');
yline(ax, 0, ':', 'Color', [0.5, 0.5, 0.5]);
effects = [screenSt175.Effect, focusedEffect];
lowerErrors = [screenSt175.Effect - screenSt175.CI_Lower, ...
    focusedEffect - focusedLower];
upperErrors = [screenSt175.CI_Upper - screenSt175.Effect, ...
    focusedUpper - focusedEffect];
bars = bar(ax, [1, 2], effects, 0.62, 'FaceColor', 'flat', ...
    'EdgeColor', 'none');
bars.CData = [0.08, 0.39, 0.72; 0.12, 0.58, 0.58];
errorbar(ax, [1, 2], effects, lowerErrors, upperErrors, ...
    'k.', 'LineWidth', 1.5, 'CapSize', 9);
xticks(ax, [1, 2]);
xticklabels(ax, {'Equal input (21 vs 24)', '1:5 mixture (3 vs 1)'});
ylabel(ax, 'Protection effect vs VPI10463 (score points)');
ylim(ax, [0, 24]);
xlim(ax, [0.45, 2.55]);
text(ax, 0.98, 0.97, ...
    'Model estimates and 95% CIs', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontSize', 7.2, ...
    'Color', [0.30, 0.30, 0.30]);
title(ax, 'Protection estimates', ...
    'FontWeight', 'bold');
formatAxes(ax);
end

function panelLabel(ax, label)
text(ax, -0.10, 1.08, label, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 15, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'top');
end

function formatAxes(ax)
set(ax, 'FontName', 'Arial', 'FontSize', 9, 'Box', 'off', ...
    'TickDir', 'out', 'LineWidth', 0.8);
end
