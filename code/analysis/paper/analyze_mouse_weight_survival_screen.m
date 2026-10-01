%% Reanalyze the mouse discovery screen using disease-course scores and survival
% This analysis deliberately excludes diarrhea because mono-colonization
% diarrhea was not recorded for ST1-68. The primary disease-course
% score uses the terminal-event-zero convention: a missing relative
% weight on a row recording death or humane-endpoint removal is assigned zero
% in the score model. That zero is an explicit composite-endpoint score value,
% not a measured body weight. Raw trajectories display observed weights only,
% with terminal events marked at the final observed weight. Survival is
% analyzed separately. The source field named `death` includes both observed
% deaths and humane-endpoint removals.

script_dir = fileparts(mfilename('fullpath'));
base_dir = fileparts(fileparts(script_dir));
project_dir = fileparts(base_dir);
addpath(genpath(fullfile(base_dir, 'analysis')));

figure_dir = fullfile(project_dir, 'results', 'diagnostic_figures', ...
    'mouse_weight_survival_reanalysis');
table_dir = fullfile(project_dir, 'results', 'tables', ...
    'mouse_weight_survival_reanalysis');
if ~exist(figure_dir, 'dir'); mkdir(figure_dir); end
if ~exist(table_dir, 'dir'); mkdir(table_dir); end

set(0, 'DefaultFigureVisible', 'off');

virulence_file = project_data_file('processed', 'mouse', 'scores', ...
    'Virulence_screen_clean_table.csv');
protection_file = project_data_file('processed', 'mouse', 'scores', ...
    'ProtectionScreen_CDI_mouse.csv');
virulence = readtable(virulence_file, 'TextType', 'string');
protection = readtable(protection_file, 'TextType', 'string');

required_variables = ["experiment", "day", "mouse", "tx", ...
    "cdiffstrain", "relweight", "death"];
validate_source_table(virulence, required_variables, "virulence");
validate_source_table(protection, required_variables, "protection");

virulence = add_animal_id(virulence);
protection = add_animal_id(protection);
validate_terminal_records(virulence, "virulence");
validate_terminal_records(protection, "protection");

% KS65 is the separate 1:5 inoculation experiment, not part of the
% cross-strain discovery screen.
protection_screen = protection(protection.experiment ~= "ks65", :);
all_screen = [virulence; protection_screen];

mono = all_screen((startsWith(all_screen.cdiffstrain, "st1.") | ...
    all_screen.cdiffstrain == "ui") & ...
    ~contains(all_screen.cdiffstrain, ".vpi") & ...
    ~contains(all_screen.cdiffstrain, ".st1.75"), :);
coinfection = protection_screen((startsWith(protection_screen.cdiffstrain, ...
    "st1.") & endsWith(protection_screen.cdiffstrain, ".vpi")) | ...
    protection_screen.cdiffstrain == "vpi", :);

% Preserve the condition universe used by the original score models. The
% final manuscript correlation is restricted below to the 21 matched ST1
% strains, but non-ST1 conditions remain in model fitting exactly as before.
score_primary = all_screen(~contains(all_screen.cdiffstrain, ".vpi") & ...
    ~contains(all_screen.cdiffstrain, ".st1.75"), :);
score_secondary = all_screen(endsWith(all_screen.cdiffstrain, ".vpi") | ...
    all_screen.cdiffstrain == "vpi", :);

[mono_weight, mono_model] = fit_terminal_event_zero_model(score_primary, ...
    "ui", "mono");
[protection_weight, protection_model] = fit_terminal_event_zero_model(score_secondary, ...
    "vpi", "protection");
mono_survival = compute_survival_comparisons(mono, "ui", "mono");
protection_survival = compute_survival_comparisons(coinfection, "vpi", ...
    "protection");

mono_join = mono_weight(:, {'Strain', 'Effect', 'CI_Lower', 'CI_Upper', ...
    'P_Value', 'Q_Value', 'N_Animals', 'N_Experiments', ...
    'N_Death_or_Humane_Endpoint_Events'});
mono_join = renamevars(mono_join, mono_join.Properties.VariableNames(2:end), ...
    {'Mono_Disease_Effect', 'Mono_Disease_CI_Lower', ...
    'Mono_Disease_CI_Upper', 'Mono_P_Value', 'Mono_Q_Value', ...
    'Mono_N_Animals', 'Mono_N_Experiments', ...
    'Mono_N_Death_or_Humane_Endpoint_Events'});
protection_join = protection_weight(:, {'Strain', 'Effect', 'CI_Lower', ...
    'CI_Upper', 'P_Value', 'Q_Value', 'N_Animals', 'N_Experiments', ...
    'N_Death_or_Humane_Endpoint_Events'});
protection_join = renamevars(protection_join, ...
    protection_join.Properties.VariableNames(2:end), ...
    {'Protection_Effect', 'Protection_CI_Lower', ...
    'Protection_CI_Upper', 'Protection_P_Value', 'Protection_Q_Value', ...
    'Protection_N_Animals', 'Protection_N_Experiments', ...
    'Protection_N_Death_or_Humane_Endpoint_Events'});
weight_joined = innerjoin(mono_join, protection_join, 'Keys', 'Strain');
weight_joined = weight_joined(startsWith(weight_joined.Strain, "ST1-"), :);

[weight_rho, weight_p] = corr(weight_joined.Mono_Disease_Effect, ...
    weight_joined.Protection_Effect, 'Type', 'Spearman', 'Rows', 'complete');

source_audit = readtable(fullfile(package_root(),'input_manifest.tsv'),'FileType','text','Delimiter','\t');
analysis_provenance = build_analysis_provenance(weight_rho, weight_p, ...
    height(weight_joined));

writetable(source_audit, fullfile(table_dir, 'source_audit.csv'));
writetable(mono_weight, fullfile(table_dir, 'mono_weight_effects.csv'));
writetable(mono_weight, fullfile(table_dir, ...
    'mono_disease_course_scores.csv'));
writetable(protection_weight, fullfile(table_dir, ...
    'coinfection_weight_effects.csv'));
writetable(protection_weight, fullfile(table_dir, ...
    'coinfection_protection_scores.csv'));
writetable(weight_joined, fullfile(table_dir, ...
    'weight_effects_joined.csv'));
writetable(weight_joined, fullfile(table_dir, ...
    'primary_terminal_event_zero_scores.csv'));
writetable(mono_survival, fullfile(table_dir, ...
    'mono_survival_comparisons.csv'));
writetable(protection_survival, fullfile(table_dir, ...
    'coinfection_survival_comparisons.csv'));
writetable(analysis_provenance, fullfile(table_dir, ...
    'analysis_provenance.csv'));

write_summary(fullfile(table_dir, 'analysis_summary.txt'), weight_joined, ...
    mono_survival, protection_survival, weight_rho, weight_p, ...
    mono_model, protection_model);
make_diagnostic_figure(figure_dir, weight_joined, coinfection, ...
    weight_rho, weight_p);

fprintf('Wrote disease-course score and survival tables to %s\n', table_dir);
fprintf('Wrote diagnostic figure to %s\n', figure_dir);

%% Local functions
function validate_source_table(tbl, required_variables, label)
variables = string(tbl.Properties.VariableNames);
assert(all(ismember(required_variables, variables)), ...
    '%s table is missing required variables.', label);
assert(all(ismember(unique(tbl.death(~ismissing(tbl.death))), [0; 1])), ...
    '%s death field contains values other than 0, 1, or missing.', label);
assert(~any(ismissing(tbl.death)), ...
    '%s death field contains missing values.', label);
end

function tbl = add_animal_id(tbl)
tbl.experiment = string(tbl.experiment);
tbl.tx = string(tbl.tx);
tbl.cdiffstrain = string(tbl.cdiffstrain);
tbl.mouse = string(tbl.mouse);
tbl.animal_id = tbl.experiment + "|" + tbl.tx + "|" + ...
    tbl.cdiffstrain + "|" + tbl.mouse;
animal_day = tbl.animal_id + "|" + string(tbl.day);
assert(numel(unique(animal_day)) == height(tbl), ...
    'Animal identity does not uniquely identify animal-day rows.');
end

function validate_terminal_records(tbl, label)
missing_weight = ismissing(tbl.relweight);
assert(all(tbl.death(missing_weight) == 1), ...
    '%s has missing weights outside terminal-event rows.', label);
assert(all(missing_weight(tbl.death == 1)), ...
    '%s has terminal-event rows with recorded weights; audit required.', label);
death_rows = find(tbl.death == 1);
for i = reshape(death_rows, 1, [])
    later_record = tbl.animal_id == tbl.animal_id(i) & tbl.day > tbl.day(i);
    assert(~any(later_record), ...
        '%s has observations after a recorded terminal event.', label);
end
end

function [effect_tbl, model] = fit_terminal_event_zero_model(raw, reference, mode)
% Fit the primary disease-course score. Day 0 is retained.
% Missing relative weights are permitted only on audited terminal-event rows
% and are assigned zero for this composite score. The score zero is not a
% measured body weight and is never drawn on the raw trajectory panels.
score_data = raw;
score_data.relweight(ismissing(score_data.relweight)) = 0;
score_data.exp_id = categorical(score_data.experiment + "_" + ...
    score_data.cdiffstrain + "_" + score_data.mouse);
score_data.cdiffstrain = categorical(score_data.cdiffstrain);
score_data.cdiffstrain = reordercats(score_data.cdiffstrain, ...
    [reference; setdiff(categories(score_data.cdiffstrain), reference)]);
if mode == "mono"
    score_data.relweight = -score_data.relweight;
end

model = fitlme(score_data, ...
    'relweight ~ cdiffstrain + (1|day) + (1|exp_id)');
coefficients = model.Coefficients;
coefficient_names = string(coefficients.Name);
keep = startsWith(coefficient_names, "cdiffstrain_");
condition = erase(coefficient_names(keep), "cdiffstrain_");
estimate = coefficients.Estimate(keep);
lower = coefficients.Lower(keep);
upper = coefficients.Upper(keep);

% The mono outcome was negated before fitting, so positive mono coefficients
% denote worse disease; positive co-infection coefficients denote protection.
Effect = estimate;
CI_Lower = lower;
CI_Upper = upper;

Condition = condition;
Strain = normalize_strain(condition, mode);
P_Value = coefficients.pValue(keep);
Q_Value = bh_fdr_local(P_Value);
N_Animals = zeros(numel(condition), 1);
N_Experiments = zeros(numel(condition), 1);
N_Score_Rows = zeros(numel(condition), 1);
N_Terminal_Event_Zero_Rows = zeros(numel(condition), 1);
N_Death_or_Humane_Endpoint_Events = zeros(numel(condition), 1);
for i = 1:numel(condition)
    one = raw(raw.cdiffstrain == condition(i), :);
    N_Animals(i) = numel(unique(one.animal_id));
    N_Experiments(i) = numel(unique(one.experiment));
    N_Score_Rows(i) = height(one);
    N_Terminal_Event_Zero_Rows(i) = sum(ismissing(one.relweight) & ...
        one.death == 1);
    ids = unique(one.animal_id);
    for j = 1:numel(ids)
        N_Death_or_Humane_Endpoint_Events(i) = ...
            N_Death_or_Humane_Endpoint_Events(i) + ...
            any(one.death(one.animal_id == ids(j)) == 1);
    end
end
Reference = repmat(reference, numel(condition), 1);
Effect_Definition = repmat("", numel(condition), 1);
if mode == "mono"
    Effect_Definition(:) = ...
        "Positive = worse terminal-event-zero disease-course score than UI";
else
    Effect_Definition(:) = ...
        "Positive = better terminal-event-zero disease-course score than VPI";
end
Model = repmat("Day 0 retained; terminal-event missing weight assigned score zero; fixed condition; random day and experiment-condition-mouse intercepts", ...
    numel(condition), 1);
Zero_Interpretation = repmat("Composite endpoint score value, not a measured body weight", ...
    numel(condition), 1);
effect_tbl = table(Strain, Condition, Reference, Effect, CI_Lower, CI_Upper, ...
    P_Value, Q_Value, N_Animals, N_Experiments, ...
    N_Score_Rows, N_Terminal_Event_Zero_Rows, ...
    N_Death_or_Humane_Endpoint_Events, Effect_Definition, ...
    Zero_Interpretation, Model);
end

function survival_tbl = compute_survival_comparisons(raw, reference, mode)
comparators = sort(unique(raw.cdiffstrain(raw.cdiffstrain ~= reference)));
n = numel(comparators);
Strain = strings(n, 1);
Condition = comparators;
Reference = repmat(reference, n, 1);
Shared_Experiments = zeros(n, 1);
Reference_N = zeros(n, 1);
Reference_Events = zeros(n, 1);
Comparator_N = zeros(n, 1);
Comparator_Events = zeros(n, 1);
Reference_Event_Fraction = nan(n, 1);
Comparator_Event_Fraction = nan(n, 1);
Event_Risk_Difference = nan(n, 1);
Chi_Square = nan(n, 1);
P_Value = nan(n, 1);

for i = 1:n
    shared = intersect(unique(raw.experiment(raw.cdiffstrain == reference)), ...
        unique(raw.experiment(raw.cdiffstrain == comparators(i))));
    observed_minus_expected = 0;
    variance_total = 0;
    ref_time = zeros(0, 1);
    ref_event = false(0, 1);
    cmp_time = zeros(0, 1);
    cmp_event = false(0, 1);
    for j = 1:numel(shared)
        one_experiment = raw.experiment == shared(j);
        [one_ref_time, one_ref_event] = survival_records( ...
            raw(one_experiment, :), reference);
        [one_cmp_time, one_cmp_event] = survival_records( ...
            raw(one_experiment, :), comparators(i));
        [oe, variance] = logrank_components(one_ref_time, one_ref_event, ...
            one_cmp_time, one_cmp_event);
        observed_minus_expected = observed_minus_expected + oe;
        variance_total = variance_total + variance;
        ref_time = [ref_time; one_ref_time]; %#ok<AGROW>
        ref_event = [ref_event; one_ref_event]; %#ok<AGROW>
        cmp_time = [cmp_time; one_cmp_time]; %#ok<AGROW>
        cmp_event = [cmp_event; one_cmp_event]; %#ok<AGROW>
    end
    Strain(i) = normalize_strain(comparators(i), mode);
    Shared_Experiments(i) = numel(shared);
    Reference_N(i) = numel(ref_event);
    Reference_Events(i) = sum(ref_event);
    Comparator_N(i) = numel(cmp_event);
    Comparator_Events(i) = sum(cmp_event);
    if Reference_N(i) > 0
        Reference_Event_Fraction(i) = Reference_Events(i) / Reference_N(i);
    end
    if Comparator_N(i) > 0
        Comparator_Event_Fraction(i) = Comparator_Events(i) / Comparator_N(i);
    end
    if mode == "mono"
        Event_Risk_Difference(i) = Comparator_Event_Fraction(i) - ...
            Reference_Event_Fraction(i);
    else
        Event_Risk_Difference(i) = Reference_Event_Fraction(i) - ...
            Comparator_Event_Fraction(i);
    end
    if variance_total > 0
        Chi_Square(i) = observed_minus_expected^2 / variance_total;
        P_Value(i) = 1 - chi2cdf(Chi_Square(i), 1);
    end
end
Q_Value = bh_fdr_local(P_Value);
Effect_Definition = repmat("", n, 1);
if mode == "mono"
    Effect_Definition(:) = ...
        "Risk difference = strain event fraction minus matched UI";
else
    Effect_Definition(:) = ...
        "Risk difference = matched VPI event fraction minus strain+VPI";
end
Inference = repmat("Log-rank test stratified by contemporaneous experiment; NA when no event-time variance", ...
    n, 1);
Event_Definition = repmat("Source death field = observed death or humane-endpoint removal", ...
    n, 1);
survival_tbl = table(Strain, Condition, Reference, Shared_Experiments, ...
    Reference_N, Reference_Events, Comparator_N, Comparator_Events, ...
    Reference_Event_Fraction, Comparator_Event_Fraction, ...
    Event_Risk_Difference, Chi_Square, P_Value, Q_Value, ...
    Effect_Definition, Event_Definition, Inference);
end

function [times, events] = survival_records(tbl, condition)
sub = tbl(tbl.cdiffstrain == condition, :);
ids = unique(sub.animal_id);
times = nan(numel(ids), 1);
events = false(numel(ids), 1);
for i = 1:numel(ids)
    one = sortrows(sub(sub.animal_id == ids(i), :), 'day');
    event_index = find(one.death == 1, 1, 'first');
    if isempty(event_index)
        times(i) = max(one.day);
    else
        times(i) = one.day(event_index);
        events(i) = true;
    end
end
end

function [observed_minus_expected, variance_total] = ...
    logrank_components(time1, event1, time2, event2)
event_times = unique([time1(event1); time2(event2)]);
observed_minus_expected = 0;
variance_total = 0;
for t = reshape(event_times, 1, [])
    risk1 = sum(time1 >= t);
    risk2 = sum(time2 >= t);
    events1 = sum(time1 == t & event1);
    events2 = sum(time2 == t & event2);
    risk_total = risk1 + risk2;
    events_total = events1 + events2;
    if risk_total <= 1 || events_total == 0
        continue;
    end
    expected1 = events_total * risk1 / risk_total;
    variance = risk1 * risk2 * events_total * ...
        (risk_total - events_total) / ...
        (risk_total^2 * (risk_total - 1));
    observed_minus_expected = observed_minus_expected + ...
        (events1 - expected1);
    variance_total = variance_total + variance;
end
end

function q = bh_fdr_local(p)
q = nan(size(p));
valid = find(~isnan(p));
if isempty(valid); return; end
[sorted_p, order] = sort(p(valid));
m = numel(sorted_p);
sorted_q = sorted_p .* m ./ (1:m)';
for i = m-1:-1:1
    sorted_q(i) = min(sorted_q(i), sorted_q(i+1));
end
sorted_q = min(sorted_q, 1);
q(valid(order)) = sorted_q;
end

function strain = normalize_strain(condition, mode)
strain = string(condition);
if mode == "protection"
    strain = erase(strain, ".vpi");
end
strain = upper(strrep(strain, '.', '-'));
end

function provenance = build_analysis_provenance(rho, p, n)
Output = ["Mono disease-course scores"; "Co-infection protection scores"; ...
    "Mono survival comparisons"; "Co-infection survival comparisons"; ...
    "Disease/protection score association"; "Focused weight timelines"];
Source = ["Virulence_screen_clean_table.csv and ProtectionScreen_CDI_mouse.csv"; ...
    "ProtectionScreen_CDI_mouse.csv"; ...
    "Virulence_screen_clean_table.csv and ProtectionScreen_CDI_mouse.csv"; ...
    "ProtectionScreen_CDI_mouse.csv"; ...
    "Joined mono disease-course and co-infection protection score tables"; ...
    "ProtectionScreen_CDI_mouse.csv"];
Method = ["Permanent terminal-event-zero score; day 0 retained; fixed condition; random day and experiment-condition-mouse intercepts"; ...
    "Permanent terminal-event-zero score; day 0 retained; fixed condition; random day and experiment-condition-mouse intercepts"; ...
    "Event/censor time from recorded death or humane-endpoint removal and last observation; experiment-stratified log-rank"; ...
    "Event/censor time from recorded death or humane-endpoint removal and last observation; experiment-stratified log-rank"; ...
    "Spearman correlation across strains"; ...
    "Individual observed weights and survivor-observed daily means; a cross marks the final measured weight before a terminal event"];
Missing_Data_Rule = ["Missing weight on recorded terminal-event row assigned score zero; zero is not measured weight"; ...
    "Missing weight on recorded terminal-event row assigned score zero; zero is not measured weight"; ...
    "Terminal events analyzed as recorded"; ...
    "Terminal events analyzed as recorded"; ...
    "Complete joined strain effects"; ...
    "No terminal weight imputation; crosses mark the last measured weight before death or humane-endpoint removal"];
Result = ["See mono_disease_course_scores.csv"; ...
    "See coinfection_protection_scores.csv"; ...
    "See mono_survival_comparisons.csv"; ...
    "See coinfection_survival_comparisons.csv"; ...
    sprintf('rho=%.6g; p=%.6g; n=%d', rho, p, n); ...
    "Figure 2C-D diagnostic timeline panels"];
provenance = table(Output, Source, Method, Missing_Data_Rule, Result);
end

function write_summary(file, weight, mono_survival, protection_survival, ...
    rho, p, mono_model, protection_model)
fid = fopen(file, 'w');
assert(fid >= 0, 'Unable to open summary file: %s', file);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Mouse disease-course score and survival reanalysis\n\n');
fprintf(fid, 'Diarrhea excluded because mono-colonization diarrhea was not recorded for ST1-68.\n');
fprintf(fid, 'Permanent primary setup: missing weight on a recorded terminal-event row is assigned zero in the composite disease-course score.\n');
fprintf(fid, 'The score zero is not a measured body weight and is never plotted as one.\n');
fprintf(fid, 'The source field named death includes both observed deaths and humane-endpoint removals.\n');
fprintf(fid, 'Trajectory crosses mark the final measured weight before a terminal event; they are not terminal-weight measurements.\n');
fprintf(fid, 'Trajectory/survival animal identity = experiment | treatment | condition | mouse number.\n');
fprintf(fid, 'Score-model random identity = experiment_condition_mouse.\n\n');
fprintf(fid, 'Score model formula: relweight score ~ condition + (1|day) + (1|experiment-condition-mouse).\n');
fprintf(fid, 'Day 0 is retained. Mono model observations: %d; protection model observations: %d.\n', ...
    mono_model.NumObservations, protection_model.NumObservations);
fprintf(fid, 'Mono disease effect vs co-infection protection effect: Spearman rho = %.6f, p = %.6g, n = %d strains.\n\n', ...
    rho, p, height(weight));

key = ["ST1-75", "ST1-68", "ST1-49"];
fprintf(fid, 'Prespecified narrative strains (terminal-event-zero score differences):\n');
for i = 1:numel(key)
    row = weight(weight.Strain == key(i), :);
    if isempty(row); continue; end
    fprintf(fid, '%s: mono disease effect %.6f (95%% CI %.6f to %.6f); protection effect %.6f (95%% CI %.6f to %.6f).\n', ...
        key(i), row.Mono_Disease_Effect, row.Mono_Disease_CI_Lower, ...
        row.Mono_Disease_CI_Upper, row.Protection_Effect, ...
        row.Protection_CI_Lower, row.Protection_CI_Upper);
end

fprintf(fid, '\nPrespecified narrative strains (survival):\n');
for i = 1:numel(key)
    mono = mono_survival(mono_survival.Strain == key(i), :);
    protection = protection_survival(protection_survival.Strain == key(i), :);
    if ~isempty(mono)
        fprintf(fid, '%s mono: %d/%d events vs %d/%d contemporaneous UI events; stratified log-rank p = %.6g.\n', ...
            key(i), mono.Comparator_Events, mono.Comparator_N, ...
            mono.Reference_Events, mono.Reference_N, mono.P_Value);
    end
    if ~isempty(protection)
        fprintf(fid, '%s + VPI: %d/%d events vs %d/%d contemporaneous VPI events; stratified log-rank p = %.6g.\n', ...
            key(i), protection.Comparator_Events, protection.Comparator_N, ...
            protection.Reference_Events, protection.Reference_N, ...
            protection.P_Value);
    end
end
end

function make_diagnostic_figure(figure_dir, weight, coinfection, rho, p)
fig = figure('Color', 'w', 'Position', [50 50 1500 1050]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile;
hold(ax, 'on');
xline(ax, 0, ':', 'Color', [0.5 0.5 0.5], ...
    'HandleVisibility', 'off');
yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], ...
    'HandleVisibility', 'off');
for i = 1:height(weight)
    plot(ax, [weight.Mono_Disease_CI_Lower(i), ...
        weight.Mono_Disease_CI_Upper(i)], ...
        [weight.Protection_Effect(i), weight.Protection_Effect(i)], '-', ...
        'Color', [0.82 0.82 0.82], 'HandleVisibility', 'off');
    plot(ax, [weight.Mono_Disease_Effect(i), ...
        weight.Mono_Disease_Effect(i)], ...
        [weight.Protection_CI_Lower(i), weight.Protection_CI_Upper(i)], '-', ...
        'Color', [0.82 0.82 0.82], 'HandleVisibility', 'off');
end
scatter(ax, weight.Mono_Disease_Effect, weight.Protection_Effect, 42, ...
    [0.32 0.32 0.32], 'filled');
label_key_strains(ax, weight);
xlabel(ax, 'Mono-colonization disease-course score (+ = worse)');
ylabel(ax, 'Co-infection protection score (+ = better)');
title(ax, sprintf('Disease-course scores: Spearman \\rho = %.2f, p = %.3g', rho, p));
box(ax, 'on');
panel_label(ax, 'A');

ax = nexttile;
plot_key_weight_effects(ax, weight);
panel_label(ax, 'B');

ax = nexttile;
timeline_scale = focused_timeline_scale(coinfection, ...
    ["st1.75.vpi", "st1.68.vpi"]);
plot_weight_timelines(ax, coinfection, "st1.75.vpi", ...
    [0.10 0.45 0.75], timeline_scale);
panel_label(ax, 'C');

ax = nexttile;
plot_weight_timelines(ax, coinfection, "st1.68.vpi", ...
    [0.88 0.38 0.12], timeline_scale);
panel_label(ax, 'D');

sgtitle(layout, ...
    'Mouse discovery screen: disease-course scores and observed trajectories', ...
    'FontWeight', 'bold', 'FontSize', 15);
exportgraphics(fig, fullfile(figure_dir, ...
    'figure2_mouse_weight_survival.png'), 'Resolution', 300);
close(fig);
end

function label_key_strains(ax, weight)
key = ["ST1-75", "ST1-68", "ST1-49"];
colors = [0.10 0.45 0.75; 0.88 0.38 0.12; 0.45 0.18 0.55];
offsets = [0.25 0.6; 0.25 -0.8; 0.25 0.6];
for i = 1:numel(key)
    index = find(weight.Strain == key(i), 1);
    if isempty(index); continue; end
    scatter(ax, weight.Mono_Disease_Effect(index), ...
        weight.Protection_Effect(index), 95, colors(i, :), 'filled', ...
        'MarkerEdgeColor', 'k');
    text(ax, weight.Mono_Disease_Effect(index) + offsets(i, 1), ...
        weight.Protection_Effect(index) + offsets(i, 2), key(i), ...
        'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', 'w');
end
end

function plot_key_weight_effects(ax, weight)
key = ["ST1-49", "ST1-68", "ST1-75"];
sub = weight(ismember(weight.Strain, key), :);
[~, order] = ismember(key, sub.Strain);
sub = sub(order, :);
y = (1:height(sub))';
hold(ax, 'on');
xline(ax, 0, ':', 'Color', [0.5 0.5 0.5], ...
    'HandleVisibility', 'off');
for i = 1:height(sub)
    plot(ax, [sub.Mono_Disease_CI_Lower(i), sub.Mono_Disease_CI_Upper(i)], ...
        [y(i)+0.12, y(i)+0.12], '-', 'Color', [0.76 0.22 0.18], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(ax, [sub.Protection_CI_Lower(i), sub.Protection_CI_Upper(i)], ...
        [y(i)-0.12, y(i)-0.12], '-', 'Color', [0.10 0.45 0.75], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
end
scatter(ax, sub.Mono_Disease_Effect, y+0.12, 55, [0.76 0.22 0.18], ...
    'filled', 'DisplayName', 'Mono disease effect');
scatter(ax, sub.Protection_Effect, y-0.12, 55, [0.10 0.45 0.75], ...
    'filled', 'DisplayName', 'Protection effect');
yticks(ax, y);
yticklabels(ax, sub.Strain);
ylim(ax, [0.5 height(sub)+0.5]);
xlabel(ax, 'Terminal-event-zero disease-course score difference');
title(ax, 'Prespecified narrative strains');
legend(ax, 'Location', 'best');
box(ax, 'on');
end

function scale = focused_timeline_scale(raw, comparators)
selected = raw(raw.cdiffstrain == "vpi" | ...
    ismember(raw.cdiffstrain, comparators), :);
values = selected.relweight(~ismissing(selected.relweight));
assert(~isempty(values), 'Focused timeline conditions have no weight data.');
scale.Data_Lower = floor((min(values) - 2) / 5) * 5;
scale.Data_Upper = ceil(max(values) / 5) * 5;
scale.Axis_Upper = scale.Data_Upper + 1;
end

function plot_weight_timelines(ax, raw, comparator, comparator_color, scale)
reference = "vpi";
shared_experiments = intersect( ...
    unique(raw.experiment(raw.cdiffstrain == reference)), ...
    unique(raw.experiment(raw.cdiffstrain == comparator)));
sub = raw(ismember(raw.experiment, shared_experiments) & ...
    (raw.cdiffstrain == reference | raw.cdiffstrain == comparator), :);
conditions = [reference, comparator];
colors = [0.25 0.25 0.25; comparator_color];
labels = ["VPI", normalize_strain(comparator, "protection") + "+VPI"];
line_handles = gobjects(2, 1);
event_color = [0.78 0.10 0.10];

hold(ax, 'on');
days_all = sort(unique(sub.day));
x_limits = [min(days_all), max(days_all)];

for g = 1:2
    group = sub(sub.cdiffstrain == conditions(g), :);
    ids = unique(group.animal_id);
    individual_color = 0.72 + 0.28 * colors(g, :);
    for i = 1:numel(ids)
        one = sortrows(group(group.animal_id == ids(i), :), 'day');
        observed = ~ismissing(one.relweight);
        plot(ax, one.day(observed), one.relweight(observed), '-', ...
            'Color', individual_color, 'LineWidth', 0.7, ...
            'HandleVisibility', 'off');
        endpoint = find(one.death == 1, 1, 'first');
        if ~isempty(endpoint)
            last_observed = find(observed, 1, 'last');
            assert(~isempty(last_observed), ...
                'Terminal-event animal has no observed weight.');
            plot(ax, one.day(last_observed), one.relweight(last_observed), ...
                'x', 'Color', event_color, 'MarkerSize', 8, ...
                'LineWidth', 1.6, 'HandleVisibility', 'off');
        end
    end

    observed_mean = nan(numel(days_all), 1);
    for d = 1:numel(days_all)
        values = group.relweight(group.day == days_all(d));
        observed_mean(d) = mean(values, 'omitnan');
    end
    line_handles(g) = plot(ax, days_all, observed_mean, '-', ...
        'Color', colors(g, :), 'LineWidth', 2.4, ...
        'DisplayName', sprintf('%s observed mean (n=%d)', ...
        labels(g), numel(ids)));
end

event_handle = plot(ax, nan, nan, 'x', 'Color', event_color, ...
    'MarkerSize', 8, 'LineWidth', 1.6, ...
    'DisplayName', 'Last observed weight before death/humane endpoint');

xlim(ax, x_limits);
ylim(ax, [scale.Data_Lower scale.Axis_Upper]);
yticks(ax, scale.Data_Lower:5:scale.Data_Upper);
xlabel(ax, 'Day');
ylabel(ax, 'Relative weight (%)');
if isscalar(shared_experiments)
    experiment_text = shared_experiments(1);
else
    experiment_text = sprintf('%d contemporaneous experiments', ...
        numel(shared_experiments));
end
title(ax, sprintf('%s versus VPI (%s)', labels(2), experiment_text));
legend(ax, [line_handles; event_handle], 'Location', 'southwest');
box(ax, 'on');
end

function panel_label(ax, label)
text(ax, 0.01, 0.99, label, 'Units', 'normalized', 'FontSize', 15, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'top');
end
