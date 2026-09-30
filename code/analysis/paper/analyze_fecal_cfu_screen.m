%% Analyze fecal C. difficile burden across the mono- and co-colonization screens
% Mono-colonization values are strain-specific CFU/g. In the co-colonization
% screen, plating reports total C. difficile CFU/g and cannot resolve the ST1
% candidate from VPI10463. Early samples provide the least survivor-conditioned
% comparison. Late samples are retained descriptively among mice with an
% observed fecal CFU measurement. Zero values are treated as non-detects and
% are never replaced by a numerical CFU value for analysis.

script_dir = fileparts(mfilename('fullpath'));
base_dir = fileparts(fileparts(script_dir));
project_dir = fileparts(base_dir);
addpath(genpath(fullfile(base_dir, 'analysis')));

figure_dir = fullfile(project_dir, 'results', 'figures');
table_dir = fullfile(project_dir, 'results', 'tables', 'main');
if ~exist(figure_dir, 'dir'); mkdir(figure_dir); end
if ~exist(table_dir, 'dir'); mkdir(table_dir); end
set(0, 'DefaultFigureVisible', 'off');

virulence_file = project_data_file('processed', 'mouse', 'scores', ...
    'Virulence_screen_clean_table.csv');
protection_file = project_data_file('processed', 'mouse', 'scores', ...
    'ProtectionScreen_CDI_mouse.csv');
score_file = fullfile(project_dir,'results','tables','mouse_weight_survival_reanalysis','primary_terminal_event_zero_scores.csv');

virulence = readtable(virulence_file, 'TextType', 'string');
protection = readtable(protection_file, 'TextType', 'string');
scores = readtable(score_file, 'TextType', 'string');
validate_source(virulence, "virulence");
validate_source(protection, "protection");

virulence.Source = repmat("Mono-colonization screen", height(virulence), 1);
protection.Source = repmat("Co-infection screen", height(protection), 1);
virulence.ExperimentID = "virulence|" + lower(virulence.experiment);
protection.ExperimentID = "protection|" + lower(protection.experiment);

schedule = build_schedule();
virulence_selected = select_scheduled_cfu(virulence, schedule, "virulence");
protection_selected = select_scheduled_cfu(protection, schedule, "protection");

candidate_order = sort_st1_labels(unique(scores.Strain));
candidate_order = candidate_order(startsWith(candidate_order, "ST1-"));
display_order = [candidate_order; "VPI10463"];

mono_v = select_mono_conditions(virulence_selected, candidate_order);
mono_p = select_mono_conditions(protection_selected, candidate_order);
mono = [mono_v; mono_p];
mono.Context = repmat("Mono-colonization", height(mono), 1);
mono.MeasurementMeaning = repmat("Strain-specific fecal C. difficile CFU/g", ...
    height(mono), 1);

co = select_coinfection_conditions(protection_selected, candidate_order);
co.Context = repmat("Co-colonization", height(co), 1);
co.MeasurementMeaning = repmat("Total fecal C. difficile CFU/g; strains unresolved", ...
    height(co), 1);

observations = [standardize_output(mono); standardize_output(co)];
observations.Detected = observations.cdiffcfu > 0 & ~isnan(observations.cdiffcfu);
observations.Log10CFU = nan(height(observations), 1);
positive = observations.Detected;
observations.Log10CFU(positive) = log10(observations.cdiffcfu(positive));

assert(all(observations.cdiffcfu(~isnan(observations.cdiffcfu)) >= 0), ...
    'Negative CFU values are not biologically interpretable.');
assert(all(ismember(candidate_order, unique(mono.Candidate))), ...
    'At least one matched ST1 strain is absent from mono-colonization CFU data.');
assert(all(ismember(candidate_order, unique(co.Candidate))), ...
    'At least one matched ST1 strain is absent from co-colonization CFU data.');

summary_table = summarize_cfu(observations, display_order);
global_tests = run_global_tests(observations);
global_tests.QValue = bh_adjust(global_tests.PValue);
association_table = run_score_associations(observations, scores, candidate_order);
source_provenance = table( ...
    ["Mono-colonization CFU source"; "Co-infection CFU source"; ...
     "Candidate disease/protection scores"], ...
    [string(virulence_file); string(protection_file); string(score_file)], ...
    ["Processed mouse table; strain-specific CFU for strain-alone conditions"; ...
     "Processed mouse table; strain-specific CFU for strain-alone controls and total CFU for candidate+VPI conditions"; ...
     "Primary terminal-event-zero score estimates used only for strain-level association tests"], ...
    'VariableNames', {'Role', 'SourceFile', 'Use'});

writetable(observations, fullfile(table_dir, ...
    'figureS1_fecal_cfu_observations.csv'));
writetable(summary_table, fullfile(table_dir, ...
    'figureS1_fecal_cfu_strain_summary.csv'));
writetable(global_tests, fullfile(table_dir, ...
    'figureS1_fecal_cfu_global_tests.csv'));
writetable(association_table, fullfile(table_dir, ...
    'figureS1_fecal_cfu_score_associations.csv'));
writetable(source_provenance, fullfile(table_dir, ...
    'figureS1_fecal_cfu_source_provenance.csv'));

figure_file = fullfile(figure_dir, 'figureS1_fecal_cfu_screen.png');
generate_figure(observations, display_order, figure_file);
write_summary(fullfile(table_dir, 'figureS1_fecal_cfu_analysis_summary.txt'), ...
    observations, summary_table, global_tests, association_table, figure_file);

fprintf('Wrote %s\n', figure_file);
disp(global_tests);
disp(association_table);

function validate_source(data, label)
required = ["experiment", "day", "mouse", "tx", "cdiffstrain", ...
    "cdiffcfu", "death"];
missing = setdiff(required, string(data.Properties.VariableNames));
assert(isempty(missing), '%s table lacks required variables: %s', ...
    label, strjoin(missing, ', '));
end

function schedule = build_schedule()
virulence_experiments = ["jwk90"; "ks1"; "ks10"; "ks11"; "ks14"; ...
    "ks15"; "ks2"; "ks25"; "ks4"; "ks54"; "ks6"; "ks8"; "qd03"];
virulence_early = ones(numel(virulence_experiments), 1);
virulence_late = 7 * ones(numel(virulence_experiments), 1);
virulence_late(virulence_experiments == "ks6") = 6;

protection_experiments = ["ks52"; "ks53"; "ks56"; "ks57"; "ks58"; ...
    "ks59"; "ks60"];
protection_early = [1; 1; 1; 1; 2; 2; 2];
protection_late = [7; 7; 6; 7; 6; 7; 7];

schedule = table( ...
    [repmat("virulence", numel(virulence_experiments), 1); ...
     repmat("protection", numel(protection_experiments), 1)], ...
    [virulence_experiments; protection_experiments], ...
    [virulence_early; protection_early], ...
    [virulence_late; protection_late], ...
    'VariableNames', {'Dataset', 'Experiment', 'EarlyDay', 'LateDay'});
end

function selected = select_scheduled_cfu(data, schedule, dataset)
data.experiment = lower(string(data.experiment));
data.cdiffstrain = lower(string(data.cdiffstrain));
data.tx = lower(string(data.tx));
subschedule = schedule(schedule.Dataset == dataset, :);
available = unique(data.experiment);
expected = subschedule.Experiment;
unexpected = setdiff(available, [expected; "ks65"]);
assert(isempty(unexpected), 'Unexpected %s experiment(s): %s', ...
    dataset, strjoin(unexpected, ', '));

pieces = cell(height(subschedule) * 2, 1);
k = 0;
for i = 1:height(subschedule)
    experiment = subschedule.Experiment(i);
    for phase = ["Early", "Late"]
        k = k + 1;
        if phase == "Early"
            target_day = subschedule.EarlyDay(i);
        else
            target_day = subschedule.LateDay(i);
        end
        block = data(data.experiment == experiment & data.day == target_day, :);
        assert(~isempty(block), '%s %s has no scheduled day-%g rows.', ...
            dataset, experiment, target_day);
        block.Phase = repmat(phase, height(block), 1);
        block.ScheduledDay = repmat(target_day, height(block), 1);
        pieces{k} = block;
    end
end
selected = vertcat(pieces{:});
end

function selected = select_mono_conditions(data, candidate_order)
raw = lower(string(data.cdiffstrain));
is_st1_alone = ~cellfun(@isempty, regexp(cellstr(raw), '^st1\.[0-9]+$', 'once'));
is_vpi_alone = raw == "vpi" | raw == "vpi10463";
selected = data(is_st1_alone | is_vpi_alone, :);
selected.Candidate = normalize_candidate(selected.cdiffstrain);
selected = selected(ismember(selected.Candidate, [candidate_order; "VPI10463"]), :);
end

function selected = select_coinfection_conditions(data, candidate_order)
raw = lower(string(data.cdiffstrain));
is_candidate_vpi = ~cellfun(@isempty, ...
    regexp(cellstr(raw), '^st1\.[0-9]+\.vpi$', 'once'));
is_vpi_alone = raw == "vpi" | raw == "vpi10463";
selected = data(is_candidate_vpi | is_vpi_alone, :);
selected.Candidate = normalize_candidate(selected.cdiffstrain);
selected = selected(ismember(selected.Candidate, [candidate_order; "VPI10463"]), :);
end

function labels = normalize_candidate(raw)
raw = lower(string(raw));
labels = strings(size(raw));
for i = 1:numel(raw)
    if raw(i) == "vpi" || raw(i) == "vpi10463"
        labels(i) = "VPI10463";
    else
        token = regexp(char(raw(i)), '^st1\.([0-9]+)', 'tokens', 'once');
        assert(~isempty(token), 'Cannot normalize condition %s.', raw(i));
        labels(i) = "ST1-" + string(token{1});
    end
end
end

function order = sort_st1_labels(labels)
labels = unique(string(labels), 'stable');
numbers = nan(size(labels));
for i = 1:numel(labels)
    token = regexp(char(labels(i)), '^ST1-([0-9]+)$', 'tokens', 'once');
    if ~isempty(token); numbers(i) = str2double(token{1}); end
end
[~, idx] = sort(numbers);
order = labels(idx);
order = order(:);
end

function output = standardize_output(data)
output = data(:, {'Context', 'Phase', 'Source', 'ExperimentID', ...
    'experiment', 'ScheduledDay', 'day', 'mouse', 'tx', 'cdiffstrain', ...
    'Candidate', 'cdiffcfu', 'death', 'MeasurementMeaning'});
output = renamevars(output, {'experiment', 'day', 'mouse', 'tx', ...
    'cdiffstrain', 'death'}, {'Experiment', 'ObservedDay', 'Mouse', ...
    'Protocol', 'SourceCondition', 'TerminalEventRow'});
end

function summary_table = summarize_cfu(observations, display_order)
contexts = ["Mono-colonization", "Co-colonization"];
phases = ["Early", "Late"];
rows = cell(0, 1);
for c = contexts
    for p = phases
        for s = display_order(:)'
            block = observations(observations.Context == c & ...
                observations.Phase == p & observations.Candidate == s, :);
            if isempty(block); continue; end
            measured = ~isnan(block.cdiffcfu);
            detected = measured & block.cdiffcfu > 0;
            positive_values = block.cdiffcfu(detected);
            if isempty(positive_values)
                median_positive = NaN;
                geometric_mean_positive = NaN;
                min_positive = NaN;
                max_positive = NaN;
            else
                median_positive = median(positive_values);
                geometric_mean_positive = 10 ^ mean(log10(positive_values));
                min_positive = min(positive_values);
                max_positive = max(positive_values);
            end
            row = table(c, p, s, height(block), sum(measured), sum(detected), ...
                sum(measured & block.cdiffcfu == 0), ...
                sum(detected) / max(sum(measured), 1), median_positive, ...
                geometric_mean_positive, min_positive, max_positive, ...
                numel(unique(block.ExperimentID)), ...
                'VariableNames', {'Context', 'Phase', 'Candidate', ...
                'RowsScheduled', 'NMeasured', 'NDetected', 'NNonDetect', ...
                'DetectionFraction', 'MedianPositiveCFU', ...
                'GeometricMeanPositiveCFU', 'MinimumPositiveCFU', ...
                'MaximumPositiveCFU', 'NExperiments'});
            rows{end+1, 1} = row; %#ok<AGROW>
        end
    end
end
summary_table = vertcat(rows{:});
end

function tests = run_global_tests(observations)
contexts = ["Mono-colonization", "Co-colonization"];
phases = ["Early", "Late"];
rows = cell(0, 1);
for c = contexts
    for p = phases
        block = observations(observations.Context == c & ...
            observations.Phase == p & observations.Detected, :);
        block.CandidateCat = categorical(block.Candidate);
        block.SourceCat = categorical(block.Source);
        block.ExperimentCat = categorical(block.ExperimentID);
        if c == "Mono-colonization"
            reduced = fitlme(block, ...
                'Log10CFU ~ SourceCat + (1|ExperimentCat)', 'FitMethod', 'ML');
            full = fitlme(block, ...
                'Log10CFU ~ SourceCat + CandidateCat + (1|ExperimentCat)', ...
                'FitMethod', 'ML');
            adjustment = "source fixed effect; experiment random intercept";
        else
            reduced = fitlme(block, ...
                'Log10CFU ~ 1 + (1|ExperimentCat)', 'FitMethod', 'ML');
            full = fitlme(block, ...
                'Log10CFU ~ CandidateCat + (1|ExperimentCat)', ...
                'FitMethod', 'ML');
            adjustment = "experiment random intercept";
        end
        comparison = compare(reduced, full, 'CheckNesting', true);
        row = table(c, p, height(block), ...
            numel(unique(block.Candidate)), comparison.LRStat(2), ...
            comparison.deltaDF(2), comparison.pValue(2), adjustment, ...
            string(full.Formula.char), ...
            p + " positive-CFU comparison; late phase is survivor-conditioned", ...
            'VariableNames', {'Context', 'Phase', 'NPositiveMeasurements', ...
            'NConditions', 'LikelihoodRatio', 'DegreesFreedom', 'PValue', ...
            'Adjustment', 'FullModel', 'InterpretationConstraint'});
        rows{end+1, 1} = row; %#ok<AGROW>
    end
end
tests = vertcat(rows{:});
end

function associations = run_score_associations(observations, scores, candidate_order)
contexts = ["Mono-colonization", "Co-colonization"];
phases = ["Early", "Late"];
rows = cell(0, 1);
for c = contexts
    for p = phases
        block = observations(observations.Context == c & ...
            observations.Phase == p & observations.Detected & ...
            ismember(observations.Candidate, candidate_order), :);
        [g, candidate, experiment] = findgroups(block.Candidate, block.ExperimentID);
        experiment_median = splitapply(@median, block.Log10CFU, g);
        by_experiment = table(candidate, experiment, experiment_median, ...
            'VariableNames', {'Strain', 'Experiment', 'MedianLog10CFU'});
        [g2, strain] = findgroups(by_experiment.Strain);
        strain_median = splitapply(@median, by_experiment.MedianLog10CFU, g2);
        strain_table = table(strain, strain_median, ...
            'VariableNames', {'Strain', 'ExperimentBalancedMedianLog10CFU'});
        joined = innerjoin(strain_table, scores, 'Keys', 'Strain');
        if c == "Mono-colonization"
            outcome = joined.Mono_Disease_Effect;
            outcome_name = "Primary mono-colonization disease effect";
        else
            outcome = joined.Protection_Effect;
            outcome_name = "Primary co-infection protection effect";
        end
        [rho, pvalue] = corr(joined.ExperimentBalancedMedianLog10CFU, ...
            outcome, 'Type', 'Spearman', 'Rows', 'complete');
        row = table(c, p, outcome_name, height(joined), rho, pvalue, ...
            "Median of within-experiment positive-log10-CFU medians; non-detects retained separately", ...
            p + " association; late phase is survivor-conditioned", ...
            'VariableNames', {'Context', 'Phase', 'Outcome', 'NStrains', ...
            'SpearmanRho', 'PValue', 'BurdenSummary', ...
            'InterpretationConstraint'});
        rows{end+1, 1} = row; %#ok<AGROW>
    end
end
associations = vertcat(rows{:});
end

function generate_figure(observations, display_order, output_file)
colors = struct('default', [0.45 0.49 0.54], ...
    'ST1_75', [0.10 0.42 0.68], 'ST1_68', [0.90 0.45 0.13], ...
    'ST1_49', [0.42 0.42 0.42], 'ST1_6', [142 90 169] ./ 255, ...
    'VPI10463', [0.74 0.18 0.20]);
fig = figure('Color', 'w', 'Position', [80 80 1800 1180]);
layout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
layout.OuterPosition = [0.035 0.095 0.95 0.88];

panel_specs = {"Mono-colonization", "Early", ...
        "A  Strain-specific burden after mono-colonization"; ...
    "Mono-colonization", "Late", ...
        "B  Late strain-specific burden among observed mice"; ...
    "Co-colonization", "Early", ...
        "C  Total burden during candidate + VPI10463 co-colonization"; ...
    "Co-colonization", "Late", ...
        "D  Late total burden among observed mice"};

for panel = 1:4
    ax = nexttile;
    hold(ax, 'on');
    context = panel_specs{panel, 1};
    phase = panel_specs{panel, 2};
    block = observations(observations.Context == context & ...
        observations.Phase == phase, :);
    y_positive = block.Log10CFU(block.Detected);
    if isempty(y_positive); positive_floor = 3; else; positive_floor = floor(min(y_positive)); end
    nd_y = positive_floor - 0.75;
    y_min = nd_y - 0.35;
    y_max = max(9, ceil(max(y_positive)) + 0.25);
    patch(ax, [0.4 numel(display_order)+0.6 numel(display_order)+0.6 0.4], ...
        [y_min y_min nd_y+0.22 nd_y+0.22], [0.94 0.94 0.94], ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');

    for i = 1:numel(display_order)
        strain = display_order(i);
        strain_block = block(block.Candidate == strain & ...
            ~isnan(block.cdiffcfu), :);
        if isempty(strain_block); continue; end
        n = height(strain_block);
        offsets = linspace(-0.22, 0.22, max(n, 2));
        offsets = offsets(1:n);
        source_shift = zeros(n, 1);
        source_shift(strain_block.Source == "Co-infection screen" & ...
            context == "Mono-colonization") = 0.035;
        x = i + offsets(:) + source_shift;
        y = strain_block.Log10CFU;
        y(~strain_block.Detected) = nd_y;
        color = strain_color(strain, colors);
        for j = 1:n
            if context == "Mono-colonization" && ...
                    strain_block.Source(j) == "Co-infection screen"
                marker = 'd';
            else
                marker = 'o';
            end
            if strain_block.Detected(j)
                scatter(ax, x(j), y(j), 35, 'Marker', marker, ...
                    'MarkerFaceColor', color, 'MarkerEdgeColor', 'w', ...
                    'LineWidth', 0.45, 'MarkerFaceAlpha', 0.78, ...
                    'HandleVisibility', 'off');
            else
                scatter(ax, x(j), y(j), 34, color, 'x', ...
                    'LineWidth', 1.25, 'HandleVisibility', 'off');
            end
        end
        positive_y = y(strain_block.Detected);
        if ~isempty(positive_y)
            med = median(positive_y);
            plot(ax, [i-0.25 i+0.25], [med med], '-', 'Color', ...
                darken(color, 0.75), 'LineWidth', 2.4, 'HandleVisibility', 'off');
        end
    end

    xline(ax, numel(display_order)-0.5, ':', 'Color', [0.60 0.60 0.60], ...
        'LineWidth', 1.0, 'HandleVisibility', 'off');
    xlim(ax, [0.4 numel(display_order)+0.6]);
    ylim(ax, [y_min y_max]);
    positive_ticks = max(positive_floor, 3):floor(y_max);
    yticks(ax, [nd_y positive_ticks]);
    ylabels = ["ND" compose('10^{%d}', positive_ticks)];
    yticklabels(ax, ylabels);
    xticks(ax, 1:numel(display_order));
    short_labels = replace(display_order, "ST1-", "");
    short_labels(display_order == "VPI10463") = "VPI";
    xticklabels(ax, short_labels);
    xtickangle(ax, 45);
    ylabel(ax, 'Fecal C. difficile (CFU/g)');
    title(ax, panel_specs{panel, 3}, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left');
    subtitle(ax, sprintf('%s sampling day varies by experiment', lower(phase)));
    grid(ax, 'on');
    ax.GridAlpha = 0.10;
    ax.XGrid = 'off';
    ax.Box = 'off';
    ax.FontName = 'Arial';
    ax.FontSize = 10;
end

annotation(fig, 'textbox', [0.045 0.008 0.91 0.055], 'String', ...
    ['Points show observed fecal CFU/g; bars show medians among positive measurements. ' ...
     'Crosses in the ND band denote samples with no recorded colonies; they are not imputed CFU values. ' ...
     'Circles denote the mono-colonization screen and diamonds denote strain-alone controls from the co-infection screen. ' ...
     'Candidate+VPI10463 panels show total C. difficile CFU and cannot resolve strain composition.'], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontName', 'Arial', 'FontSize', 10);
export_vector_asset(fig, output_file);
exportgraphics(fig, output_file, 'Resolution', 300);
close(fig);
end

function color = strain_color(strain, colors)
switch strain
    case "ST1-75"; color = colors.ST1_75;
    case "ST1-68"; color = colors.ST1_68;
    case "ST1-49"; color = colors.ST1_49;
    case "ST1-6"; color = colors.ST1_6;
    case "VPI10463"; color = colors.VPI10463;
    otherwise; color = colors.default;
end
end

function out = darken(color, factor)
out = max(0, min(1, color * factor));
end

function adjusted = bh_adjust(pvalues)
pvalues = pvalues(:);
adjusted = nan(size(pvalues));
valid = ~isnan(pvalues);
[sorted_p, order] = sort(pvalues(valid));
m = numel(sorted_p);
sorted_q = sorted_p .* m ./ (1:m)';
sorted_q = flipud(cummin(flipud(sorted_q)));
sorted_q = min(sorted_q, 1);
valid_indices = find(valid);
adjusted(valid_indices(order)) = sorted_q;
end

function write_summary(filename, observations, summary_table, tests, associations, figure_file)
fid = fopen(filename, 'w');
assert(fid >= 0, 'Unable to create %s.', filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Fecal CFU screen analysis\n');
fprintf(fid, '=========================\n\n');
fprintf(fid, 'Figure: %s\n', figure_file);
fprintf(fid, 'Selected observed rows: %d\n', height(observations));
fprintf(fid, 'Measured CFU values: %d\n', sum(~isnan(observations.cdiffcfu)));
fprintf(fid, 'Source-recorded non-detects: %d\n\n', ...
    sum(observations.cdiffcfu == 0, 'omitnan'));

for phase = ["Early", "Late"]
    st1 = summary_table(summary_table.Context == "Mono-colonization" & ...
        summary_table.Phase == phase & startsWith(summary_table.Candidate, "ST1-"), :);
    fprintf(fid, '%s mono-colonization: %d/%d ST1 strains had at least one positive recovery; %d/%d had all measured samples positive.\n', ...
        phase, sum(st1.NDetected > 0), height(st1), ...
        sum(st1.NDetected == st1.NMeasured), height(st1));
end

st168 = summary_table(summary_table.Candidate == "ST1-68" & ...
    summary_table.Context == "Mono-colonization", :);
fprintf(fid, 'ST1-68 mono-colonization summary\n');
for i = 1:height(st168)
    fprintf(fid, '%s: %d/%d detected; positive range %.6g to %.6g CFU/g; %d experiments.\n', ...
        st168.Phase(i), st168.NDetected(i), st168.NMeasured(i), ...
        st168.MinimumPositiveCFU(i), st168.MaximumPositiveCFU(i), ...
        st168.NExperiments(i));
end

fprintf(fid, '\nGlobal positive-load tests\n');
for i = 1:height(tests)
    fprintf(fid, '%s %s: LR=%.6g, df=%g, p=%.6g, BH q=%.6g, n=%d positive measurements.\n', ...
        tests.Context(i), tests.Phase(i), tests.LikelihoodRatio(i), ...
        tests.DegreesFreedom(i), tests.PValue(i), tests.QValue(i), ...
        tests.NPositiveMeasurements(i));
end

fprintf(fid, '\nAssociations with manuscript scores\n');
for i = 1:height(associations)
    fprintf(fid, '%s %s versus %s: rho=%.6g, p=%.6g, n=%d strains.\n', ...
        associations.Context(i), associations.Phase(i), ...
        associations.Outcome(i), associations.SpearmanRho(i), ...
        associations.PValue(i), associations.NStrains(i));
end

fprintf(fid, ['\nInterpretation constraints: mono-colonization CFU is strain-specific. ' ...
    'Candidate+VPI10463 CFU is total C. difficile burden and cannot identify ' ...
    'which strain contributed the colonies. Late values describe mice with an ' ...
    'observed sample and are conditioned by terminal events and other missingness. ' ...
    'Zero CFU values remain non-detects and were excluded from log10-load models.\n']);
clear cleanup
end
