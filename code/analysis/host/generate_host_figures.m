root=package_root(); tableDir=fullfile(root,'results','tables','main');
figureDir=fullfile(root,'results','figures');
rag=readtable(package_input('KS11_rechallenge_relweight.xlsx'),'TextType','string');
rng(1,'twister'); microbiome=compute_residual_microbiome_stats([]);
validate_rag1(rag);validate_microbiome(microbiome);
[ragEffects,~]=fit_rag1_protection_model(rag);
writetable(ragEffects,fullfile(tableDir,'figure5_rag1_protection_model_effects.csv'));
writetable(microbiome,fullfile(tableDir,'figure5_residual_microbiome_summary.csv'));
make_figure5(fullfile(figureDir,'figure5_host_microbiota_accessory_context.png'),rag,ragEffects,microbiome,[],[],[],[]);
make_figureS3(fullfile(figureDir,'figureS3_contextual_analyses.png'),microbiome);

function validate_rag1(rag)
required = ["group", "mouse_id", "day", "relweight"];
assert(all(ismember(required, string(rag.Properties.VariableNames))), ...
    'RAG1 table is missing a required column.');
groups = ["uninfected_b6", "st175_b6", "st175_rag1ko"];
observedGroups = sort(unique(string(rag.group)));
expectedGroups = sort(groups(:));
assert(isequal(observedGroups(:), expectedGroups), ...
    'Unexpected RAG1 challenge groups.');
for group = groups
    sub = rag(string(rag.group) == group, :);
    assert(numel(unique(sub.mouse_id)) == 10, ...
        'Each RAG1 challenge group must contain ten mice.');
end
assert(all(rag.day >= 0 & rag.day <= 7) && ...
    all(rag.relweight > 0 | isnan(rag.relweight)), ...
    'RAG1 trajectories contain an unexpected day or weight value.');

uncolonized = rag(string(rag.group) == "uninfected_b6", :);
finalScheduledDay = max(uncolonized.day);
uncolonizedMice = unique(uncolonized.mouse_id);
nTerminal = 0;
for i = 1:numel(uncolonizedMice)
    one = uncolonized(uncolonized.mouse_id == uncolonizedMice(i), :);
    observedDays = one.day(~isnan(one.relweight));
    if ~isempty(observedDays) && max(observedDays) < finalScheduledDay
        nTerminal = nTerminal + 1;
    end
end
assert(nTerminal == 9, ...
    'Expected nine terminated previously uncolonized B6 trajectories.');
end


function validate_microbiome(microbiome)
alpha = microbiome(microbiome.Analysis == "Residual_Shannon", :);
bray = microbiome(microbiome.Analysis == "Residual_BrayCurtis", :);
genera = microbiome(microbiome.Analysis == "Residual_Top_Genus", :);
assert(height(alpha) == 4 && isequal(alpha.Day', [0 1 2 7]), ...
    'Expected residual Shannon summaries for days 0, 1, 2, and 7.');
assert(height(bray) == 4 && ~isempty(genera), ...
    'Residual microbiome summary is incomplete.');
assert(min([alpha.Q_Value; bray.Q_Value; genera.Q_Value], [], ...
    'omitnan') >= 0.05, ...
    'A residual-microbiome comparison now survives correction; re-review text.');
end


function [effects, model] = fit_rag1_protection_model(rag)
% Use the same mixed-effects disease-course model structure as Figures 1 and
% 2. The available RAG1 source lacks explicit terminal-event rows, so the
% model retains observed relative weights and leaves all subsequent values
% missing rather than assigning terminal-event score zeros.
groups = ["uninfected_b6", "st175_b6", "st175_rag1ko"];
assert(isequal(sort(unique(string(rag.group)))', sort(groups)), ...
    'Unexpected RAG1 groups for the protection model.');
Score = rag.relweight;
Condition = categorical(string(rag.group), ...
    ["uninfected_b6", "st175_b6", "st175_rag1ko"]);
DayGroup = categorical(string(rag.day));
Animal = categorical(string(rag.group) + "_" + string(rag.mouse_id));
modelData = table(Score, Condition, DayGroup, Animal);
model = fitlme(modelData, ...
    'Score ~ Condition + (1|DayGroup) + (1|Animal)', ...
    'FitMethod', 'REML');
coef = model.Coefficients;
wildType = coef(string(coef.Name) == "Condition_st175_b6", :);
ragKo = coef(string(coef.Name) == "Condition_st175_rag1ko", :);
assert(height(wildType) == 1 && height(ragKo) == 1, ...
    'RAG1 protection model is missing an expected condition coefficient.');

Comparison = ["ST1-75-colonized WT versus previously uncolonized WT"; ...
    "ST1-75-colonized RAG1 KO versus previously uncolonized WT"];
Effect_Label = ["ST1-75-colonized WT"; "ST1-75-colonized RAG1 KO"];
Effect = [wildType.Estimate; ragKo.Estimate];
CI_Lower = [wildType.Lower; ragKo.Lower];
CI_Upper = [wildType.Upper; ragKo.Upper];
P_Value = [wildType.pValue; ragKo.pValue];
N_Animals = [20; 20];
Model = repmat("relative weight ~ condition + (1|day) + (1|animal)", 2, 1);
effects = table(Effect_Label, Comparison, Effect, CI_Lower, CI_Upper, ...
    P_Value, N_Animals, Model);
end


function make_figure5(outputFile, rag, ragEffects, ...
    microbiome, ~, ~, adaptiveFile, innateFile)
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.2 0.2 14.5 8.8]);

axA = axes(fig,'Position',[0.055 0.59 0.34 0.34]);
plot_rag1(axA, rag);
title(axA, 'Weight after VPI10463 challenge', ...
    'FontSize', 11, 'FontWeight', 'bold');
panel_label(axA, 'A');

axB = axes(fig,'Position',[0.465 0.59 0.19 0.34]);
plot_rag1_protection_bars(axB, ragEffects);
title(axB, {'Protection in WT and', 'RAG1-deficient mice'}, ...
    'FontSize', 11, 'FontWeight', 'bold');
panel_label(axB, 'B');

axC = axes(fig,'Position',[0.055 0.14 0.40 0.30]);
plot_draft39_vector_flow(axC, 'adaptive', 'Adaptive immune-cell fractions');
panel_label(axC, 'C');

axD = axes(fig,'Position',[0.55 0.14 0.41 0.30]);
plot_draft39_vector_flow(axD, 'innate', 'Innate immune-cell fractions');
panel_label(axD, 'D');

axE = axes(fig,'Position',[0.735 0.59 0.245 0.34]);
plot_residual_microbiome(axE, microbiome);
title(axE, {'Residual microbiota', 'after colonization'}, ...
    'FontSize', 11, 'FontWeight', 'bold');
panel_label(axE, 'E');

export_draft39_vector_asset(fig, outputFile);
exportgraphics(fig, outputFile, 'Resolution', 300);
close(fig);
end


function plot_rag1_protection_bars(ax, effects)
assert(height(effects) == 2, 'Expected two RAG1 protection-model effects.');
colors = [0.08 0.38 0.74; 0.48 0.16 0.60];
hold(ax, 'on');
for i = 1:height(effects)
    bar(ax, i, effects.Effect(i), 0.58, 'FaceColor', colors(i, :), ...
        'EdgeColor', [0.18 0.18 0.18], 'LineWidth', 0.6);
    errorbar(ax, i, effects.Effect(i), ...
        effects.Effect(i) - effects.CI_Lower(i), ...
        effects.CI_Upper(i) - effects.Effect(i), 'k', ...
        'LineWidth', 1.4, 'CapSize', 8);
end
yline(ax, 0, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.0);
set(ax, 'XTick', 1:2, 'XTickLabel', {'WT', 'RAG1 KO'}, ...
    'XLim', [0.45 2.55], ...
    'YLim', [0 14], ...
    'FontName', 'Arial', 'FontSize', 9, 'Box', 'off', ...
    'TickDir', 'out', 'YGrid', 'on', 'GridAlpha', 0.08);
xlabel(ax, 'ST1-75-colonized genotype');
ylabel(ax, 'Protection effect (percentage points)');
end


function plot_rag1(ax, rag)
rag.group = string(rag.group);
groups = ["uninfected_b6", "st175_b6", "st175_rag1ko"];
labels = ["Previously uncolonized WT", "ST1-75-colonized WT", ...
    "ST1-75-colonized RAG1 KO"];
colors = [0.42 0.42 0.42; 0.08 0.38 0.74; 0.48 0.16 0.60];
hold(ax, 'on');
finalScheduledDay = max(rag.day);
for i = 1:numel(groups)
    sub = rag(rag.group == groups(i), :);
    mice = unique(sub.mouse_id);
    days = sort(unique(sub.day));
    pale = 0.79 * [1 1 1] + 0.21 * colors(i, :);
    for m = 1:numel(mice)
        one = sortrows(sub(sub.mouse_id == mice(m), :), 'day');
        observed = ~isnan(one.relweight);
        plot(ax, one.day(observed), one.relweight(observed), '-', ...
            'Color', pale, ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
        lastObserved = find(observed, 1, 'last');
        if groups(i) == "uninfected_b6" && ...
                one.day(lastObserved) < finalScheduledDay
            plot(ax, one.day(lastObserved), one.relweight(lastObserved), ...
                'x', 'Color', [0.42 0.42 0.42], 'MarkerSize', 7, ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    means = nan(numel(days), 1);
    for d = 1:numel(days)
        means(d) = mean(sub.relweight(sub.day == days(d)), 'omitnan');
    end
    plot(ax, days, means, '-o', 'Color', colors(i, :), ...
        'LineWidth', 2.6, 'MarkerFaceColor', colors(i, :), ...
        'MarkerSize', 4.5, 'DisplayName', labels(i));
end
plot(ax, nan, nan, 'x', 'Color', [0.42 0.42 0.42], ...
    'MarkerSize', 7, 'LineWidth', 1.5, ...
    'DisplayName', 'Death');
yline(ax, 100, ':', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 0.9, 'HandleVisibility', 'off');
xlabel(ax, 'Day after VPI10463 challenge');
ylabel(ax, 'Relative weight (%)');
xlim(ax, [0 7]);
ylim(ax, [60 110]);
legend(ax, 'Location', 'southwest', 'Box', 'off', 'FontSize', 8.5);
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off', ...
    'TickDir', 'out', 'XGrid', 'on', 'GridAlpha', 0.08);
end


function plot_residual_microbiome(ax, microbiome)
alpha = microbiome(microbiome.Analysis == "Residual_Shannon", :);
hold(ax, 'on');
errorbar(ax, alpha.Day - 0.06, alpha.ST1_75_Mean, alpha.ST1_75_SD, ...
    '-o', 'Color', [0.08 0.38 0.74], 'LineWidth', 2.0, ...
    'MarkerFaceColor', [0.08 0.38 0.74], 'MarkerSize', 4.5, ...
    'DisplayName', 'ST1-75');
errorbar(ax, alpha.Day + 0.06, alpha.ST1_12_Mean, alpha.ST1_12_SD, ...
    '-o', 'Color', [0.88 0.35 0.10], 'LineWidth', 2.0, ...
    'MarkerFaceColor', [0.88 0.35 0.10], 'MarkerSize', 4.5, ...
    'DisplayName', 'ST1-12');
for i = 1:height(alpha)
    top = max(alpha.ST1_75_Mean(i) + alpha.ST1_75_SD(i), ...
        alpha.ST1_12_Mean(i) + alpha.ST1_12_SD(i));
    text(ax, alpha.Day(i), top + 0.055, ...
        sprintf('q=%.2g', alpha.Q_Value(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end
text(ax, 0.98, 0.06, 'All Shannon q values ≥ 0.24', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'FontSize', 8.2, 'Color', [0.30 0.30 0.30]);
xlabel(ax, 'Day after colonization');
ylabel(ax, 'Residual Shannon index');
xlim(ax, [-0.35 7.35]);
ylim(ax, [0 2.12]);
legend(ax, 'Location', 'northwest', 'Orientation', 'horizontal', ...
    'Box', 'off', 'FontSize', 8.5);
set(ax, 'FontName', 'Arial', 'FontSize', 9.5, 'Box', 'off', ...
    'TickDir', 'out', 'XGrid', 'on', 'GridAlpha', 0.10);
end


function make_figureS3(outputFile, microbiome)
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.2 0.2 12.0 9.2]);

annotation(fig, 'textbox', [0.05 0.94 0.90 0.045], ...
    'String', 'Expanded residual-microbiota comparisons', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontName', 'Arial', 'FontSize', 15, 'FontWeight', 'bold');

% Shannon diversity is shown once in Figure 5E.
axA = axes(fig, 'Position', [0.15 0.58 0.75 0.30]);
plot_bray_curtis_summary(axA, microbiome);
panel_label(axA, 'A');
set_panel_label_position(axA);

axB = axes(fig, 'Position', [0.23 0.10 0.65 0.35]);
plot_genus_summary(axB, microbiome);
panel_label(axB, 'B');
set_panel_label_position(axB);

export_draft39_vector_asset(fig, outputFile);
exportgraphics(fig, outputFile, 'Resolution', 300);
close(fig);
end

function set_panel_label_position(ax)
for label=findall(ax,'Type','text')'
    if ischar(label.String) && ismember(label.String,{'A','B'})
        label.Position=[-0.02 1.10 0];
    end
end
end


function plot_bray_curtis_summary(ax, microbiome)
bray = microbiome(microbiome.Analysis == "Residual_BrayCurtis", :);
[~, order] = sort(bray.Day);
bray = bray(order, :);
r2 = str2double(extractBetween(bray.Feature, "R2 ", ", F"));
assert(all(isfinite(r2)) && height(bray) == 4, ...
    'Unable to recover the four verified Bray-Curtis R2 summaries.');

bars = bar(ax, 1:height(bray), r2, 0.62, ...
    'FaceColor', [0.30 0.60 0.67], 'EdgeColor', 'none');
bars.FaceAlpha = 0.90;
hold(ax, 'on');
for i = 1:height(bray)
    text(ax, i, r2(i) + 0.025, ...
        sprintf('p=%.3g\nq=%.2g', bray.P_Value(i), bray.Q_Value(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontName', 'Arial', 'FontSize', 8.2);
end
yline(ax, 0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);
set(ax, 'XTick', 1:height(bray), ...
    'XTickLabel', compose('Day %d', bray.Day), ...
    'FontName', 'Arial', 'FontSize', 9.5, 'Box', 'off', ...
    'TickDir', 'out', 'YGrid', 'on', 'GridAlpha', 0.10);
ylabel(ax, 'Bray-Curtis group effect (R^2)');
ylim(ax, [0 0.46]);
title(ax, 'Residual community-composition differences', ...
    'FontSize', 11.5, 'FontWeight', 'bold');
text(ax, 0.98, 0.04, 'Permutation tests; Benjamini-Hochberg-adjusted q values', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'FontName', 'Arial', 'FontSize', 7.8, 'Color', [0.30 0.30 0.30]);
end


function plot_genus_summary(ax, microbiome)
genera = microbiome(microbiome.Analysis == "Residual_Top_Genus", :);
[~, order] = sort(abs(genera.Effect_Size), 'ascend');
genera = genera(order, :);
effects = genera.Effect_Size;
colors = repmat([0.08 0.38 0.74], height(genera), 1);
colors(effects < 0, :) = repmat([0.88 0.35 0.10], ...
    sum(effects < 0), 1);

hold(ax, 'on');
for i = 1:height(genera)
    barh(ax, i, effects(i), 0.66, 'FaceColor', colors(i, :), ...
        'EdgeColor', 'none');
    text(ax, 43.5, i, sprintf('q=%.2g', genera.Q_Value(i)), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', 'FontSize', 7.6, ...
        'Color', [0.25 0.25 0.25]);
end
xline(ax, 0, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);
labels = replace(genera.Feature, "unclassified_Muribaculaceae_ASV_4", ...
    "Muribaculaceae ASV 4");
labels = replace(labels, "ClostridiumSensuStricto1", ...
    "Clostridium sensu stricto 1");
labels = replace(labels, "LachnospiraceaeNK4A136Group", ...
    "Lachnospiraceae NK4A136");
set(ax, 'YTick', 1:height(genera), 'YTickLabel', labels, ...
    'FontName', 'Arial', 'FontSize', 8.2, 'Box', 'off', ...
    'TickDir', 'out', 'XGrid', 'on', 'GridAlpha', 0.10);
xlim(ax, [-37 45]);
xlabel(ax, {'Mean relative-abundance difference', ...
    '(ST1-75 - ST1-12, percentage points)'});
title(ax, 'Day-1 genus-level contrasts', ...
    'FontSize', 11.5, 'FontWeight', 'bold');
end


function panel_label(ax, label)
text(ax, -0.06, 1.085, label, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 15, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', 'Color', [0.05 0.05 0.05], ...
    'BackgroundColor', 'none');
end
