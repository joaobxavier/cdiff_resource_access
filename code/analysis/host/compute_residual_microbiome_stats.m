function microbiome_tbl = compute_residual_microbiome_stats(base_dir)
tbl = readtable(project_data_file('processed', '16s_sequencing', 'tblAbund.xls'), ...
    'VariableNamingRule', 'preserve');
tbl.Initial_infection = string(tbl.Initial_infection);
tbl.Initial_antibiotics = string(tbl.Initial_antibiotics);
genus_cols = tbl.Properties.VariableNames(6:end);
residual_cols = genus_cols(~strcmp(genus_cols, 'Clostridioides'));
residual = double(tbl{:, residual_cols});
residual_total = sum(residual, 2);
residual_frac = residual ./ residual_total;
residual_frac(residual_total == 0, :) = 0;

days = [0 1 2 7];
Analysis = strings(0, 1);
Day = zeros(0, 1);
Feature = strings(0, 1);
ST1_75_Mean = zeros(0, 1);
ST1_12_Mean = zeros(0, 1);
ST1_75_SD = zeros(0, 1);
ST1_12_SD = zeros(0, 1);
Effect_Size = zeros(0, 1);
P_Value = zeros(0, 1);
N_ST1_75 = zeros(0, 1);
N_ST1_12 = zeros(0, 1);

for d = days
    idx = tbl.Day == d & tbl.Initial_antibiotics == "mnvc" & ...
        (tbl.Initial_infection == "ST1_75" | tbl.Initial_infection == "ST1_12");
    X = residual_frac(idx, :);
    groups = tbl.Initial_infection(idx);
    shannon = shannon_rows(X);
    idx75 = groups == "ST1_75";
    idx12 = groups == "ST1_12";

    [Analysis, Day, Feature, ST1_75_Mean, ST1_12_Mean, ST1_75_SD, ST1_12_SD, ...
        Effect_Size, P_Value, N_ST1_75, N_ST1_12] = append_stat_row( ...
        Analysis, Day, Feature, ST1_75_Mean, ST1_12_Mean, ST1_75_SD, ST1_12_SD, ...
        Effect_Size, P_Value, N_ST1_75, N_ST1_12, "Residual_Shannon", d, ...
        "Residual microbiota Shannon", shannon(idx75), shannon(idx12), "rank");

    if sum(idx75) >= 2 && sum(idx12) >= 2
        [pseudo_f, p_perm, r2] = permanova_two_group(X, groups, 999);
        [Analysis, Day, Feature, ST1_75_Mean, ST1_12_Mean, ST1_75_SD, ST1_12_SD, ...
            Effect_Size, P_Value, N_ST1_75, N_ST1_12] = append_stat_row( ...
            Analysis, Day, Feature, ST1_75_Mean, ST1_12_Mean, ST1_75_SD, ST1_12_SD, ...
            Effect_Size, P_Value, N_ST1_75, N_ST1_12, "Residual_BrayCurtis", d, ...
            sprintf('Residual Bray-Curtis R2 %.3f, F %.2f', r2, pseudo_f), ...
            nan(sum(idx75), 1), nan(sum(idx12), 1), p_perm);
    end
end

key_day = 1;
idx_key = tbl.Day == key_day & tbl.Initial_antibiotics == "mnvc" & ...
    (tbl.Initial_infection == "ST1_75" | tbl.Initial_infection == "ST1_12");
X_key = residual_frac(idx_key, :);
groups_key = tbl.Initial_infection(idx_key);
[~, top_ord] = sort(mean(X_key, 1, 'omitnan'), 'descend');
top_n = min(10, numel(top_ord));
for k = 1:top_n
    col = top_ord(k);
    idx75 = groups_key == "ST1_75";
    idx12 = groups_key == "ST1_12";
    vals75 = 100 * X_key(idx75, col);
    vals12 = 100 * X_key(idx12, col);
    [Analysis, Day, Feature, ST1_75_Mean, ST1_12_Mean, ST1_75_SD, ST1_12_SD, ...
        Effect_Size, P_Value, N_ST1_75, N_ST1_12] = append_stat_row( ...
        Analysis, Day, Feature, ST1_75_Mean, ST1_12_Mean, ST1_75_SD, ST1_12_SD, ...
        Effect_Size, P_Value, N_ST1_75, N_ST1_12, "Residual_Top_Genus", key_day, ...
        string(residual_cols{col}), vals75, vals12, "rank");
end

Q_Value = bh_fdr(P_Value);
microbiome_tbl = table(Analysis, Day, Feature, ST1_75_Mean, ST1_12_Mean, ...
    ST1_75_SD, ST1_12_SD, Effect_Size, P_Value, Q_Value, N_ST1_75, N_ST1_12);
end


function q = bh_fdr(p)
q = nan(size(p));
idx = find(~isnan(p));
if isempty(idx)
    return;
end
[p_sorted, ord] = sort(p(idx));
m = numel(p_sorted);
q_sorted = p_sorted(:) .* m ./ (1:m)';
q_sorted = flipud(cummin(flipud(q_sorted)));
q(idx(ord)) = min(q_sorted, 1);
end


function shannon = shannon_rows(X)
shannon = zeros(size(X, 1), 1);
for i = 1:size(X, 1)
    p = X(i, :);
    p = p(p > 0);
    shannon(i) = -sum(p .* log(p));
end
end


function [Analysis, Day, Feature, ST1_75_Mean, ST1_12_Mean, ST1_75_SD, ST1_12_SD, ...
    Effect_Size, P_Value, N_ST1_75, N_ST1_12] = append_stat_row( ...
    Analysis, Day, Feature, ST1_75_Mean, ST1_12_Mean, ST1_75_SD, ST1_12_SD, ...
    Effect_Size, P_Value, N_ST1_75, N_ST1_12, analysis, day, feature, vals75, vals12, test_spec)

vals75 = vals75(~isnan(vals75));
vals12 = vals12(~isnan(vals12));
mean75 = mean(vals75, 'omitnan');
mean12 = mean(vals12, 'omitnan');
sd75 = std(vals75, 0, 'omitnan');
sd12 = std(vals12, 0, 'omitnan');
effect = mean75 - mean12;

if isnumeric(test_spec)
    p = test_spec;
elseif numel(vals75) >= 2 && numel(vals12) >= 2
    p = ranksum(vals75, vals12);
else
    p = nan;
end

Analysis(end+1, 1) = analysis;
Day(end+1, 1) = day;
Feature(end+1, 1) = feature;
ST1_75_Mean(end+1, 1) = mean75;
ST1_12_Mean(end+1, 1) = mean12;
ST1_75_SD(end+1, 1) = sd75;
ST1_12_SD(end+1, 1) = sd12;
Effect_Size(end+1, 1) = effect;
P_Value(end+1, 1) = p;
N_ST1_75(end+1, 1) = numel(vals75);
N_ST1_12(end+1, 1) = numel(vals12);
end


function [pseudo_f, p_value, r2] = permanova_two_group(X, groups, n_perm)
valid = all(isfinite(X), 2) & (groups == "ST1_75" | groups == "ST1_12");
X = X(valid, :);
groups = groups(valid);
if numel(unique(groups)) < 2 || numel(groups) < 4
    pseudo_f = nan;
    p_value = nan;
    r2 = nan;
    return;
end

D = bray_curtis_matrix(X);
[pseudo_f, r2] = permanova_stat(D, groups);
rng(42);
extreme = 1;
for i = 1:n_perm
    perm_groups = groups(randperm(numel(groups)));
    [f_perm, ~] = permanova_stat(D, perm_groups);
    if f_perm >= pseudo_f
        extreme = extreme + 1;
    end
end
p_value = extreme / (n_perm + 1);
end


function [pseudo_f, r2] = permanova_stat(D, groups)
n = size(D, 1);
H = eye(n) - ones(n) / n;
G = -0.5 * H * (D .^ 2) * H;
group_levels = unique(groups);
X = zeros(n, numel(group_levels));
for i = 1:numel(group_levels)
    X(:, i) = groups == group_levels(i);
end
H_group = X / (X' * X) * X';
ss_total = trace(G);
ss_group = trace(H_group * G);
ss_resid = ss_total - ss_group;
df_group = numel(group_levels) - 1;
df_resid = n - numel(group_levels);
pseudo_f = (ss_group / df_group) / (ss_resid / df_resid);
r2 = ss_group / ss_total;
end


function D = bray_curtis_matrix(X)
n = size(X, 1);
D = zeros(n);
for i = 1:n
    for j = (i+1):n
        denom = sum(X(i, :) + X(j, :));
        if denom == 0
            d = 0;
        else
            d = sum(abs(X(i, :) - X(j, :))) / denom;
        end
        D(i, j) = d;
        D(j, i) = d;
    end
end
end


function x = numeric_column(x)
if iscell(x) || isstring(x)
    x = str2double(string(x));
else
    x = double(x);
end
end
