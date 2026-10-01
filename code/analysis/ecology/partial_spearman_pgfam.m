function result = partial_spearman_pgfam(x, y, controls)
%PARTIAL_SPEARMAN_PGFAM Average-rank residual correlation with correct df.
% Invalid designs are exposed rather than silently dropping rows/covariates.
% The status also permits explicitly counted degenerate bootstrap draws.

if nargin < 3
    controls = zeros(numel(x), 0);
end
x = x(:); y = y(:);
n = numel(x);
assert(numel(y) == n && size(controls, 1) == n);
k = size(controls, 2);
df = n - k - 2;
result = struct('Rho', NaN, 'PValue', NaN, 'TStatistic', NaN, ...
    'N', n, 'Covariates', k, 'DF', df, 'Valid', false, 'Reason', "");
if df <= 0 || any(~isfinite([x; y; controls(:)]))
    result.Reason = "Insufficient degrees of freedom or nonfinite input";
    return
end
ranks = tiedrank([x, y, controls]);
design = [ones(n, 1), ranks(:, 3:end)];
if rank(design) ~= k + 1
    result.Reason = "Rank-deficient covariate design";
    return
end
residuals = ranks(:, 1:2) - design * (design \ ranks(:, 1:2));
residuals = residuals - mean(residuals, 1);
norms = sqrt(sum(residuals.^2, 1));
if any(norms < 1e-10)
    result.Reason = "Constant predictor or response after adjustment";
    return
end
rho = dot(residuals(:, 1), residuals(:, 2)) / prod(norms);
rho = max(-1, min(1, rho));
if abs(rho) == 1
    t = sign(rho) * Inf;
else
    t = rho * sqrt(df / (1 - rho^2));
end
result.Rho = rho;
result.TStatistic = t;
result.PValue = 2 * tcdf(-abs(t), df);
result.Valid = true;
end
