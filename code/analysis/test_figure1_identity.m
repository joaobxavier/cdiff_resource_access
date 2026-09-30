function test_figure1_identity()
% Regression check for the Figure 1 mouse identities and numerical results.
root = package_root();
validate_inputs();
manifest = readtable(fullfile(root,'expected','regression_manifest.tsv'), ...
    'FileType','text','Delimiter','\t','TextType','string');
keep = contains(manifest.Actual, "figure1_") | ...
    contains(manifest.Actual, "coinfection_survival_comparisons");
manifest = manifest(keep,:);
rows = cell(height(manifest),4);
for i = 1:height(manifest)
    a = readtable(fullfile(root,manifest.Actual(i)), ...
        'Delimiter',',','TextType','string','VariableNamingRule','preserve');
    e = readtable(fullfile(root,manifest.Expected(i)), ...
        'Delimiter',',','TextType','string','VariableNamingRule','preserve');
    assert(height(a)==height(e));
    count = 0; maximum = 0;
    for field = string(e.Properties.VariableNames)
        if isnumeric(e.(field))
            actual = a.(field); expected = e.(field);
            assert(isequaln(isnan(actual),isnan(expected)));
            use = isfinite(expected);
            delta = abs(actual(use)-expected(use));
            assert(all(delta <= 1e-7 + 1e-7*abs(expected(use))));
            count = count + sum(use);
            if ~isempty(delta); maximum = max(maximum,max(delta)); end
        end
    end
    rows(i,:) = {manifest.Actual(i),"PASS",count,maximum};
end
checks = cell2table(rows,'VariableNames', ...
    {'Actual','Status','NumericValues','MaximumDifference'});
out = fullfile(root,'results','verification','figure1_identity');
if ~isfolder(out); mkdir(out); end
writetable(checks,fullfile(out,'numerical_regression.csv'));
copyfile(fullfile(root,'results','tables','main', ...
    'figure1_ks65_mouse_identity.csv'),out);
disp(checks);
fprintf('Figure 1 mouse-identity and numerical regression checks passed.\n');
end
