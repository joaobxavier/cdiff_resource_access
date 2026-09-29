function growth = analyze_biolog_plate(plateFile, biologNamesFile, opts)
%   Usage:
%       g = analyze_biolog_plate('ST1-20_PM1.xlsx', 'Biolog_names_sorted.xlsx');
%       g = analyze_biolog_plate('ST1-75_PM1.xlsx', 'Biolog_names_sorted.xlsx', ...
%                                struct('Plot',true,'Write',true, ...
%                                       'OutDir','results/biolog'));

  arguments
        plateFile        (1,1) string
        biologNamesFile  (1,1) string
        opts.Plot        (1,1) logical = false
        opts.Write       (1,1) logical = false
        opts.OutDir      (1,1) string  = pwd
        opts.TimeCutoff  (1,1) double  = 24
        opts.Alpha       (1,1) double  = 0.05
        opts.NumTests    (1,1) double  = 95
    end
 
    refLevel   = "Negative Control";
    strainName = parseStrainName(plateFile);
 
    % Load and reshape
    data = readtable(plateFile, 'VariableNamingRule', 'preserve');
    data = stack(data, 2:97, ...
        'NewDataVariableName', 'Absorbance', ...
        'IndexVariableName',   'Well');
    data.Properties.VariableNames{1} = 'Time';
    data = data(data.Time < opts.TimeCutoff, :);
 
    % Attach metabolite names
    biolog      = readtable(biologNamesFile);
    biolog.Well = categorical(biolog.Well);
    data.Well   = categorical(data.Well);
    data        = join(data, biolog, 'Keys', 'Well');
    data        = rmmissing(data);
 
    data.Names = categorical(data.Names);
    data.Names = reordercats(data.Names, ...
        [refLevel, setdiff(categories(data.Names), refLevel)']);
    data = sortrows(data, {'Names', 'Well', 'Time'});
 
    met = categories(removecats(data.Names));
    assert(numel(met) >= 2, 'Fewer than two metabolite groups found.');
 
    [g, ~]  = findgroups(data.Well);
    t0       = splitapply(@(a,t) a(find(t==min(t),1)), ...
                          data.Absorbance, data.Time, g);
    data.AbsorbanceminusT0 = max(data.Absorbance - t0(g), 0);
 
    ODmodel = fitlme(data, 'AbsorbanceminusT0 ~ Time*Names + (Time|Well)');
 
    mdlResults = table( ...
        strrep(string(ODmodel.CoefficientNames(:)), 'Metabolites_', ''), ...
        ODmodel.Coefficients.Estimate, ...
        ODmodel.Coefficients.Upper, ...
        ODmodel.Coefficients.Lower, ...
        ODmodel.Coefficients.pValue, ...
        'VariableNames', {'Metabolites','impact','upper','lower','pvalue'});
 
    % Bonferroni correction
    thresh = opts.Alpha / opts.NumTests;
    sig    = mdlResults(mdlResults.pvalue <= thresh & mdlResults.impact > 0, :);
    sig(sig.Metabolites == "(Intercept)", :) = [];      
    sig.ParsedName = extractAfter(sig.Metabolites, '_');
 
    % Binary growth table 
    growth           = table();
    growth.Metabolites = string(met);
    growth.Growth      = double(ismember(growth.Metabolites, sig.ParsedName));
    fprintf('%s: %d / %d metabolites support growth.\n', ...
        strainName, sum(growth.Growth), height(growth));
 
    %  Control QC
    checkControl(growth, "a-D-Glucose",      1, strainName);
    checkControl(growth, "Negative Control", 0, strainName);
    if opts.Plot
        plotControls(data, "a-D-Glucose", "Negative Control", strainName);
    end
 
    % Visualize
    if opts.Plot
        plotGrid(data, met, refLevel, strainName);
    end
 
    % write 
    if opts.Write
        if ~exist(opts.OutDir, 'dir'); mkdir(opts.OutDir); end
        outFile = fullfile(opts.OutDir, strainName + ".xlsx");
        writetable(growth, outFile);
        fprintf('Wrote %s\n', outFile);
    end
end
 
function checkControl(growth, name, expected, strainName)
 
    idx = growth.Metabolites == name;
    if ~any(idx)
        warning('%s: control "%s" not found in results.', strainName, name);
        return;
    end
    actual = growth.Growth(idx);
    if actual == expected
        fprintf('  QC ok:   %-16s = %d (expected %d)\n', name, actual, expected);
    else
        warning('%s: QC FAIL: %s = %d, expected %d.', ...
            strainName, name, actual, expected);
    end
    
end
% 
function plotControls(data, posName, negName, strainName)
% Plot positive control (glucose) vs negative control curves side by side.
    figure('Name', strainName + " controls", 'Position', [150 150 700 300]);
 
    names = [posName, negName];
    cols  = {'r', 'b'};
    for k = 1:2
        subplot(1, 2, k);
        m = data.Names == names(k);
        if ~any(m)
            title(names(k) + " (not found)", 'Interpreter', 'none');
            continue;
        end
        t  = data.Time(m);
        od = data.AbsorbanceminusT0(m);
        [t, s] = sort(t); od = od(s);
        plot(t, od, cols{k}, 'LineWidth', 2);
        ylim([0 1]); xlabel('Time'); ylabel('OD - T0');
        title(names(k), 'Interpreter', 'none');
    end
    sgtitle(strainName + " controls", 'Interpreter', 'none');
end
 
function name = parseStrainName(plateFile)
    [~, base] = fileparts(plateFile);
    tok = regexp(base, '^(.*?)(_PM\d+)?$', 'tokens', 'once');
    name = string(tok{1});
end
 
function plotGrid(data, met, refLevel, strainName)
    ncMask = data.Names == refLevel;
    tNC    = data.Time(ncMask);
    odNC   = data.AbsorbanceminusT0(ncMask);
    [tNC, k] = sort(tNC); odNC = odNC(k);
 
    figure('Name', strainName, 'Position', [100 100 1400 900]);
    for i = 1:numel(met)
        subplot(8, 12, i);
        m  = data.Names == met{i};
        t  = data.Time(m);
        od = data.AbsorbanceminusT0(m);
        [t, k] = sort(t); od = od(k);
 
        plot(t, od, 'r', 'LineWidth', 2); hold on;
        plot(tNC, odNC, 'b', 'LineWidth', 2);
        title(met{i}, 'Interpreter', 'none');
    end
    set(findall(gcf, 'Type', 'axes'), 'YLim', [0 1]);
    sgtitle(strainName, 'Interpreter', 'none');
end