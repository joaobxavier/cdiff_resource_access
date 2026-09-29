function mapping = label_ks65_mice(axC, axD, workbook)
% Match explicit KS65 mouse IDs and relative weights before labeling C.
% The Biometrics sheet names mice in row 7; rows 8:10 contain relative
% weights. Rows 2:5 independently link the same IDs to raw body weights.
% No identity is inferred from plotted order or trajectory similarity.
cells = readcell(workbook, 'Sheet', 'Biometrics');
ids = cell2mat(cells(7, 2:6));
rawIds = cell2mat(cells(2, 2:6));
days = cell2mat(cells(8:10, 1));
relative = cell2mat(cells(8:10, 2:6));
raw = cell2mat(cells(3:5, 2:6));
assert(isequal(ids, [1 2 3 4 5]) && isequal(ids, rawIds));
assert(isequal(days, [0; 2; 3]));
assert(max(abs(relative - 100 .* raw ./ raw(1,:)), [], 'all') < 1e-10);
lines = findall(axC, 'Type', 'line');
beforeX = get(lines, 'XData');
beforeY = get(lines, 'YData');
labelY = [102 91.57 97];
Mouse = (3:5)';
WeightColumn = Mouse + 1;
for k = 1:3
    id = Mouse(k);
    col = ids == id;
    matches = false(size(lines));
    for j = 1:numel(lines)
        matches(j) = isequal(lines(j).XData(:), days) && ...
            numel(lines(j).YData) == 3 && ...
            max(abs(lines(j).YData(:) - relative(:,col))) < 1e-10;
    end
    assert(sum(matches) == 1, 'Mouse %d: weight match is not unique.', id);
    h = lines(matches);
    q = findall(axD, 'Type', 'line', 'DisplayName', sprintf('Mouse %d', id));
    assert(isscalar(q), 'Mouse %d: qPCR identity is unavailable.', id);
    h.Color = q.Color;
    % IDs beside the last observed point; label displacement changes no data.
    text(axC, 3.08, labelY(k), string(id), 'Color', q.Color, ...
        'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle', 'Tag', 'KS65MouseIdentity');
end
xlim(axC, [0 3.38]);
assert(isequaln(beforeX, get(lines, 'XData')) && ...
    isequaln(beforeY, get(lines, 'YData')), 'Observed weights changed.');
mapping = table(Mouse, WeightColumn, ...
    relative(1,3:5)', relative(2,3:5)', relative(3,3:5)', ...
    'VariableNames', {'Mouse','BiometricsColumn','Day0','Day2','Day3'});
end
