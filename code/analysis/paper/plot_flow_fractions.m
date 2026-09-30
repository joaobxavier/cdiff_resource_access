function plot_flow_fractions(ax, kind, panelTitle)
% Rebuild the previously embedded flow-cytometry bars from exact gated data.
data = load(fullfile(getappdata(0, 'CdiffVectorOutput'), ['flow_',kind,'.mat']));
groups = {data.UI_data, data.Avirulent_data};
colors = [0.4 0.6 0.9; 0.9 0.5 0.3];
hold(ax, 'on');
n = size(groups{1}, 2);
bars = gobjects(2, 1);
for i = 1:2
    x = (1:n) - 0.35 + (i-0.5)*0.35;
    bars(i) = bar(ax, x, mean(groups{i},1), 0.35, ...
        'FaceColor', colors(i,:), 'FaceAlpha', 0.85, 'EdgeColor', 'none');
    errorbar(ax, x, mean(groups{i},1), std(groups{i},0,1), ...
        'k.', 'LineWidth', 1, 'CapSize', 4, 'HandleVisibility', 'off');
    for j = 1:n
        scatter(ax, repmat(x(j),3,1), groups{i}(:,j), 15, 'k', ...
            'filled', 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
    end
end
set(ax, 'XTick', 1:n, 'XTickLabel', data.prettyTypeNames, ...
    'XTickLabelRotation', 45, 'TickLabelInterpreter', 'tex', ...
    'FontName', 'Arial', 'FontSize', 9, 'Box', 'on', ...
    'YGrid', 'on', 'XGrid', 'on', 'GridAlpha', 0.12);
if strcmp(kind, 'adaptive'); ylim(ax,[0 0.6]); else; ylim(ax,[0 0.4]); end
xlabel(ax, [upper(kind(1)),kind(2:end),' immune cell type']);
ylabel(ax, 'Fraction of retained events');
legend(ax, bars, {'Previously uncolonized', 'ST1-75-colonized'}, ...
    'Location','northeast','FontSize',8);
title(ax, panelTitle, 'FontSize',11,'FontWeight','bold');
end
