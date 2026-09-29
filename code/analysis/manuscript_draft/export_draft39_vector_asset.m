function export_draft39_vector_asset(fig, originalName, varargin)
% Force true vector export and retain a MATLAB figure for future editing.
out = getappdata(0, 'Draft39VectorOutput');
[~, name] = fileparts(originalName);
% Flat per-bar colors are rasterized by MATLAB even in vector mode.
for b = findall(fig, 'Type', 'bar')'
    if isequal(b.FaceColor, 'flat') && strcmp(b.BarLayout, 'grouped')
        ax = b.Parent; hold(ax,'on'); oldChildren = ax.Children;
        replacements = gobjects(numel(b.XData),1);
        x = b.XEndPoints; y = b.YData; colors = b.CData;
        width = b.BarWidth * min(diff(unique(b.XData)));
        for j = 1:numel(x)
            color = colors(min(j,size(colors,1)),:);
            xx = x(j) + width/2*[-1 1 1 -1];
            yy = [b.BaseValue b.BaseValue y(j) y(j)];
            if strcmp(b.Horizontal,'on'); [xx,yy] = deal(yy,xx); end
            replacements(j) = patch(ax,xx,yy,color,'EdgeColor',b.EdgeColor,'FaceAlpha',b.FaceAlpha);
        end
        index = find(oldChildren==b);
        ax.Children = [oldChildren(1:index-1); replacements; b; oldChildren(index+1:end)];
        delete(b);
    end
end
% Replace the colorbar gradient with editable colored rectangles.
for cb = findall(fig,'Type','colorbar')'
    pos = cb.Position; lim = cb.Limits; ticks = cb.Ticks;
    labels = cb.TickLabels; label = cb.Label.String; fs = cb.FontSize;
    parent = cb.Axes; parentPos = parent.Position;
    cmap = colormap(parent); delete(cb); parent.Position = parentPos;
    ax = axes(fig,'Position',pos,'XLim',[0 1],'YLim',lim, ...
        'XTick',[],'YTick',ticks,'YTickLabel',labels,'YAxisLocation','right', ...
        'FontSize',fs,'Box','on'); hold(ax,'on');
    edges = linspace(lim(1),lim(2),size(cmap,1)+1);
    for j=1:size(cmap,1)
        patch(ax,[0 1 1 0],edges([j j j+1 j+1]),cmap(j,:),'EdgeColor','none');
    end
    ylabel(ax,label);
end
% MATLAB imagesc represents the 95-cell chemical-class strip as an image.
% Convert its indexed cells to individual colored rectangles without changing
% the values, colormap, axes, or layout.
for im = findall(fig, 'Type', 'image')'
    values = im.CData;
    if ismatrix(values) && size(values,2) == 1 && size(values,1) == 95
        ax = im.Parent;
        cmap = colormap(ax);
        limits = ax.CLim;
        hold(ax,'on');
        for row = 1:95
            index = min(size(cmap,1), max(1, floor((double(values(row))-limits(1)) / diff(limits) * size(cmap,1))+1));
            patch(ax, [0.5 1.5 1.5 0.5], row+[-0.5 -0.5 0.5 0.5], ...
                cmap(index,:), 'EdgeColor','none');
        end
        delete(im);
    end
end
assert(isempty(findall(fig, 'Type', 'image')), ...
    'Raster image remains in %s; rebuild it before exporting.', name);
for obj = findall(fig,'Type','patch')'
    if isequal(obj.FaceColor,'flat') && size(obj.FaceVertexCData,1)==1 && size(obj.FaceVertexCData,2)==3
        obj.FaceColor = obj.FaceVertexCData;
    end
end
drawnow;
for t = findall(fig,'Type','text')'
    if ischar(t.String) && contains(t.String,'All Shannon q values')
        t.String = 'All Shannon q values >= 0.24';
    end
end
if strcmp(name,'figure2_ranked_mouse_screen')
    for ax = findall(fig,'Type','axes')'
        s = string(ax.YLabel.String);
        if isscalar(s) && contains(s,'score') && ax.Position(3)>.7
            ax.YLabel.String = replace(s," (+",string(newline)+"(+");
        end
    end
end
exportgraphics(fig, fullfile(out, [name, '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
print(fig, fullfile(out, [name, '.svg']), '-dsvg', '-vector');
savefig(fig, fullfile(out, [name, '.fig']));
fprintf('Exported %s\n', name);
end
