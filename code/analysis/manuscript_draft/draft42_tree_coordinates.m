function data=draft42_tree_coordinates(analysisFile)
% Draw the supported tree with unsupported internal edges collapsed to zero.
% Original fitted branch lengths remain in the analysis tree for distances.
tree=phytreeread(analysisFile);
names=string(get(tree,'LeafNames'));tree=reroot(tree,find(names=="VPI"));
names=string(get(tree,'LeafNames'));tree=prune(tree,find(names=="VPI"));
names=string(get(tree,'NodeNames'));n=get(tree,'NumLeaves');
dist=get(tree,'Distances');collapsed=false(size(dist));
for i=n+1:numel(names)-1
    support=sscanf(names(i),'%f/%f');
    if numel(support)~=2||support(1)<80||support(2)<95
        dist(i)=0;collapsed(i)=true;
    end
end
tree=phytree(get(tree,'Pointers'),dist,cellstr(names));
data.CollapsedInternalEdges=sum(collapsed);
data.DisplayTree=tree;
f=figure('Visible','off','Color','w');h=plot(tree,'Type','square','Orientation','top');
data.BranchX=arrayfun(@(p)p.XData,h.BranchLines,'UniformOutput',false);
data.BranchY=arrayfun(@(p)p.YData,h.BranchLines,'UniformOutput',false);
x=arrayfun(@(p)p.XData,h.LeafDots,'UniformOutput',false);
y=arrayfun(@(p)p.YData,h.LeafDots,'UniformOutput',false);
data.LeafX=horzcat(x{:});data.LeafY=horzcat(y{:});
tipNames=string({h.terminalNodeLabels.String})';positions=vertcat(h.terminalNodeLabels.Position);
data.SourceYDir=h.axes.YDir;close(f);
[data.TipX,order]=sort(positions(:,1));
data.TipOrder=replace(tipNames(order),'.','-');
branchY=horzcat(data.BranchY{:});span=range(branchY);
data.SourceYLim=[min(branchY)-.08*span,max(branchY)+.28*span];
data.LabelBaselineY=data.SourceYLim(2)-.03*span;
assert(numel(data.TipOrder)==21 && numel(unique(data.TipOrder))==21);
end
