function summarize_final_tree_tests
% Recalculate all tree-dependent statistics; retain non-tree analyses unchanged.
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
project=fileparts(repo);
root=fullfile(package_root(),'results','phylogeny','tree_tests');
traits=readtable(fullfile(project,'results', 'tables', 'main','figure4_phylogeny_aligned_traits.csv'),'TextType','string');
tests=["masked_new22","unmasked_new22","masked_st1_same_columns","masked_new22_common","masked_old22_common"];
rows=cell(0,6);splitResults=cell(numel(tests),1);
for k=1:numel(tests)
    path=fullfile(root,tests(k)+'.treefile');
    tree=readUniqueTree(char(path));
    names=string(get(tree,'LeafNames'));
    isST1=startsWith(names,"ST1.");
    d=squareform(pdist(tree,'Nodes','leaves'));d=d(isST1,isST1);
    candidates=replace(names(isST1),'.','-');
    [candidates,order]=sort(candidates);d=d(order,order);
    [ok,ix]=ismember(candidates,traits.Strain);assert(all(ok));t=traits(ix,:);
    upper=triu(true(21),1);
    values=[t.Total_Breadth,t.Protection_Effect,t.Mono_Disease_Effect,t.Amino_Acid_Relative_Breadth_Advantage];
    labels=["Total PM1 call breadth","Co-infection protection effect","Mono-colonization disease effect","Amino-acid relative access"];
    rng(1);
    for j=1:4
        if j==4,rng(29);end
        val=values(:,j);delta=abs(val-val');
        rho=corr(d(upper),delta(upper),'Type','Spearman');null=zeros(10000,1);
        for r=1:10000
            p=val(randperm(21));pd=abs(p-p');
            null(r)=corr(d(upper),pd(upper),'Type','Spearman');
        end
        pvalue=(1+sum(abs(null)>=abs(rho)))/10001;
        focal=d(candidates=="ST1-75",candidates=="ST1-68");
        rows(end+1,:)={tests(k),labels(j),rho,pvalue,focal,10000}; %#ok<AGROW>
    end
    if any(~isST1),tree=prune(tree,find(~isST1));end
    splitResults{k}=splitTable(tree);
    writetable(splitResults{k},fullfile(root,tests(k)+'_splits.csv'));
end
result=cell2table(rows,'VariableNames',{'Tree','Trait','Rho','PValue','FocalPatristicDistance','Permutations'});
writetable(result,fullfile(root,'tree_trait_sensitivity.csv'));disp(result);
primary=splitResults{1};StrongSplits=zeros(5,1);SharedStrongWithPrimary=StrongSplits;
for k=1:5
    s=splitResults{k};strong=s.SH>=80&s.UFB>=95;
    StrongSplits(k)=sum(strong);
    SharedStrongWithPrimary(k)=sum(ismember(s.Split(strong),primary.Split(primary.SH>=80&primary.UFB>=95)));
end
writetable(table(tests',StrongSplits,SharedStrongWithPrimary,'VariableNames',{'Tree','StrongSplits','SharedStrongWithPrimary'}),fullfile(root,'supported_split_comparison.csv'));
end

function t=splitTable(tree)
names=string(get(tree,'LeafNames')); n=numel(names);
p=get(tree,'Pointers'); labels=string(get(tree,'NodeNames'));
desc=false(2*n-1,n);desc(1:n,:)=eye(n)>0;
keys=strings(n-1,1);sh=nan(n-1,1);uf=sh;valid=false(n-1,1);
for k=1:n-1
    node=n+k;
    desc(node,:)=desc(p(k,1),:) | desc(p(k,2),:);
    mask=desc(node,:);
    if sum(mask)<2 || sum(~mask)<2, continue; end
    % Canonical unrooted split: the side without the alphabetically first tip.
    [~,order]=sort(names);
    if mask(order(1)),mask=~mask;end
    keys(k)=join(sort(names(mask)),',');valid(k)=true;
    values=sscanf(labels(node),'%f/%f');
    if numel(values)==2,sh(k)=values(1);uf(k)=values(2);end
end
t=table(keys(valid),sh(valid),uf(valid),'VariableNames',{'Split','SH','UFB'});
[~,keep]=unique(t.Split);t=t(keep,:);
end

function tree=readUniqueTree(path)
% MATLAB requires unique branch names, unlike Newick support labels.
text=fileread(path);
[starts,ends,~,~,tokens]=regexp(text,'\)([0-9.]+/[0-9.]+):');
for k=numel(starts):-1:1
    text=[text(1:starts(k)-1) sprintf(')%s_node%d:',tokens{k}{1},k) text(ends(k)+1:end)];
end
text=strrep(text,'ST1_','ST1.'); text=strrep(text,'VPI10463','VPI');
temporary=strrep(path,'.treefile','.analysis.nwk');
fid=fopen(temporary,'w');assert(fid>=0);fprintf(fid,'%s',text);fclose(fid);
tree=phytreeread(temporary);
end

