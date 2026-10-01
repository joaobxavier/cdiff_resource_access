function test_code_consistency()
% Check statistical definitions and current figure composition independently.
root=package_root(); folder=fullfile(root,'results','tables','main');
t=readtable(fullfile(folder,'figure4_phylogeny_aligned_traits.csv'),'TextType','string');
calls=readtable(fullfile(folder,'figure4_pm1_call_matrix_plotted_order.csv'),'TextType','string');
c=readtable(fullfile(folder,'figure4_class_associations.csv'),'TextType','string');
rawP=nan(height(c),1); adjustedP=rawP; adjustedRho=rawP;
matrix=zeros(height(calls),height(t));
for k=1:height(t), matrix(:,k)=calls.(replace(t.Strain(k),'-','_')); end
for j=1:height(c)
    x=sum(matrix(calls.Chemical_Class==c.Category(j),:),1)';
    [rawRho,rawP(j)]=corr(x,t.Protection_Effect,'Type','Spearman');
    [adjustedRho(j),adjustedP(j)]=partialcorr(x,t.Protection_Effect, ...
        t.Total_Breadth,'Type','Spearman');
    assert(abs(rawRho-c.Spearman_Rho(j))<1e-12);
end
actualP=c.Protection_Partial_Rank_Adjusted_For_Total_Breadth_P_Value;
assert(max(abs(adjustedRho-c.Protection_Partial_Rank_Adjusted_For_Total_Breadth_Rho))<1e-12);
assert(max(abs(actualP-adjustedP))<1e-12,'Class-adjusted p-values must use one-control degrees of freedom.');
assert(max(abs(rawP-c.P_Value))<1e-12);
[sorted,ix]=sort(adjustedP); q=sorted.*height(c)./(1:height(c))';
q=flipud(cummin(flipud(q))); restored=nan(size(q)); restored(ix)=min(q,1);
assert(max(abs(restored-c.Protection_Partial_Rank_Adjusted_For_Total_Breadth_Q_Value))<1e-12);
out=fullfile(root,'results','verification'); if ~isfolder(out),mkdir(out);end
writetable(table(c.Category,adjustedRho,adjustedP,actualP,restored, ...
    repmat(height(t)-3,height(c),1),'VariableNames', ...
    {'Category','NativePartialRho','NativePartialP','OutputPartialP','BH_Q','DF'}), ...
    fullfile(out,'class_partial_rank_native_checks.csv'));
% The 81-plate file retains the plain residual-correlation p only as a
% diagnostic; the separately named control-df column is the inferential test.
p=readtable(fullfile(root,'results','tables','pm1_equal_replication','all_81_plate_selections.csv'));
df=height(t)-3;
definition=2*tcdf(-abs(p.DiseaseAdjustedRankRho).*sqrt(df./(1-p.DiseaseAdjustedRankRho.^2)),df);
assert(max(abs(definition-p.AdjustedRankP_ControlDF))<1e-12);
% Check saved artwork, not just filenames, for the new Figure 4C/F and S2D.
for name=["figure4_phylogeny_resource_breadth","figureS2_pm1_detail_and_generalization_limits"]
    figureFile=fullfile(root,'results','vectors',name+'.fig');
    if ~isfile(figureFile)
        figureFile=fullfile(root,'results','figures','editable_vectors',name+'.fig');
    end
    f=openfig(figureFile,'invisible');
    cleaner=onCleanup(@()close(f));
    ax=findall(f,'Type','axes'); titles=strings(numel(ax),1);
    for j=1:numel(ax),titles(j)=join(string(ax(j).Title.String),' ');end
    if startsWith(name,'figure4')
        assert(any(contains(titles,'Overlap with VPI10463')) && ...
            any(contains(titles,'VPI nutrient overlap')) && ...
            any(contains(titles,'Additional resource access')));
        assert(~any(contains(titles,'Shared and strain-specific PM1 calls')));
    else
        assert(numel(ax)==4 && any(contains(titles,'Protection by resource class')));
        assert(~any(contains(titles,'What the association')));
    end
    clear cleaner
end
fprintf('PASS: eight native class-adjusted tests (18 df), 81 plate p-values, Figure 4/S2 composition.\n');
end
