function assess_pm1_equal_replication(base, classFile, out)
% Package defaults; explicit paths also support isolated verification.
if nargin==0
    base=fullfile(package_root(),'results','tables');
    classFile=package_input('moas.xlsx');
    out=fullfile(base,'pm1_equal_replication');
end
% Exhaustive one-plate sensitivity; never modifies authoritative inputs.
if ~isfolder(out), mkdir(out); end
plateDir = fullfile(base,'biolog');
pm = readtable(fullfile(plateDir,'Biolog_growth_matrix.csv'),'TextType','string','VariableNamingRule','preserve');
pm = pm(pm.Metabolites ~= "Negative Control",:);
cl = readtable(classFile,'TextType','string','VariableNamingRule','preserve');
[ok,ix] = ismember(pm.Metabolites,cl.Chemical); assert(all(ok));
aa = regexprep(cl.MoA(ix),'^C-Source, ','') == "amino acid";
assert(sum(aa)==16 && height(pm)==95);
names = replace(string(pm.Properties.VariableNames(2:end)),"_","-");
unionCalls = double(pm{:,2:end}); assert(sum(unionCalls,'all')==700);
files = dir(fullfile(plateDir,'*_calls.csv')); assert(numel(files)==30);
calls = zeros(95,30); strain = strings(30,1); source = strings(30,1);
for k=1:30
    source(k)=string(files(k).name);
    strain(k)=string(regexp(files(k).name,'^(ST1-\d+|VPI)','match','once'));
    t=readtable(fullfile(plateDir,files(k).name),'TextType','string');
    [ok,ix]=ismember(pm.Metabolites,t.Metabolites); assert(all(ok));
    calls(:,k)=t.Growth(ix);
end
assert(all(ismember(calls(:),[0 1])));
for j=1:numel(names)
    assert(isequal(max(calls(:,strain==names(j)),[],2),unionCalls(:,j)), ...
        'Plate union does not match authoritative matrix.');
end
plateCounts=table(strain,source,sum(calls,1)',sum(calls(aa,:),1)', ...
    'VariableNames',{'Strain','PlateFile','TotalCalls','AminoAcidCalls'});
writetable(plateCounts,fullfile(out,'plate_call_counts.csv'));
s=readtable(fullfile(base,'mouse_weight_survival_reanalysis','primary_terminal_event_zero_scores.csv'),'TextType','string');
[ok,col]=ismember(s.Strain,names); assert(all(ok) && height(s)==21);
vpi=find(names=="VPI");
y=s.Protection_Effect; d=s.Mono_Disease_Effect;
baseline=metrics(sum(unionCalls(aa,col),1)'-sum(unionCalls(aa,vpi)),y,d);
assert(abs(baseline(1)-0.702371354675232)<1e-12);
assert(abs(baseline(3)-0.470057708101081)<1e-12);
metricNames={'RawRho','RawP','DiseaseAdjustedRankRho','DiseaseAdjustedRankP', ...
    'AminoCoefficient','DiseaseOnlyAdjustedR2','DiseasePlusAminoAdjustedR2','AdjustedRankP_ControlDF'};
writetable(array2table(baseline,'VariableNames',metricNames),fullfile(out,'pooled_baseline.csv'));
rep=["ST1-12","ST1-68","ST1-75","VPI"];
indices=cell(1,4);
for j=1:4, indices{j}=find(strain==rep(j)); assert(numel(indices{j})==3); end
[a,b,c,e]=ndgrid(1:3,1:3,1:3,1:3); choices=[a(:),b(:),c(:),e(:)];
result=zeros(81,numel(metricNames)); predictors=zeros(81,21);
for k=1:81
    selected=zeros(1,numel(names));
    for j=1:numel(names), selected(j)=find(strain==names(j),1); end
    for j=1:4, selected(names==rep(j))=indices{j}(choices(k,j)); end
    x=sum(calls(aa,selected(col)),1)'-sum(calls(aa,selected(vpi)));
    predictors(k,:)=x'; result(k,:)=metrics(x,y,d);
end
t=[array2table(choices,'VariableNames',{'ST1_12_Plate','ST1_68_Plate','ST1_75_Plate','VPI_Plate'}), ...
    array2table(result,'VariableNames',metricNames)];
writetable(t,fullfile(out,'all_81_plate_selections.csv'));
writetable(array2table(predictors,'VariableNames',cellstr(replace(s.Strain,'-','_'))), ...
    fullfile(out,'all_81_relative_amino_predictors.csv'));
summary=table(string(metricNames)',baseline',min(result,[],1)',median(result,1)',max(result,[],1)', ...
    'VariableNames',{'Metric','Pooled','Minimum','Median','Maximum'});
writetable(summary,fullfile(out,'sensitivity_summary.csv'));
disp(plateCounts(ismember(strain,[rep,"ST1-6"]),:)); disp(summary);
fprintf('Raw p<0.05: %d/81; residual-correlation diagnostic p<0.05: %d/81\n',nnz(result(:,2)<.05),nnz(result(:,4)<.05));
fprintf('Adjusted p with one control degrees of freedom <0.05: %d/81\n',nnz(result(:,8)<.05));
fprintf('Distinct rank predictor patterns: %d\n',size(unique(tiedrank(predictors')','rows'),1));
end

function z=metrics(x,y,d)
[r,p]=corr(x,y,'Type','Spearman');
n=numel(y); design=[ones(n,1),tiedrank(d)];
R=eye(n)-design*pinv(design);
[ar,ap]=corr(R*tiedrank(x),R*tiedrank(y)); % rho matches Figure 4E; ap is an unadjusted residual-correlation diagnostic.
m0=fitlm(d,y); m1=fitlm([d,x],y);
controlP=2*tcdf(-abs(ar)*sqrt((n-3)/(1-ar^2)),n-3);
z=[r,p,ar,ap,m1.Coefficients.Estimate(3),m0.Rsquared.Adjusted,m1.Rsquared.Adjusted,controlP];
end
