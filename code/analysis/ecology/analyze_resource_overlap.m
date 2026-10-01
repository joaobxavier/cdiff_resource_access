function associations = analyze_resource_overlap()
% Nutrient and PGFam overlap: complete strain-level associations and sensitivity.
root=package_root(); out=fullfile(root,'results','analyses','resource_overlap');
if ~isfolder(out), mkdir(out); end
pg=calculate_pgfam_overlap();
base=fullfile(root,'results','tables'); plateDir=fullfile(base,'biolog');
s=sortrows(readtable(fullfile(base,'mouse_weight_survival_reanalysis', ...
    'primary_terminal_event_zero_scores.csv'),'TextType','string'),'Strain');
pm=readtable(fullfile(plateDir,'Biolog_growth_matrix.csv'),'TextType','string','VariableNamingRule','preserve');
pm=pm(pm.Metabolites~="Negative Control",:);
cl=readtable(package_input('moas.xlsx'),'TextType','string','VariableNamingRule','preserve');
[ok,ix]=ismember(pm.Metabolites,cl.Chemical); assert(all(ok));
aa=regexprep(cl.MoA(ix),'^C-Source, ','')=="amino acid";
calls=double(pm{:,2:end}); names=replace(string(pm.Properties.VariableNames(2:end)),'_','-');
[ok,col]=ismember(s.Strain,names); assert(all(ok) && numel(unique(col))==21);
assert(isequal(s.Strain,pg.Strain) && height(s)==21);
assert(height(pm)==95 && sum(aa)==16 && sum(calls,'all')==700 && all(ismember(calls(:),[0 1])));
vpi=find(names=="VPI"); assert(isscalar(vpi));
[overlap,shared,den]=nutrient_overlap(calls,col,vpi,aa);
assert(isequal(den,[23 3 20]));
assert(all(ismember(shared(:,2),0:3)));
relativeAA=sum(calls(aa,col),1)'-sum(calls(aa,vpi));
y=s.Protection_Effect; disease=s.Mono_Disease_Effect;
inputs=table(s.Strain,pg.GenomeID,y,disease,pg.VPIOverlapFraction, ...
    shared(:,1),shared(:,2),shared(:,3),overlap(:,1),overlap(:,2),overlap(:,3), ...
    sum(calls(:,col),1)',sum(calls(:,col),1)'-shared(:,1),relativeAA, ...
    'VariableNames',{'Strain','GenomeID','ProtectionEffect','DiseaseEffect','PGFamOverlap', ...
    'SharedAll','SharedAmino','SharedNonAmino','OverlapAll','OverlapAmino','OverlapNonAmino', ...
    'TotalAccess','AdditionalAccess','RelativeAminoAccess'});
writetable(inputs,fullfile(out,'strain_aligned_inputs.csv'));
writetable(table(["All";"Amino acids";"Other substrates"],den', ...
    'VariableNames',{'ResourceClass','VPIPositiveCalls'}),fullfile(out,'overlap_denominators.csv'));
labels=["PGFamRaw";"PGFamDiseaseAdjusted";"AminoAccessDiseasePGFamAdjusted"; ...
    "NutrientAllRaw";"NutrientAllDiseaseAdjusted";"NutrientAminoRaw"; ...
    "NutrientAminoDiseaseAdjusted";"NutrientNonAminoRaw";"NutrientNonAminoDiseaseAdjusted"];
x={pg.VPIOverlapFraction,pg.VPIOverlapFraction,relativeAA, ...
    overlap(:,1),overlap(:,1),overlap(:,2),overlap(:,2),overlap(:,3),overlap(:,3)};
z={zeros(21,0),disease,[disease,pg.VPIOverlapFraction], ...
    zeros(21,0),disease,zeros(21,0),disease,zeros(21,0),disease};
rho=zeros(9,1); p=rho; df=rho; nr=rho; np=rho;
for j=1:9
    r=partial_spearman_pgfam(x{j},y,z{j}); assert(r.Valid);
    rho(j)=r.Rho; p(j)=r.PValue; df(j)=r.DF;
    if isempty(z{j}), [nr(j),np(j)]=corr(x{j},y,'Type','Spearman');
    else, [nr(j),np(j)]=partialcorr(x{j},y,z{j},'Type','Spearman'); end
end
assert(max(abs(rho-nr))<1e-12 && max(abs(p-np))<1e-12);
assert(isequal(df,[19;18;17;19;18;19;18;19;18]));
q=[bh_adjust(p(1:3));bh_adjust(p(4:9))];
originalRng=rng; cleanup=onCleanup(@()rng(originalRng)); %#ok<NASGU>
rng(20261001,'twister'); B=5000; boot=nan(B,9);
for b=1:B
    selected=randi(21,21,1);
    for j=1:9
        r=partial_spearman_pgfam(x{j}(selected),y(selected),z{j}(selected,:));
        if r.Valid, boot(b,j)=r.Rho; end
    end
end
limits=nan(9,2); valid=sum(isfinite(boot),1)';
for j=1:9, limits(j,:)=prctile(boot(isfinite(boot(:,j)),j),[2.5 97.5]); end
associations=table(labels,repmat(21,9,1),df,rho,p,q, ...
    [repmat("PGFam three-test family",3,1);repmat("Nutrient six-test family",6,1)], ...
    limits(:,1),limits(:,2),valid,B-valid, ...
    'VariableNames',{'Comparison','N','DF','Rho','PValue','BH_Q','CorrectionFamily', ...
    'Bootstrap_CI_Lower','Bootstrap_CI_Upper','ValidBootstrapDraws','InvalidBootstrapDraws'});
writetable(associations,fullfile(out,'rank_associations.csv'));
writetable(array2table(boot,'VariableNames',cellstr(labels)),fullfile(out,'bootstrap_coefficients.csv'));
writetable(table(labels,rho,nr,p,np,'VariableNames', ...
    {'Comparison','ImplementedRho','NativeRho','ImplementedP','NativeP'}),fullfile(out,'native_matlab_crosschecks.csv'));
removed=strings(189,1); comparison=removed; lr=nan(189,1); lp=lr; ldf=lr; lv=false(189,1);
for k=1:21
    keep=true(21,1); keep(k)=false;
    for j=1:9
        row=(k-1)*9+j; r=partial_spearman_pgfam(x{j}(keep),y(keep),z{j}(keep,:));
        removed(row)=s.Strain(k); comparison(row)=labels(j);
        lr(row)=r.Rho; lp(row)=r.PValue; ldf(row)=r.DF; lv(row)=r.Valid;
    end
end
writetable(table(removed,comparison,lr,lp,ldf,lv,'VariableNames', ...
    {'RemovedStrain','Comparison','Rho','PValue','DF','Valid'}),fullfile(out,'leave_one_strain_out.csv'));
baselineRaw=partial_spearman_pgfam(relativeAA,y,zeros(21,0));
baselineAdj=partial_spearman_pgfam(relativeAA,y,disease);
assert(abs(baselineRaw.Rho-.702371354675232)<1e-12);
assert(abs(baselineAdj.Rho-.470057708101081)<1e-12);
assert(abs(baselineAdj.PValue-.0364940843264458)<1e-12);
% Reference coefficients check source alignment; estimates are computed above.
assert(abs(rho(4)-.29103)<1e-5 && abs(rho(6)-.57934)<1e-5 && abs(rho(7)-.41453)<1e-5);
equal_plate_sensitivity(plateDir,pm,names,calls,col,vpi,aa,y,disease,s.Strain,out);
disp(associations);
end

function [fraction,shared,den]=nutrient_overlap(calls,col,vpi,aa)
classes={true(size(aa)),aa,~aa}; shared=zeros(numel(col),3); den=zeros(1,3);
for j=1:3
    mask=classes{j}; target=logical(calls(mask,vpi)); den(j)=sum(target);
    shared(:,j)=sum(logical(calls(mask,col)) & target,1)';
end
fraction=shared./den; % Zero denominators intentionally produce NaN.
end

function q=bh_adjust(p)
[sorted,order]=sort(p); adjusted=sorted.*numel(p)./(1:numel(p))';
for j=numel(p)-1:-1:1, adjusted(j)=min(adjusted(j),adjusted(j+1)); end
q=zeros(size(p)); q(order)=min(adjusted,1);
end

function equal_plate_sensitivity(plateDir,pm,names,unionCalls,col,vpi,aa,y,disease,strainNames,out)
files=dir(fullfile(plateDir,'*_calls.csv')); assert(numel(files)==30);
plateCalls=zeros(95,30); strains=strings(30,1); source=strains;
for k=1:30
    source(k)=string(files(k).name);
    strains(k)=string(regexp(files(k).name,'^(ST1-\d+|VPI)','match','once'));
    t=readtable(fullfile(plateDir,files(k).name),'TextType','string');
    [ok,ix]=ismember(pm.Metabolites,t.Metabolites); assert(all(ok));
    plateCalls(:,k)=t.Growth(ix);
end
assert(all(ismember(plateCalls(:),[0 1])));
for j=1:numel(names), assert(isequal(max(plateCalls(:,strains==names(j)),[],2),unionCalls(:,j))); end
rep=["ST1-12","ST1-68","ST1-75","VPI"]; indices=cell(1,4);
for j=1:4, indices{j}=find(strains==rep(j)); assert(numel(indices{j})==3); end
[a,b,c,d]=ndgrid(1:3,1:3,1:3,1:3); choices=[a(:),b(:),c(:),d(:)];
rows=cell(486,11); allPredictors=cell(81*21,6);
for k=1:81
    selected=zeros(1,numel(names));
    for j=1:numel(names), selected(j)=find(strains==names(j),1); end
    for j=1:4, selected(names==rep(j))=indices{j}(choices(k,j)); end
    [fraction,~,den]=nutrient_overlap(plateCalls(:,selected),col,vpi,aa);
    for i=1:21
        allPredictors((k-1)*21+i,:)={k,strainNames(i),fraction(i,1),fraction(i,2),fraction(i,3),den(2)};
    end
    for j=1:3
        for adjusted=0:1
            controls=zeros(21,0); if adjusted, controls=disease; end
            r=partial_spearman_pgfam(fraction(:,j),y,controls);
            row=(k-1)*6+(j-1)*2+adjusted+1;
            rows(row,:)={k,j,logical(adjusted),den(j),r.Rho,r.PValue,r.DF,r.Valid,r.Reason, ...
                source(selected(vpi)),strjoin(source(selected(ismember(names,rep(1:3)))), ';')};
        end
    end
end
t=cell2table(rows,'VariableNames',{'Selection','ResourceClass','DiseaseAdjusted','VPIPositiveCalls', ...
    'Rho','PValue','DF','Valid','Reason','VPIPlate','ReplicatedST1Plates'});
writetable(t,fullfile(out,'equal_plate_all_81_selections.csv'));
writetable(cell2table(allPredictors,'VariableNames', ...
    {'Selection','Strain','OverlapAll','OverlapAmino','OverlapNonAmino','VPIAminoCalls'}), ...
    fullfile(out,'equal_plate_all_predictors.csv'));
summary=cell(6,7);
for j=1:3
    for adjusted=0:1
        keep=t.ResourceClass==j & t.DiseaseAdjusted==logical(adjusted); values=t.Rho(keep & t.Valid);
        summary((j-1)*2+adjusted+1,:)={j,logical(adjusted),numel(values),81-numel(values), ...
            min(values),median(values),max(values)};
    end
end
writetable(cell2table(summary,'VariableNames',{'ResourceClass','DiseaseAdjusted','ValidSelections', ...
    'InvalidSelections','MinimumRho','MedianRho','MaximumRho'}),fullfile(out,'equal_plate_summary.csv'));
end
