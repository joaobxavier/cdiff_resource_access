function assess_phylogeny_adjusted_amino_access(mode,nBoot)
% Assess whether amino-acid access adds information beyond disease + ancestry.
% V = S + v*((1-w)*I + w*C), C = root-shared branch covariance,
% normalized to mean diagonal one. S is the full protection-score covariance.
% w is a variance-mixture fraction, NOT Pagel's lambda. Fit by ML; explicitly
% evaluate boundaries. Null parametric bootstrap refits both covariance models.
% Write model, covariance, influence and bootstrap results; no extra figure.
if nargin<1,mode="prepare";end
if nargin<2,nBoot=1999;end
localRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
project=fileparts(localRoot);addpath(fullfile(localRoot,'analysis'));
out=fullfile(project,'results','analyses','phylogeny_adjusted_amino_access');
if ~isfolder(out),mkdir(out);end
cache=fullfile(out,'analysis_state.mat');
if mode=="prepare"
    t=readtable(fullfile(project,'results', 'tables', 'main', ...
        'figure4_phylogeny_aligned_traits.csv'),'TextType','string');
    t=sortrows(t,'Strain');assert(height(t)==21);
    y=t.Protection_Effect;d=t.Mono_Disease_Effect;
    x=t.Amino_Acid_Relative_Breadth_Advantage;
    X0=[ones(21,1),d];X1=[X0,x];
    S=recoverCovariance(t.Strain,y);
    treeRoot=fullfile(package_root(),'results','phylogeny','tree_tests');
    keys=["masked_new22","unmasked_new22","masked_new22_common","masked_old22_common"];
    Cs=cell(5,1);rows=cell(0,13);
    for k=1:4
        Cs{k}=treeCovariance(fullfile(treeRoot,keys(k)+'.analysis.nwk'),t.Strain,false);
        [f0,f1,row]=compare(y,X0,X1,S,Cs{k},"fit",keys(k));
        rows(end+1,:)=row; %#ok<AGROW>
        if k==1,primary0=f0;primary1=f1;end
    end
    Cs{5}=treeCovariance(fullfile(treeRoot,'masked_new22.analysis.nwk'),t.Strain,true);
    [~,~,row]=compare(y,X0,X1,S,Cs{5},"fit","unsupported_branches_collapsed");rows(end+1,:)=row;
    C=Cs{1};
    [~,~,row]=compare(y,X0,X1,S,C,0,"no_phylogenetic_component");rows(end+1,:)=row;
    [~,~,row]=compare(y,X0,X1,S,C,1,"forced_Brownian_component");rows(end+1,:)=row;
    [~,~,row]=compare(y,X0,X1,zeros(21),C,"fit","without_score_sampling_covariance");rows(end+1,:)=row;
    [~,~,row]=compare(y,X0,X1,zeros(21),C,1,"strict_Brownian_without_sampling_error");rows(end+1,:)=row;
    results=cell2table(rows,'VariableNames',{'Analysis','N','AA_Effect','CI_Lower','CI_Upper', ...
        'Approx_LRT_P','LRT','Delta_AICc','Null_Phylo_Fraction','Full_Phylo_Fraction', ...
        'Null_Residual_Variance','Full_Residual_Variance','Approx_Wald_P'});
    writetable(results,fullfile(out,'model_comparisons.csv'));disp(results);
    influence=cell(21,13);prediction=zeros(21,4);
    for k=1:21
        keep=true(21,1);keep(k)=false;
        [r,f,influence(k,:)]=compare(y(keep),X0(keep,:),X1(keep,:),S(keep,keep), ...
            C(keep,keep),"fit","without_"+t.Strain(k));
        prediction(k,1)=X0(k,:)*r.Beta;
        prediction(k,2)=X1(k,:)*f.Beta;
        % Conditional predictions of the held-out estimated score, including
        % phylogenetic covariance and shared-control sampling covariance.
        vr=S(k,keep)+r.Variance*r.Weight*C(k,keep);
        vf=S(k,keep)+f.Variance*f.Weight*C(k,keep);
        prediction(k,3)=prediction(k,1)+vr*(r.V\(y(keep)-X0(keep,:)*r.Beta));
        prediction(k,4)=prediction(k,2)+vf*(f.V\(y(keep)-X1(keep,:)*f.Beta));
    end
    influence=cell2table(influence,'VariableNames',results.Properties.VariableNames);
    writetable(influence,fullfile(out,'leave_one_strain_out.csv'));
    predictions=table(t.Strain,y,prediction(:,1),prediction(:,2),prediction(:,3),prediction(:,4), ...
        'VariableNames',{'Strain','Observed_Protection','Disease_Mean','Disease_AA_Mean', ...
        'Disease_Conditional','Disease_AA_Conditional'});
    writetable(predictions,fullfile(out,'held_out_predictions.csv'));
    fprintf('LOSO RMSE [mean0 mean1 conditional0 conditional1]: ');disp(sqrt(mean((prediction-y).^2)));
    % Numerical invariance to reversing the strain order.
    rev=21:-1:1;check=fitModel(y(rev),X1(rev,:),S(rev,rev),C(rev,rev),"fit");
    assert(max(abs(check.Beta-primary1.Beta))<1e-6);
    assert(abs(check.NLL-primary1.NLL)<1e-7);
    assert(min(eig(S))>0 && min(eig(C))>-1e-10);
    writetable(t,fullfile(out,'strain_inputs.csv'));
    writematrix(S,fullfile(out,'protection_sampling_covariance.csv'));
    writematrix(C,fullfile(out,'phylogenetic_shared_branch_covariance.csv'));
    save(cache,'t','y','d','x','X0','X1','S','C','Cs','primary0','primary1','results','influence','prediction');
elseif mode=="bootstrap"
    a=load(cache);rng(20260919,'twister');
    z=randn(numel(a.y),nBoot);nullY=a.X0*a.primary0.Beta+chol(a.primary0.V,'lower')*z;
    nullLRT=zeros(nBoot,1);nullEffect=zeros(nBoot,1);
    observed=2*(a.primary0.NLL-a.primary1.NLL);
    for b=1:nBoot
        r=fitModel(nullY(:,b),a.X0,a.S,a.C,"fit");
        f=fitModel(nullY(:,b),a.X1,a.S,a.C,"fit");
        assert(f.NLL<=r.NLL+1e-6,'Optimizer failed nested likelihood check.');
        nullLRT(b)=max(0,2*(r.NLL-f.NLL));nullEffect(b)=f.Beta(3);
        if mod(b,100)==0,fprintf('Bootstrap %d/%d completed\n',b,nBoot);end
    end
    p=(1+sum(nullLRT>=observed-1e-10))/(nBoot+1);
    mcse=sqrt(p*(1-p)/(nBoot+1));
    save(fullfile(out,'bootstrap.mat'),'nullLRT','nullEffect','nBoot','p','mcse','observed');
    writetable(table(observed,nBoot,p,mcse),fullfile(out,'bootstrap_summary.csv'));
    fprintf('Null-refitted parametric bootstrap p=%.6f (MC SE %.6f), B=%d\n',p,mcse,nBoot);
else
    error('Unknown mode.');
end
end

function S=recoverCovariance(strains,expected)
raw=readtable(project_data_file('processed','mouse','scores','ProtectionScreen_CDI_mouse.csv'),'TextType','string');
raw=raw(raw.experiment~="ks65" & (endsWith(raw.cdiffstrain,".vpi") | raw.cdiffstrain=="vpi"),:);
missing=ismissing(raw.relweight);assert(all(raw.death(missing)==1));
raw.relweight(missing)=0;
raw.exp_id=categorical(string(raw.experiment)+"_"+string(raw.cdiffstrain)+"_"+string(raw.mouse));
raw.cdiffstrain=categorical(raw.cdiffstrain);
raw.cdiffstrain=reordercats(raw.cdiffstrain,["vpi";setdiff(string(categories(raw.cdiffstrain)),"vpi")]);
model=fitlme(raw,'relweight ~ cdiffstrain + (1|day) + (1|exp_id)');
coef=model.Coefficients;names=string(coef.Name);keep=startsWith(names,"cdiffstrain_st1.");
names=upper(replace(erase(erase(names(keep),"cdiffstrain_"),".vpi"),'.','-'));
[ok,idx]=ismember(strains,names);assert(all(ok));indices=find(keep);indices=indices(idx);
assert(max(abs(coef.Estimate(indices)-expected))<1e-8);
S=model.CoefficientCovariance(indices,indices);S=(S+S')/2;
end

function C=treeCovariance(file,strains,collapse)
tree=phytreeread(file);names=string(get(tree,'LeafNames'));
tree=reroot(tree,find(names=="VPI"));names=string(get(tree,'LeafNames'));
tree=prune(tree,find(names=="VPI"));
n=get(tree,'NumLeaves');p=get(tree,'Pointers');len=get(tree,'Distances');
labels=string(get(tree,'NodeNames'));desc=false(2*n-1,n);desc(1:n,:)=eye(n)>0;
for k=1:n-1,desc(n+k,:)=desc(p(k,1),:) | desc(p(k,2),:);end
original=len;
if collapse
    for k=n+1:2*n-2
        support=sscanf(labels(k),'%f/%f');
        if numel(support)~=2 || support(1)<80 || support(2)<95,len(k)=0;end
    end
end
len(end)=0;C=double(desc)'*diag(len)*double(desc);
if ~collapse
    dist=squareform(pdist(tree,'Nodes','leaves'));
    derived=diag(C)+diag(C)'-2*C;
    assert(max(abs(dist-derived),[],'all')<1e-12,'Branch covariance does not reproduce patristic distances.');
end
assert(all(original>=0));
names=replace(string(get(tree,'LeafNames')),'.','-');[ok,ix]=ismember(strains,names);assert(all(ok));
C=C(ix,ix);C=C/mean(diag(C));C=(C+C')/2;
end

function [r,f,row]=compare(y,X0,X1,S,C,weight,label)
r=fitModel(y,X0,S,C,weight);f=fitModel(y,X1,S,C,weight);
assert(f.NLL<=r.NLL+1e-6);
lrt=max(0,2*(r.NLL-f.NLL));p=chi2cdf(lrt,1,'upper');
critical=tinv(.975,numel(y)-size(X1,2));
se=sqrt(f.BetaCov(end,end));ci=f.Beta(end)+[-1,1]*critical*se;
wald=2*tcdf(-abs(f.Beta(end)/se),numel(y)-size(X1,2));
row={label,numel(y),f.Beta(end),ci(1),ci(2),p,lrt,f.AICc-r.AICc, ...
    r.Weight,f.Weight,r.Variance,f.Variance,wald};
end

function fit=fitModel(y,X,S,C,weight)
% Profile the residual scale, then the mixture weight; boundaries explicit.
if isnumeric(weight)
    [value,v]=profileWeight(weight,y,X,S,C);w=weight;nv=1;
else
    grid=linspace(0,1,11);values=zeros(size(grid));scales=values;
    for j=1:numel(grid),[values(j),scales(j)]=profileWeight(grid(j),y,X,S,C);end
    [value,i]=min(values);w=grid(i);v=scales(i);
    for j=2:numel(grid)-1
        if values(j)<=values(j-1) && values(j)<=values(j+1)
            [wc,vc]=fminbnd(@(ww)profileWeight(ww,y,X,S,C),grid(j-1),grid(j+1), ...
                optimset('Display','off','TolX',1e-6));
            if vc<value,value=vc;w=wc;[~,v]=profileWeight(w,y,X,S,C);end
        end
    end
    nv=2;
end
[~,beta,bc,V]=objective(v,w,y,X,S,C);
k=size(X,2)+nv;n=numel(y);
fit=struct('NLL',value,'Beta',beta,'BetaCov',bc,'V',V, ...
    'Weight',w,'Variance',v,'AICc',2*value+2*k+2*k*(k+1)/(n-k-1));
end

function [nll,v]=profileWeight(w,y,X,S,C)
if ~any(S,'all')
    K=(1-w)*eye(numel(y))+w*C;L=chol(K,'lower');
    xx=L\X;yy=L\y;b=xx\yy;v=sum((yy-xx*b).^2)/numel(y);
    nll=objective(v,w,y,X,S,C);return
end
upper=max(100*var(y),1);
[v,nll]=fminbnd(@(s)objective(s,w,y,X,S,C),0,upper,optimset('Display','off','TolX',1e-7));
zero=objective(0,w,y,X,S,C);
if zero<nll,v=0;nll=zero;end
assert(v<upper*.99,'Residual variance reached search bound.');
end

function [nll,b,bc,V]=objective(v,w,y,X,S,C)
V=S+v*((1-w)*eye(numel(y))+w*C);L=chol(V,'lower');
xx=L\X;yy=L\y;b=xx\yy;r=yy-xx*b;
nll=.5*(numel(y)*log(2*pi)+2*sum(log(diag(L)))+r'*r);
if nargout>2,bc=(xx'*xx)\eye(size(X,2));end
end
