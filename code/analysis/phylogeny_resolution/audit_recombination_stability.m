function audit_recombination_stability
% Distinguish branch-length nonconvergence from instability of excluded sites.
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(package_root(),'results','phylogeny');
work=package_phylogeny_work();
Iteration=(1:20)';UnionMaskedBases=zeros(20,1);ChangedSincePrevious=nan(20,1);
previous=[];masks=false(20,4191339);
for i=1:20
    path=fullfile(work,sprintf('coordinate_qc_st1.iteration_%d.tre.gff',i));
    if i==20,path=fullfile(out,'st1_qc_recombination.recombination_predictions.gff');end
    lines=splitlines(string(fileread(path)));mask=false(1,4191339);
    for j=1:numel(lines)
        if startsWith(lines(j),'#')||strlength(lines(j))==0,continue;end
        f=split(lines(j),sprintf('\t'));mask(str2double(f(4)):str2double(f(5)))=true;
    end
    UnionMaskedBases(i)=sum(mask);masks(i,:)=mask;
    if ~isempty(previous),ChangedSincePrevious(i)=sum(mask~=previous);end
    previous=mask;
    copyfile(path,fullfile(out,sprintf('gubbins_iteration_%02d.gff',i)));
end
t=table(Iteration,UnionMaskedBases,ChangedSincePrevious);disp(t);
writetable(t,fullfile(out,'gubbins_mask_stability.csv'));
fprintf('Final five masks identical: %d\n',all(all(masks(16:20,:)==masks(20,:))));
save(fullfile(out,'gubbins_iteration_masks.mat'),'masks');
end
