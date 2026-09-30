function prepare_final_tree_tests
% Apply the union of inferred ST1 recombination intervals to every taxon.
% Complete-column alignments retain invariant sites for ML branch lengths.
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
root=fullfile(package_root(),'results','phylogeny');
load(fullfile(root,'coordinate_qc.mat'),'names','seq','qcMask');
gff=splitlines(string(fileread(fullfile(root,'st1_qc_recombination.recombination_predictions.gff'))));
mask=false(1,size(seq,2));events=0;
for i=1:numel(gff)
    if startsWith(gff(i),'#')||strlength(gff(i))==0,continue;end
    f=split(gff(i),sprintf('\t'));mask(str2double(f(4)):str2double(f(5)))=true;events=events+1;
end
valid=ismember(seq,'ACGT');st1=1:21;new22=1:22;old22=[1:21 23];
allComplete=all(valid,1);newComplete=all(valid(new22,:),1);
clean=newComplete & ~mask;common=allComplete & ~mask;
out=fullfile(root,'tree_tests');if ~isfolder(out),mkdir(out);end
writeFasta(fullfile(out,'masked_new22.fasta'),names(new22),seq(new22,clean));
writeFasta(fullfile(out,'unmasked_new22.fasta'),names(new22),seq(new22,newComplete));
writeFasta(fullfile(out,'masked_st1_same_columns.fasta'),names(st1),seq(st1,clean));
writeFasta(fullfile(out,'masked_new22_common.fasta'),names(new22),seq(new22,common));
oldNames=names(old22);oldNames(end)="VPI10463";
writeFasta(fullfile(out,'masked_old22_common.fasta'),oldNames,seq(old22,common));
pair=valid(names=="ST1_75",:)&valid(names=="ST1_68",:);
diff=seq(names=="ST1_75",:)~=seq(names=="ST1_68",:);
s=struct('ReferencePositions',size(seq,2),'RecombinationEvents',events, ...
    'UnionRecombinationBases',sum(mask),'QCExcludedBases',sum(qcMask), ...
    'UnmaskedNew22Sites',sum(newComplete),'MaskedNew22Sites',sum(clean), ...
    'MaskedCommonOutgroupSites',sum(common), ...
    'UnmaskedST1VariableSites',sum(any(seq(st1,:)~=seq(1,:),1)&newComplete), ...
    'MaskedST1VariableSites',sum(any(seq(st1,:)~=seq(1,:),1)&clean), ...
    'FocalPairCallable',sum(pair),'FocalPairSNPs',sum(pair&diff), ...
    'FocalPairNonrecombinantCallable',sum(pair&~mask), ...
    'FocalPairNonrecombinantSNPs',sum(pair&~mask&diff), ...
    'FocalPairTreeSitesSNPs',sum(clean&diff));
disp(s);fid=fopen(fullfile(out,'alignment_summary.json'),'w');fprintf(fid,'%s\n',jsonencode(s,PrettyPrint=true));fclose(fid);
pairRows=nchoosek(1:21,2);SNPs=zeros(size(pairRows,1),1);Callable=SNPs;
for i=1:numel(SNPs)
    a=pairRows(i,1);b=pairRows(i,2);keep=valid(a,:)&valid(b,:)&~mask;
    Callable(i)=sum(keep);SNPs(i)=sum(seq(a,keep)~=seq(b,keep));
end
writetable(table(names(pairRows(:,1)),names(pairRows(:,2)),Callable,SNPs,'VariableNames',{'Strain1','Strain2','CallableSites','NonrecombinantSNPs'}),fullfile(out,'pairwise_nonrecombinant_snps.csv'));
save(fullfile(out,'column_masks.mat'),'mask','newComplete','clean','common','s');
end
function writeFasta(path,names,seq)
fid=fopen(path,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));
for i=1:numel(names),fprintf(fid,'>%s\n%s\n',names(i),seq(i,:));end
end
