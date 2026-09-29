function finalize_coordinate_alignment
% Mask duplicated reference segments and independently discordant base calls.
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
root=fullfile(package_root(),'results','phylogeny');
work=package_phylogeny_work();bin=getenv('CDIFF_PHYLO_BIN');
load(fullfile(root,'coordinate_alignment.mat'),'names','seq','refseq');
prefix=string(fullfile(work,'reference_self'));
if ~isfile(prefix+'.filtered.delta')
    cmd=sprintf('"%s/nucmer" --maxmatch --nosimplify -t 4 -p %s %s/reference.fasta %s/reference.fasta > %s.log 2>&1',bin,prefix,work,work,prefix);
    [status,msg]=system(cmd);assert(status==0,'%s',msg);
    cmd=sprintf('"%s/delta-filter" -i 95 -l 500 %s.delta > %s.filtered.delta',bin,prefix,prefix);
    [status,msg]=system(cmd);assert(status==0,'%s',msg);
end
lines=splitlines(string(fileread(prefix+'.filtered.delta')));repeatMask=false(1,numel(refseq));
n=3;
while n<=numel(lines)
    line=strtrim(lines(n));n=n+1;
    if strlength(line)==0||startsWith(line,'>'),continue;end
    v=sscanf(line,'%d');assert(numel(v)==7);
    if v(1)~=v(3)||v(2)~=v(4)
        repeatMask(v(1):v(2))=true;
        repeatMask(min(v(3:4)):max(v(3:4)))=true;
    end
    while str2double(lines(n))~=0,n=n+1;end
    n=n+1;
end
copyfile(prefix+'.filtered.delta',fullfile(root,'reference_self.filtered.delta'));
ska=fastaread(fullfile(root,'st1_k17.aln'));sn=string({ska.Header});
discordantMask=false(size(repeatMask));
for i=1:21
    a=upper(ska(sn==names(i)).Sequence);b=seq(i,:);
    discordantMask=discordantMask|(ismember(a,'ACGT')&ismember(b,'ACGT')&a~=b);
end
qcMask=repeatMask|discordantMask;
seq(:,qcMask)='N';
writeFasta(fullfile(root,'coordinate_qc_all23.fasta'),names,seq);
writeFasta(fullfile(work,'coordinate_qc_st1.fasta'),names(1:21),seq(1:21,:));
copyfile(fullfile(work,'coordinate_qc_st1.fasta'),fullfile(root,'coordinate_qc_st1.fasta'));
fprintf('Repeat mask: %d bp; discordant sites: %d; combined: %d\n',sum(repeatMask),sum(discordantMask),sum(qcMask));
valid=ismember(seq,'ACGT');
writetable(table(names,sum(valid,2),mean(valid,2),'VariableNames',{'Strain','CallableBases','CallableFraction'}),fullfile(root,'coordinate_qc_callability.csv'));
save(fullfile(root,'coordinate_qc.mat'),'names','seq','refseq','repeatMask','discordantMask','qcMask');
end
function writeFasta(path,names,seq)
fid=fopen(path,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));
for i=1:numel(names),fprintf(fid,'>%s\n%s\n',names(i),seq(i,:));end
end
