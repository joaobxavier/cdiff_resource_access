function prepare_genome_coordinate_check
% Stage accession-verified references and a repeat-masked whole-genome alignment.
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
root=fullfile(package_root(),'results','phylogeny');
in=fullfile(package_root(),'input','genomes'); out=root;
if ~isfolder(out),mkdir(out);end
refs=["r20291","vpi"]; accessions=["GCF_000027105.1","GCF_015238635.1"];
for i=1:2
    fasta=dir(fullfile(in,accessions(i)+'*.fna'));
    assert(isscalar(fasta)); f=fastaread(fullfile(fasta.folder,fasta.name));
    assert(isscalar(f));
    if i==1,assert(contains(f.Header,'R20291'));else,assert(contains(f.Header,'ATCC 43255'));end
    copyfile(fullfile(fasta.folder,fasta.name),fullfile(out,refs(i)+"_reference.fasta"));
end
% SKA's file-list parser and Gubbins helpers require whitespace-free paths.
% Stage immutable copies at a path without the project directory's spaces.
work=package_phylogeny_work();
if ~isfolder(work),mkdir(work);end
f=dir(fullfile(in,'ST1_*.fasta'));assert(numel(f)==21);
listPath=fullfile(work,'assemblies.txt');fid=fopen(listPath,'w');
for i=1:numel(f)
    [~,name]=fileparts(f(i).name);
    copyfile(fullfile(f(i).folder,f(i).name),fullfile(work,f(i).name));
    fprintf(fid,'%s\t%s\n',name,fullfile(work,f(i).name));
end
fclose(fid);
copyfile(fullfile(out,'r20291_reference.fasta'),fullfile(work,'reference.fasta'));
copyfile(fullfile(out,'vpi_reference.fasta'),fullfile(work,'VPI10463.fasta'));
copyfile(fullfile(in,'GCA_001995155.1_VPI10463.fasta'),fullfile(work,'VPI10463_old.fasta'));
bin=getenv('CDIFF_PHYLO_BIN');
cmd=sprintf('PATH="%s:$PATH" "%s/generate_ska_alignment.py" --reference %s/reference.fasta --input %s/assemblies.txt --out %s/st1_k17.aln --k 17 --threads 4 --no-cleanup > %s/ska.log 2>&1',bin,bin,work,work,work,work);
if ~isfile(fullfile(work,'st1_k17.aln'))
    [status,msg]=system(cmd);assert(status==0,'%s; see ska.log',msg);
end
copyfile(listPath,fullfile(out,'assemblies_staging_manifest.txt'));
copyfile(fullfile(work,'st1_k17.aln'),fullfile(out,'st1_k17.aln'));
copyfile(fullfile(work,'ska.log'),fullfile(out,'ska.log'));
a=fastaread(fullfile(out,'st1_k17.aln'));seq=upper(char({a.Sequence}));
assert(size(seq,1)==21 && size(seq,2)==4191339);
names=string({a.Header})'; valid=ismember(seq,'ACGT');
t=table(names,sum(valid,2),mean(valid,2),'VariableNames',{'Strain','CallableBases','CallableFraction'});
writetable(t,fullfile(out,'ska_callability.csv'));disp(t);
fprintf('Complete ST1 columns: %d; variable complete columns: %d\n',sum(all(valid,1)),sum(all(valid,1)&any(seq~=seq(1,:),1)));
end
