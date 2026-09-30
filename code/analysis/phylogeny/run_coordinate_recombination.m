function run_coordinate_recombination(mode)
% Run Gubbins on whole-genome coordinate alignment, ST1 strains only.
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(package_root(),'results','phylogeny');
work=package_phylogeny_work();bin=getenv('CDIFF_PHYLO_BIN');
old=pwd;c=onCleanup(@()cd(old));cd(work);
if nargin==0,mode='initial';end
if strcmp(mode,'final')
    prefix='st1_qc_recombination';input='coordinate_qc_st1.fasta';iterations=20;
else
    prefix='st1_recombination';input='mummer_st1.fasta';iterations=5;
end
logFile=[prefix '.console.log'];
cmd=sprintf('PATH="%s:$PATH" "%s/run_gubbins.py" --prefix %s --tree-builder iqtree --model GTRGAMMA --threads 4 --seed 20260919 --iterations %d --no-cleanup %s > %s 2>&1',bin,bin,prefix,iterations,input,logFile);
if ~isfile([prefix '.final_tree.tre'])
    [status,msg]=system(cmd);assert(status==0,'%s; see gubbins.log',msg);
end
files=dir([prefix '.*']);
for i=1:numel(files),if ~files(i).isdir,copyfile(files(i).name,fullfile(out,files(i).name));end,end
fprintf('Gubbins completed; results copied to %s\n',out);
end
