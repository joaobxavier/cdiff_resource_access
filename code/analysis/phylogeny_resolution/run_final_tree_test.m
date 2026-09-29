function run_final_tree_test(name)
% Full-alignment ML with ModelFinder and dual branch-support measures.
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
root=fullfile(package_root(),'results','phylogeny');
folder=fullfile(root,'tree_tests');
exe=fullfile(getenv('CDIFF_PHYLO_BIN'),'iqtree2');
input=fullfile(folder,[name '.fasta']);prefix=fullfile(folder,name);
% IQ-TREE writes a provisional tree before bootstrap refinement finishes.
% A treefile alone is not evidence that an interrupted run completed.
if isfile([prefix '.treefile']) && isfile([prefix '.log']) && ...
        contains(fileread([prefix '.log']),'Analysis results written to:')
    fprintf('%s already complete.\n',name);return;
end
outgroup='';if ~contains(name,'st1_same'),outgroup='-o VPI10463';end
cmd=sprintf('"%s" -s "%s" -m MFP -B 1000 -bnni -alrt 1000 -T 4 -seed 20260919 %s --prefix "%s" > "%s.console.log" 2>&1',exe,input,outgroup,prefix,prefix);
fprintf('Running %s\n',name);[status,msg]=system(cmd);assert(status==0,'%s',msg);
fprintf('Completed %s\n',name);
end
