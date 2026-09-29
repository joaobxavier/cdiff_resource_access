function rebuild_phylogeny()
bin=getenv('CDIFF_PHYLO_BIN');
assert(isfolder(bin),'Set CDIFF_PHYLO_BIN to the genomics environment bin directory.');
for exe={'nucmer','delta-filter','iqtree2','run_gubbins.py','generate_ska_alignment.py'}
    assert(isfile(fullfile(bin,exe{1})),'Missing executable %s.',exe{1});
end
prepare_genome_coordinate_check();
build_mummer_coordinate_alignment();
finalize_coordinate_alignment();
run_coordinate_recombination('final');
audit_recombination_stability();
prepare_final_tree_tests();
for name=["masked_new22","unmasked_new22","masked_st1_same_columns","masked_new22_common","masked_old22_common"]
    run_final_tree_test(char(name));
end
% Normalize duplicate support labels for MATLAB without fitting to traits.
folder=fullfile(package_root(),'results','phylogeny','tree_tests');
files=dir(fullfile(folder,'*.treefile'));
for i=1:numel(files)
    text=fileread(fullfile(folder,files(i).name));
    [starts,ends,~,~,tokens]=regexp(text,'\)([0-9.]+/[0-9.]+):');
    for k=numel(starts):-1:1
        text=[text(1:starts(k)-1) sprintf(')%s_node%d:',tokens{k}{1},k) text(ends(k)+1:end)];
    end
    text=strrep(text,'ST1_','ST1.');text=strrep(text,'VPI10463','VPI');
    target=fullfile(folder,strrep(files(i).name,'.treefile','.analysis.nwk'));
    fid=fopen(target,'w');assert(fid>=0);fprintf(fid,'%s',text);fclose(fid);
end
end
