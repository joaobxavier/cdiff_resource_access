function run_all(mode)
% Rebuild the manuscript analysis without access to the author's workspace.
if nargin<1,mode='standard';end
root=fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root,'code')));
assert(ismember(string(mode),["standard","full","tests"]),'Unknown run mode.');
for folder={'results','results/tables/main','results/figures','results/vectors','results/logs'}
    if ~isfolder(fullfile(root,folder{1})),mkdir(fullfile(root,folder{1}));end
end
setappdata(0,'Draft39VectorOutput',fullfile(root,'results','vectors'));
set(0,'DefaultFigureVisible','off');
if strcmp(mode,'tests'),test_release();return;end
retire_superseded_figure_exports();
validate_inputs();
check_release_code();
if strcmp(mode,'full')
    rebuild_phylogeny();
else
    copyfile(fullfile(root,'cache','phylogeny'),fullfile(root,'results','phylogeny'));
end
stages={
    'mouse','analyze_mouse_weight_survival_screen';
    'qpcr','rebuild_qpcr';
    'fecal_cfu','analyze_fecal_cfu_screen';
    'figure1','generate_theory_led_thirty_ninth_draft_figure1';
    'figure2','generate_theory_led_thirty_ninth_draft_figure2';
    'biolog','rebuild_pm1';
    'pm1_equal_replication','assess_pm1_equal_replication';
    'intracellular','analyze_intracellular_matched_batches';
    'intracellular_plsda','plot_intracellular_plsda_significant_background';
    'extracellular','analyze_extracellular_biolog_concordance';
    'figure3','generate_draft39_gcms_figure';
    'pm1_traits','prepare_package_pm1';
    'pm1_models','compare_simple_pm1_models_for_paper';
    'pm1_specificity','analyze_amino_acid_breadth_specificity';
    'figure4_S2','generate_draft42_figure4';
    'tree_sensitivity','summarize_final_tree_tests';
    'ancestry_models','prepare_ancestry';
    'ancestry_bootstrap','bootstrap_ancestry';
    'adaptive_flow','flow_adaptive';
    'innate_flow','flow_innate';
    'figure5_S3','generate_host_figures'};
status=cell(0,3);
for i=1:size(stages,1)
    fprintf('\n=== %s (%d/%d) ===\n',stages{i,1},i,size(stages,1));
    diary(fullfile(root,'results','logs',[stages{i,1} '.log']));
    start=tic;
    try
        execute_stage(stages{i,2});
        status(end+1,:)={stages{i,1},'PASS',toc(start)}; %#ok<AGROW>
    catch err
        status(end+1,:)={stages{i,1},'FAIL',toc(start)}; %#ok<AGROW>
        writetable(cell2table(status,'VariableNames',{'Stage','Status','Seconds'}),fullfile(root,'results','run_status.csv'));
        report=getReport(err,'extended','hyperlinks','off');
        fid=fopen(fullfile(root,'results','last_error.txt'),'w');fprintf(fid,'%s',report);fclose(fid);
        fprintf('%s\n',report);diary off;
        rethrow(err);
    end
    diary off;
    writetable(cell2table(status,'VariableNames',{'Stage','Status','Seconds'}),fullfile(root,'results','run_status.csv'));
end
test_release();
end

function execute_stage(name)
% Local script clears cannot erase the master's progress state.
eval([name ';']);
end
