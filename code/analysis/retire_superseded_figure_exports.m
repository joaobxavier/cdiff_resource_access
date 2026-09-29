function retire_superseded_figure_exports()
% Preserve prior output files when rerunning a Draft 43 results directory.
root=package_root();
stems={'figureS2_phylogenetic_context', ...
    'figureS3_pm1_detail_and_generalization_limits','figureS4_contextual_analyses'};
folders={'figures','vectors','figures/editable_vectors'};
archive=fullfile(root,'results','historical_figure_exports', ...
    char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
for i=1:numel(folders)
    for j=1:numel(stems)
        for ext={'.png','.pdf','.svg','.fig'}
            source=fullfile(root,'results',folders{i},[stems{j},ext{1}]);
            if isfile(source)
                destination=fullfile(archive,folders{i});
                if ~isfolder(destination);mkdir(destination);end
                movefile(source,fullfile(destination,[stems{j},ext{1}]));
            end
        end
    end
end
end
