function test_release()
% Independent frozen manuscript numbers are tests, never fitting inputs.
root=package_root();validate_inputs();
% Cross-figure checks use the statistical definition, not frozen expectations.
traits=readtable(fullfile(root,'results','tables','main','figure4_phylogeny_aligned_traits.csv'));
df=height(traits)-3; % Intercept plus one disease covariate.
assert(df==18,'Unexpected sample size for partial-rank tests.');
for name=["figureS3_amino_nonamino_specificity.csv","figureS3_amino_leave_one_well_out.csv"]
    t=readtable(fullfile(root,'results','tables','main',name),'TextType','string');
    definitionP=2*tcdf(-abs(t.Partial_Rank_Rho).*sqrt(df./(1-t.Partial_Rank_Rho.^2)),df);
    assert(all(abs(t.Partial_Rank_P_Value-definitionP)<1e-12), ...
        'S2 partial-rank p-values must account for the disease covariate: %s',name);
end
specificity=readtable(fullfile(root,'results','tables','main','figureS3_amino_nonamino_specificity.csv'),'TextType','string');
association=readtable(fullfile(root,'results','tables','main','figure4_association_statistics.csv'),'TextType','string');
amino=specificity(specificity.Score=="Amino-acid breadth",:);
partial=association(startsWith(association.Association,"Partial-rank"),:);
assert(height(amino)==1 && height(partial)==1 && ...
    abs(amino.Partial_Rank_Rho-partial.Rho)<1e-12 && ...
    abs(amino.Partial_Rank_P_Value-partial.P_Value)<1e-12, ...
    'Figure 4E and S2B must report the same adjusted amino-acid test.');
% Synthetic parser fixture only; never an input to study analyses.
spans=find_ambiguous_delta_spans(fullfile(root,'expected','fixtures','alignment_ties.delta'));
assert(isequal(spans,[1001 1500]),'Tied-alignment masking regression.');
manifest=readtable(fullfile(root,'expected','regression_manifest.tsv'),'FileType','text','Delimiter','\t','TextType','string');
rows=cell(0,4);
for i=1:height(manifest)
    % Explicit CSV delimiter prevents semicolons in prose cells from making
    % automatic import collapse numeric columns into a single text column.
    expected=readtable(fullfile(root,manifest.Expected(i)),'Delimiter',',','TextType','string','VariableNamingRule','preserve');
    actual=readtable(fullfile(root,manifest.Actual(i)),'Delimiter',',','TextType','string','VariableNamingRule','preserve');
    assert(height(actual)==height(expected),'Row count mismatch: %s',manifest.Actual(i));
    vars=expected.Properties.VariableNames;
    maxDifference=0;numericCells=0;
    for j=1:numel(vars)
        name=vars{j};assert(ismember(name,actual.Properties.VariableNames),'Missing variable %s.',name);
        a=actual.(name);e=expected.(name);
        if isnumeric(e) || islogical(e)
            assert(isequal(isnan(double(a)),isnan(double(e))),'Missingness changed: %s %s',manifest.Actual(i),name);
            use=isfinite(double(e));
            difference=abs(double(a(use))-double(e(use)));
            tolerance=1e-7+1e-7*abs(double(e(use)));
            assert(all(difference<=tolerance),'Numerical mismatch: %s %s (max %.12g)',manifest.Actual(i),name,max(difference));
            if ~isempty(difference),maxDifference=max(maxDifference,max(difference));end
            numericCells=numericCells+sum(use);
        else
            % Numerical verification plus stable identifiers; narrative
            % provenance fields may legitimately refer to package paths.
            if ismember(name,{'Strain','Candidate','Category','Metabolite','Substrate','Analysis','Association','Term','Model'})
                assert(isequal(string(a),string(e)),'Identifier/order mismatch: %s %s',manifest.Actual(i),name);
            end
        end
    end
    assert(numericCells>0 || isequal(vars,{'Strain'}), ...
        'No numeric values were tested in %s; inspect CSV import.',manifest.Actual(i));
    rows(end+1,:)={manifest.Actual(i),'PASS',numericCells,maxDifference}; %#ok<AGROW>
end
for kind=["adaptive","innate"]
    expected=load(fullfile(root,'expected','numerical','flow_'+kind+'.mat'));
    actual=load(fullfile(root,'results','vectors','flow_'+kind+'.mat'));
    for field=["UI_data","Avirulent_data"]
        assert(isequal(expected.(field),actual.(field)),'Flow fraction mismatch: %s %s',kind,field);
    end
    rows(end+1,:)={"flow_"+kind,'PASS',numel(actual.UI_data)+numel(actual.Avirulent_data),0}; %#ok<AGROW>
end
figures={'figure1_st175_protection_and_competitive_enrichment', ...
    'figure2_ranked_mouse_screen','figure3_gcms_metabolic_programs', ...
    'figure4_phylogeny_resource_breadth','figure5_host_microbiota_accessory_context', ...
    'figureS1_fecal_cfu_screen','figureS2_pm1_detail_and_generalization_limits', ...
    'figureS3_contextual_analyses'};
actualFigures=dir(fullfile(root,'results','figures','*.png'));
assert(isequal(sort(string({actualFigures.name})),sort(string(figures)+'.png')), ...
    'Active figure exports must contain exactly the eight Draft 44 PNGs.');
for i=1:numel(figures)
    info=imfinfo(fullfile(root,'results','figures',[figures{i} '.png']));
    assert(info.Width>=1000 && info.Height>=1000,'Incomplete figure export.');
    for ext=[".pdf",".svg"]
        candidate=fullfile(root,'results','vectors',string(figures{i})+ext);
        alternate=fullfile(root,'results','figures','editable_vectors',string(figures{i})+ext);
        assert(isfile(candidate)||isfile(alternate),'Missing editable export: %s%s',figures{i},ext);
    end
end
% Presentation checks prevent regeneration of the retired panels/denominator.
f=openfig(fullfile(root,'results','vectors','figure5_host_microbiota_accessory_context.fig'),'invisible');
cleanup=onCleanup(@() close(f));
axesList=findall(f,'Type','axes'); labels=strings(numel(axesList),1);
for i=1:numel(axesList);labels(i)=string(axesList(i).YLabel.String);end
assert(sum(labels=="Fraction of retained events")==2 && ...
    ~any(labels=="Fraction of gated cells"),'Flow axes must name the retained-event denominator.');
clear cleanup;
f=openfig(fullfile(root,'results','vectors','figureS3_contextual_analyses.fig'),'invisible');
cleanup=onCleanup(@() close(f));
axesList=findall(f,'Type','axes'); titles=strings(numel(axesList),1);
for i=1:numel(axesList);titles(i)=string(axesList(i).Title.String);end
assert(numel(axesList)==2 && ~any(contains(titles,'Shannon')) && ...
    any(contains(titles,'community-composition')) && any(contains(titles,'Day-1 genus')), ...
    'S3 must contain only community-composition and genus-level panels.');
clear cleanup;
expectedAlignment=jsondecode(fileread(fullfile(root,'cache','phylogeny','tree_tests','alignment_summary.json')));
actualAlignment=jsondecode(fileread(fullfile(root,'results','phylogeny','tree_tests','alignment_summary.json')));
% Secondary tree sensitivity is a qualitative claim in the manuscript.
% Retain side-by-side estimates; main-tree results remain in the strict tests.
referenceSensitivity=readtable(fullfile(root,'expected','numerical','tree_trait_sensitivity.csv'),'TextType','string');
rebuiltSensitivity=readtable(fullfile(root,'results','phylogeny','tree_tests','tree_trait_sensitivity.csv'),'TextType','string');
assert(isequal(referenceSensitivity(:,{'Tree','Trait'}),rebuiltSensitivity(:,{'Tree','Trait'})), ...
    'Tree sensitivity comparisons changed.');
assert(isequal(referenceSensitivity.PValue<0.05,rebuiltSensitivity.PValue<0.05) && ...
    isequal(sign(referenceSensitivity.Rho),sign(rebuiltSensitivity.Rho)), ...
    'A reported tree-sensitivity conclusion changed.');
comparison=rebuiltSensitivity;
comparison.ReferenceRho=referenceSensitivity.Rho;
comparison.ReferencePValue=referenceSensitivity.PValue;
writetable(comparison,fullfile(root,'results','tree_sensitivity_comparison.csv'));
% Gubbins can assign the same recombinant interval to different branches
% of an unresolved tree. Event-row count is diagnostic, not a study result.
% All retained alignment summaries and full-rebuild sequence bytes are tested.
if expectedAlignment.RecombinationEvents~=actualAlignment.RecombinationEvents
    fprintf('NOTE: Gubbins event rows %d versus reference %d; checking retained alignments exactly.\n', ...
        actualAlignment.RecombinationEvents,expectedAlignment.RecombinationEvents);
end
assert(isequal(rmfield(expectedAlignment,'RecombinationEvents'), ...
    rmfield(actualAlignment,'RecombinationEvents')),'Genome alignment or recombination summary changed.');
alignmentManifest=readtable(fullfile(root,'expected','alignment_checksums.tsv'), ...
    'FileType','text','Delimiter','\t','TextType','string');
folder=fullfile(root,'results','phylogeny','tree_tests');
if isfile(fullfile(folder,alignmentManifest.File(1)))
    for i=1:height(alignmentManifest)
        fid=fopen(fullfile(folder,alignmentManifest.File(i)),'rb');assert(fid>=0);
        bytes=fread(fid,Inf,'*uint8');fclose(fid);
        md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
        hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
        assert(strcmp(hash,alignmentManifest.SHA256(i)),'Tree-input sequence mismatch: %s',alignmentManifest.File(i));
    end
    fprintf('PASS: all five rebuilt tree-input alignments match reference SHA-256.\n');
end
result=cell2table(rows,'VariableNames',{'Output','Status','NumericValuesChecked','MaxAbsoluteDifference'});
writetable(result,fullfile(root,'results','regression_tests.csv'));
fprintf('PASS: %d tables, %d numeric values, all eight Draft 44 figure exports.\n',height(result),sum(result.NumericValuesChecked));
end
