% Produce shared traits/calls in the original resampling order.
prepare_figure4;
resamplingOrder=[53 62 75 68 12 49 57 11 65 6 23 27 20 69 63 67 58 2 25 19 66];
[ok,ix]=ismember("ST1-"+string(resamplingOrder)',nicheTree.Strain);assert(all(ok));
traits=nicheTree(ix,:);
writetable(traits,fullfile(tableDir,'figure4_phylogeny_aligned_traits_and_breadth.csv'));
calls=table(pm1.Metabolites(rowOrder),substrateClass(rowOrder),vpiCalls(rowOrder), ...
    'VariableNames',{'Substrate','Chemical_Class','VPI10463'});
for i=1:height(traits)
    name=replace(traits.Strain(i),'-','_');
    calls.(name)=st1Calls(rowOrder,ix(i));
end
writetable(calls,fullfile(tableDir,'figure4_pm1_call_matrix_plotted_order.csv'));
amino=sum(double(calls{calls.Chemical_Class=="amino acid",4:end}),1)';
names=["Total PM1 breadth";"Candidate-private calls";"Amino acid breadth"];
features=[traits.Total_Breadth,traits.ST1_Private,amino];rows=cell(0,6);
for j=1:3
    for i=1:height(traits)
        keep=true(height(traits),1);keep(i)=false;
        y=traits.Protection_Effect(keep);d=traits.Mono_Disease_Effect(keep);x=features(keep,j);
        Z=[ones(sum(keep),1),tiedrank(d)];
        yr=tiedrank(y)-Z*(Z\tiedrank(y));xr=tiedrank(x)-Z*(Z\tiedrank(x));
        model=fitlm([d,x],y);
        rows(end+1,:)={names(j),traits.Strain(i),sum(keep),corr(xr,yr),model.Coefficients.Estimate(3),model.Coefficients.pValue(3)}; %#ok<SAGROW>
    end
end
influence=cell2table(rows,'VariableNames',{'Feature','Removed_Strain','N_Strains','Partial_Rank_Rho','OLS_Feature_Estimate','OLS_Feature_P_Value'});
out=fullfile(package_root(),'results','analyses','pm1_models','tables');
if ~isfolder(out),mkdir(out);end
writetable(influence,fullfile(out,'highlighted_scalar_influence_analysis.csv'));
