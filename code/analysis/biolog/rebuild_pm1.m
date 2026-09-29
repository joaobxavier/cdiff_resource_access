function rebuild_pm1()
% Original plate-level coefficient rule; union positive calls across repeats.
root=package_root();out=fullfile(root,'results','tables','biolog');
if ~isfolder(out),mkdir(out);end
files=dir(fullfile(root,'input','biolog','**','*PM1*.xlsx'));
assert(numel(files)==30,'Expected 30 raw PM1 plates.');
mapping=package_input('Biolog_names_sorted.xlsx');
matrix=table(); strains=strings(0,1);
for i=1:numel(files)
    plate=fullfile(files(i).folder,files(i).name);
    growth=analyze_biolog_plate(plate,mapping);
    name=string(regexp(files(i).name,'^(ST1-\d+|VPI)','match','once'));
    name=replace(name,'-','_');
    if i==1,matrix=table(growth.Metabolites,'VariableNames',{'Metabolites'});end
    [ok,ix]=ismember(matrix.Metabolites,growth.Metabolites);assert(all(ok));
    calls=growth.Growth(ix);
    if any(strains==name),matrix.(name)=max(matrix.(name),calls);
    else,matrix.(name)=calls;strains(end+1)=name;end %#ok<AGROW>
    writetable(growth,fullfile(out,replace(files(i).name,'.xlsx','_calls.csv')));
end
assert(sum(matrix{:,2:end},'all')==700,'PM1 reconstruction changed.');
writetable(matrix,fullfile(out,'Biolog_growth_matrix.csv'));
% Supplied combined matrix is used only here, as an independent expectation.
expected=readtable(fullfile(root,'expected','source','Biolog_growth_matrix.xlsx'),'VariableNamingRule','preserve','TextType','string');
[ok,ix]=ismember(string(expected.Metabolites),matrix.Metabolites);assert(all(ok));
for k=2:width(expected)
    name=expected.Properties.VariableNames{k};
    assert(isequal(double(expected.(name)),matrix.(name)(ix)),'Mismatch in %s.',name);
end
fprintf('PASS: all 700 calls match the supplied combined matrix.\n');
end
