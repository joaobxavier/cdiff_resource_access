function check_release_code()
% Run MATLAB's analyzer on every distributed script, retaining all warnings.
root=package_root();files=dir(fullfile(root,'code','**','*.m'));
files=[files;dir(fullfile(root,'run_all.m'))];
rows=cell(0,4);
for k=1:numel(files)
    name=fullfile(files(k).folder,files(k).name);
    issues=checkcode(name,'-id');
    for j=1:numel(issues)
        rows(end+1,:)={erase(name,[root filesep]),issues(j).line,issues(j).id,issues(j).message}; %#ok<AGROW>
    end
end
t=cell2table(rows,'VariableNames',{'File','Line','ID','Message'});
writetable(t,fullfile(root,'results','static_analysis.csv'));
% Parse errors are release blockers; style/performance warnings are retained.
assert(~any(ismember(string(t.ID),["PARSE","SCERR","ENDCT","INVRT"])), ...
    'MATLAB parser errors found; inspect results/static_analysis.csv.');
fprintf('Analyzed %d MATLAB files; %d diagnostics retained.\n',numel(files),height(t));
end
