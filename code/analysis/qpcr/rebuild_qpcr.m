function rebuild_qpcr()
% Recalculate Kevin's KS65 proportions from original quantity estimates.
% The same mass and 50/2 volume conversion applies to both strains in a
% pellet, so it cancels in their relative proportion. No new calibration fit.
source=package_input('KS65.xlsx');
days=[0.5 1 3];Day=[];Mouse=[];ST175=[];VPI=[];Copies75=[];CopiesVPI=[];
for k=1:3
    sheet="Day "+string(days(k));
    if k==2
        % Kevin's plotted day-1 values use CFU/PBS pellet rows 13:17.
        quantities=readmatrix(source,'Sheet',sheet,'Range','L13:M17');
        masses=readmatrix(source,'Sheet',sheet,'Range','B13:C17');
    else
        quantities=readmatrix(source,'Sheet',sheet,'Range','F4:G8');
        masses=readmatrix(source,'Sheet',sheet,'Range','B4:C8');
    end
    assert(isequal(size(quantities),[5 2]));
    delta=masses(:,2)-masses(:,1);
    copies=25*quantities./delta;
    percentages=100*quantities./sum(quantities,2);
    Day=[Day;repmat(days(k),5,1)];Mouse=[Mouse;(1:5)']; %#ok<AGROW>
    ST175=[ST175;percentages(:,1)];VPI=[VPI;percentages(:,2)]; %#ok<AGROW>
    Copies75=[Copies75;copies(:,1)];CopiesVPI=[CopiesVPI;copies(:,2)]; %#ok<AGROW>
end
result=table(Day,Mouse,ST175,VPI,'VariableNames',{'Day','Mouse','ST1-75','VPI'});
assert(sum(isfinite(ST175(Mouse>=3)))==8);
assert(isnan(ST175(Day==1 & Mouse==5)));
out=fullfile(package_root(),'results','tables','main');
writetable(result,fullfile(out,'qpcr_fractions.csv'));
writetable(table(Day,Mouse,Copies75,CopiesVPI),fullfile(out,'qpcr_genomic_copies_per_gram.csv'));
expected=readtable(fullfile(package_root(),'expected','source','qPCR section.xlsx'),'VariableNamingRule','preserve');
for i=1:height(expected)
    row=Day==expected.Day(i)&Mouse==expected.Mouse(i);
    if any(row) && expected.Mouse(i)>=3
        observed=ST175(row); target=expected.('ST1-75')(i);
        assert((isnan(observed)&&isnan(target)) || abs(observed-target)<1e-8);
    end
end
fprintf('PASS: eight measured mixed-infection qPCR fractions reconstructed; missing day 1 preserved.\n');
end
