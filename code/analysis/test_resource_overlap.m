function test_resource_overlap()
% Verify every new analysis row, source mapping and model grid outcome.
root=package_root(); out=fullfile(root,'results','analyses','resource_overlap');
t=readtable(fullfile(out,'rank_associations.csv'),'TextType','string');
expected=readtable(fullfile(root,'expected','resource_overlap','rank_associations.csv'),'TextType','string');
assert(isequal(t.Comparison,expected.Comparison));
vars={'Rho','PValue','BH_Q','Bootstrap_CI_Lower','Bootstrap_CI_Upper', ...
    'N','DF','ValidBootstrapDraws','InvalidBootstrapDraws'};
for j=1:numel(vars)
    assert(max(abs(t.(vars{j})-expected.(vars{j})))<1e-11,'Resource overlap regression failed: %s',vars{j});
end
assert(isequal(t.DF,[19;18;17;19;18;19;18;19;18]));
assert(all(t.ValidBootstrapDraws+t.InvalidBootstrapDraws==5000));
assert(all(t.InvalidBootstrapDraws([6 7])==1));
s=readtable(fullfile(out,'strain_aligned_inputs.csv'),'TextType','string');
expectedInputs=readtable(fullfile(root,'expected','resource_overlap','strain_aligned_inputs.csv'),'TextType','string');
assert(isequal(s.Strain,expectedInputs.Strain) && height(s)==21);
assert(isequal(s.GenomeID,expectedInputs.GenomeID));
assert(max(abs(s{:,3:end}-expectedInputs{:,3:end}),[],'all')<1e-11);
assert(max(abs(s.OverlapAll-s.SharedAll/23))<1e-14 && max(abs(s.OverlapAmino-s.SharedAmino/3))<1e-14);
assert(max(abs(s.OverlapNonAmino-s.SharedNonAmino/20))<1e-14);
assert(max(abs(sort(unique(s.OverlapAmino))-[0;1/3;2/3;1]))<1e-14);
loo=readtable(fullfile(out,'leave_one_strain_out.csv'),'TextType','string');
assert(height(loo)==189 && all(loo.Valid));
plates=readtable(fullfile(out,'equal_plate_all_81_selections.csv'),'TextType','string');
assert(height(plates)==486 && all(plates.Valid));
assert(numel(unique(plates.Selection))==81);
native=readtable(fullfile(out,'native_matlab_crosschecks.csv'));
assert(max(abs(native.ImplementedRho-native.NativeRho))<1e-12);
assert(max(abs(native.ImplementedP-native.NativeP))<1e-12);
lambdaGrid=linspace(0,82,501); comparisons=0;
for i=1:height(s)
    shared=s.SharedAll(i); st1Private=s.TotalAccess(i)-shared; vpiPrivate=23-shared;
    current=resource_competition_outcomes(shared,s.TotalAccess(i),23,lambdaGrid);
    old=legacy_outcomes(shared,st1Private,vpiPrivate,lambdaGrid);
    assert(isequal(current,old)); comparisons=comparisons+numel(current);
end
% Boundary cases, including zero breadth and equality to loss.
for totals={[0 0 0],[0 1 1],[1 1 1],[1 2 3],[2 5 2]}
    v=totals{1};
    assert(isequal(resource_competition_outcomes(v(1),v(2),v(3),[0 1 2 3 5]), ...
        legacy_outcomes(v(1),v(2)-v(1),v(3)-v(1),[0 1 2 3 5])));
end
fprintf('PASS: 9 associations, 189 strain deletions, 81 plate selections, %d unchanged model-grid outcomes.\n',comparisons);
end

function outcome=legacy_outcomes(shared,st1Private,vpiPrivate,lambdaGrid)
outcome = strings(size(lambdaGrid));
for k = 1:numel(lambdaGrid)
    lambda = lambdaGrid(k);
    st1Feasible = (shared + st1Private) > lambda;
    vpiFeasible = (shared + vpiPrivate) > lambda;
    if st1Feasible && ~vpiFeasible
        outcome(k) = "ST1 excludes VPI";
    elseif vpiFeasible && ~st1Feasible
        outcome(k) = "VPI excludes ST1";
    elseif ~st1Feasible && ~vpiFeasible
        outcome(k) = "neither feasible";
    else
        st1Invasion = st1Private - ...
            lambda * vpiPrivate / (shared + vpiPrivate);
        vpiInvasion = vpiPrivate - ...
            lambda * st1Private / (shared + st1Private);
        if st1Invasion > 0 && vpiInvasion > 0
            outcome(k) = "coexistence";
        elseif st1Invasion > 0 && vpiInvasion <= 0
            outcome(k) = "ST1 excludes VPI";
        elseif st1Invasion <= 0 && vpiInvasion > 0
            outcome(k) = "VPI excludes ST1";
        else
            outcome(k) = "priority";
        end
    end
end

end
