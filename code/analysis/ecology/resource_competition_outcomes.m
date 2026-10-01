function outcome = resource_competition_outcomes(overlap,candidateTotal,pathogenTotal,lambdaGrid)
% Binary-access model parameterized by overlap and total resource access.
assert(overlap>=0 && overlap<=min(candidateTotal,pathogenTotal));
assert(all([overlap,candidateTotal,pathogenTotal]==fix([overlap,candidateTotal,pathogenTotal])));
candidateAdditional=candidateTotal-overlap;
pathogenAdditional=pathogenTotal-overlap;
outcome=strings(size(lambdaGrid));
for k=1:numel(lambdaGrid)
    lambda=lambdaGrid(k);
    candidateFeasible=candidateTotal>lambda; pathogenFeasible=pathogenTotal>lambda;
    if candidateFeasible && ~pathogenFeasible, outcome(k)="ST1 excludes VPI";
    elseif pathogenFeasible && ~candidateFeasible, outcome(k)="VPI excludes ST1";
    elseif ~candidateFeasible && ~pathogenFeasible, outcome(k)="neither feasible";
    else
        candidateInvasion=candidateAdditional-lambda*pathogenAdditional/pathogenTotal;
        pathogenInvasion=pathogenAdditional-lambda*candidateAdditional/candidateTotal;
        if candidateInvasion>0 && pathogenInvasion>0, outcome(k)="coexistence";
        elseif candidateInvasion>0 && pathogenInvasion<=0, outcome(k)="ST1 excludes VPI";
        elseif candidateInvasion<=0 && pathogenInvasion>0, outcome(k)="VPI excludes ST1";
        else, outcome(k)="priority";
        end
    end
end
end
