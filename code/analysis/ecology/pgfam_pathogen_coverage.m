function result = pgfam_pathogen_coverage(pathogenFamilies, residentFamilies)
%PGFAM_PATHOGEN_COVERAGE Pathogen-centered overlap of unique assigned PGFams.
% A single resident is the one-member-community case of Spragge et al.'s
% community-union metric. Gene copy numbers and product descriptions do not
% change the weight of a family. Empty/missing assignments are not families.

pathogenFamilies = cleanFamilies(pathogenFamilies);
residentFamilies = cleanFamilies(residentFamilies);
if isempty(pathogenFamilies)
    error('PGFam:EmptyPathogen', 'The pathogen has no assigned PGFam families.');
end
covered = ismember(pathogenFamilies, residentFamilies);
result = struct('PathogenFamilies', numel(pathogenFamilies), ...
    'ResidentFamilies', numel(residentFamilies), ...
    'CoveredFamilies', sum(covered), ...
    'UncoveredFamilies', sum(~covered), ...
    'OverlapFraction', sum(covered) / numel(pathogenFamilies));
end

function families = cleanFamilies(values)
families = string(values(:));
families = families(~ismissing(families) & strlength(families) > 0);
families = unique(families);
end
