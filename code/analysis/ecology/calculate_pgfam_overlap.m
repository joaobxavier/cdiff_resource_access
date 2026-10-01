function results = calculate_pgfam_overlap()
%CALCULATE_PGFAM_OVERLAP Spragge's metric for 21 one-strain residents.
% Use the user-supplied VPI10463 PGFam export (1496.1555, MUJV01 assembly).
% This is the alternative VPI assembly in the paper, not the main-tree genome.
% This function calculates family sets; analyze_resource_overlap tests their
% associations with protection.

projectRoot = package_root();
outputRoot = fullfile(projectRoot,'results','analyses','resource_overlap');
if ~isfolder(outputRoot), mkdir(outputRoot); end
vpiFile = package_input('VPI10463_BVBRC_features.csv');
st1File = package_input('ST1_BVBRC_features.csv');
mappingFile = package_input('strain_names_mapping.csv');
assemblyFile = package_input('GCA_001995155.1_VPI10463.fasta');
inputFiles = string({vpiFile; st1File; mappingFile; assemblyFile});
before = hashInputs(inputFiles);
vpi = readStrings(vpiFile);
st1 = readStrings(st1File);
mapping = readStrings(mappingFile);
mapping = mapping(startsWith(mapping.Var2, "ST1-"), :);
required = {'Genome', 'Genome ID', 'Accession', 'Feature ID', 'Annotation', ...
    'Feature Type', 'Start', 'End', 'PATRIC cross-genus families (PGfams)'};
assert(all(ismember(required, vpi.Properties.VariableNames)));
assert(all(ismember(required, st1.Properties.VariableNames)));
assert(all(vpi.("Genome ID") == "1496.1555"), 'Unexpected VPI genome ID.');
assert(all(vpi.Genome == "Clostridioides difficile strain VPI 10463"));
assert(all(vpi.Annotation == "PATRIC") && all(st1.Annotation == "PATRIC"));
assert(numel(unique(vpi.("Feature ID"))) == height(vpi), 'Duplicate VPI feature IDs.');
assert(numel(unique(st1.("Feature ID"))) == height(st1), 'Duplicate ST1 feature IDs.');
assert(height(mapping) == 21 && numel(unique(mapping.Var1)) == 21);
assert(isequal(sort(unique(st1.("Genome ID"))), sort(mapping.Var1)));

% Verify the exported contig names and coordinates against the local FASTA.
% This verifies assembly-coordinate compatibility, not a new genome assembly.
lines = readlines(assemblyFile);
headerRows = find(startsWith(lines, '>'));
stops = [headerRows(2:end) - 1; numel(lines)];
accessions = strings(numel(headerRows), 1);
lengths = zeros(numel(headerRows), 1);
for k = 1:numel(headerRows)
    token = regexp(char(lines(headerRows(k))), '^>(\S+)', 'tokens', 'once');
    assert(~isempty(token));
    accessions(k) = regexprep(string(token{1}), '\.\d+$', '');
    lengths(k) = sum(strlength(strip(lines(headerRows(k)+1:stops(k)))));
end
assert(isequal(sort(accessions), ["MUJV01000001"; "MUJV01000002"]));
assert(isequal(sort(unique(vpi.Accession)), sort(accessions)));
[known, loc] = ismember(vpi.Accession, accessions);
assert(all(known));
starts = str2double(vpi.Start);
ends = str2double(vpi.End);
valid = isfinite(starts) & isfinite(ends) & starts >= 1 & ends >= 1 & ...
    starts == fix(starts) & ends == fix(ends) & max(starts, ends) <= lengths(loc);
assert(all(valid), 'Exported VPI coordinates exceed their assembly.');

familyColumn = 'PATRIC cross-genus families (PGfams)';
vpiCDS = vpi.("Feature Type") == "CDS";
vpiPG = vpi.(familyColumn);
vpiAssigned = vpiCDS & ~ismissing(vpiPG) & vpiPG ~= "";
assert(all(startsWith(vpiPG(vpiAssigned), "PGF_")));
vpiFamilies = unique(vpiPG(vpiAssigned));
st1CDS = st1.("Feature Type") == "CDS";
st1PG = st1.(familyColumn);
st1Assigned = st1CDS & ~ismissing(st1PG) & st1PG ~= "";
assert(all(startsWith(st1PG(st1Assigned), "PGF_")));
template = struct('Strain', "", 'GenomeID', "", 'VPIAssembly', "GCA_001995155.1", ...
    'VPIGenomeID', "1496.1555", 'ST1UniquePGFams', 0, 'VPIUniquePGFams', 0, ...
    'SharedPGFams', 0, 'VPIOnlyPGFams', 0, 'ST1OnlyPGFams', 0, ...
    'VPIOverlapFraction', 0);
rows = repmat(template, height(mapping), 1);
for k = 1:height(mapping)
    selected = st1Assigned & st1.("Genome ID") == mapping.Var1(k);
    residentFamilies = unique(st1PG(selected));
    r = pgfam_pathogen_coverage(vpiFamilies, residentFamilies);
    % Independent direct-set check protects the pathogen denominator.
    shared = intersect(vpiFamilies, residentFamilies);
    assert(r.CoveredFamilies == numel(shared));
    assert(r.OverlapFraction == numel(shared) / numel(vpiFamilies));
    rows(k).Strain = mapping.Var2(k);
    rows(k).GenomeID = mapping.Var1(k);
    rows(k).ST1UniquePGFams = r.ResidentFamilies;
    rows(k).VPIUniquePGFams = r.PathogenFamilies;
    rows(k).SharedPGFams = r.CoveredFamilies;
    rows(k).VPIOnlyPGFams = r.UncoveredFamilies;
    rows(k).ST1OnlyPGFams = r.ResidentFamilies - r.CoveredFamilies;
    rows(k).VPIOverlapFraction = r.OverlapFraction;
end
results = sortrows(struct2table(rows), 'Strain');
assert(all(results.SharedPGFams + results.VPIOnlyPGFams == results.VPIUniquePGFams));
assert(all(results.SharedPGFams + results.ST1OnlyPGFams == results.ST1UniquePGFams));
qc = table("VPI10463", "1496.1555", "GCA_001995155.1", height(vpi), sum(vpiCDS), ...
    sum(vpiAssigned), sum(vpiCDS) - sum(vpiAssigned), numel(vpiFamilies), ...
    sum(vpiAssigned) / sum(vpiCDS), numel(accessions), all(valid), ...
    'VariableNames', {'Strain', 'GenomeID', 'Assembly', 'Features', 'CDS', ...
    'AssignedCDS', 'UnassignedCDS', 'UniquePGFams', 'AssignedCDSFraction', ...
    'ContigsWithFeatures', 'AllFeatureCoordinatesValid'});
after = hashInputs(inputFiles);
assert(isequal(before, after), 'A supplied source file changed during analysis.');
writetable(results, fullfile(outputRoot, 'ST1_VPI1496_1555_PGFam_overlap.csv'));
writetable(qc, fullfile(outputRoot, 'VPI1496_1555_export_QC.csv'));
writetable(table(vpiFamilies, 'VariableNames', {'PGFam'}), ...
    fullfile(outputRoot, 'VPI1496_1555_unique_assigned_PGFams.csv'));
writetable(table(accessions, lengths, 'VariableNames', {'Accession', 'SequenceLength'}), ...
    fullfile(outputRoot, 'VPI1496_1555_assembly_contig_check.csv'));
writetable(before, fullfile(outputRoot, 'ST1_VPI1496_1555_input_hashes_before.csv'));
writetable(after, fullfile(outputRoot, 'ST1_VPI1496_1555_input_hashes_after.csv'));
summary = struct('VPIGenomeID', '1496.1555', 'VPIAssembly', 'GCA_001995155.1', ...
    'VPIAssemblyRole', 'PGFam reference and alternative phylogenetic outgroup; not the primary outgroup', ...
    'ExportSHA256', char(before.SHA256(1)), ...
    'VPIUniqueAssignedPGFams', numel(vpiFamilies), 'ST1Genomes', height(results), ...
    'MinimumOverlapFraction', min(results.VPIOverlapFraction), ...
    'MaximumOverlapFraction', max(results.VPIOverlapFraction), ...
    'DistinctOverlapCounts', numel(unique(results.SharedPGFams)), ...
    'InputFilesUnchanged', true, 'ProtectionAssociations', 'Calculated separately by analyze_resource_overlap', ...
    'MainTreeAssemblyReplaced', false, ...
    'AnnotationSnapshot', 'Date/version not supplied in the feature CSVs');
fid = fopen(fullfile(outputRoot, 'ST1_VPI1496_1555_summary.json'), 'w');
assert(fid >= 0);
fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(summary, 'PrettyPrint', true));
disp(qc);
disp(results(:, {'Strain', 'SharedPGFams', 'VPIUniquePGFams', 'VPIOverlapFraction'}));
end

function data = readStrings(file)
options = detectImportOptions(file, 'VariableNamingRule', 'preserve');
options = setvartype(options, options.VariableNames, 'string');
data = readtable(file, options);
end

function hashes = hashInputs(files)
sha256 = strings(numel(files), 1);
for k = 1:numel(files)
    fid = fopen(files(k), 'rb');
    assert(fid >= 0);
    bytes = fread(fid, Inf, '*uint8');
    fclose(fid);
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(bytes);
    sha256(k) = lower(string(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, [])));
end
[~,names,ext] = arrayfun(@fileparts,files);
hashes = table(names+ext, sha256, 'VariableNames', {'File', 'SHA256'});
end
