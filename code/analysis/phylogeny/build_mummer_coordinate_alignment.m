function build_mummer_coordinate_alignment
% Whole-genome, reference-coordinate alignment from study assemblies.
% No reference-base filling: all bases are copied from aligned query sequence.
% Only unique one-to-one matches >=500 bp and >=95% identity are retained;
% reference overlaps, ambiguous bases, and +/-5 bp around indels are masked.
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
root=fullfile(package_root(),'results','phylogeny');
work=package_phylogeny_work();bin=getenv('CDIFF_PHYLO_BIN');
ref=fastaread(fullfile(work,'reference.fasta'));refseq=upper(ref.Sequence);
files=dir(fullfile(work,'ST1_*.fasta'));
% macOS can match lower-case Gubbins products with this wildcard on resume.
files=files(~cellfun(@isempty,regexp({files.name},'^ST1_[0-9]+\.fasta$','once')));
assert(numel(files)==21,'Expected exactly 21 numbered ST1 source assemblies.');
names=replace(string({files.name})','.fasta','');
names=[names;"VPI10463";"VPI10463_old"];
seq=repmat('N',numel(names),numel(refseq));
aligned=zeros(numel(names),1);mismatches=aligned;overlaps=aligned;indelmask=aligned;
ambiguousmask=aligned;
for i=1:numel(names)
    query=fullfile(work,names(i)+'.fasta');prefix=fullfile(work,names(i));
    delta=prefix+'.one.delta';
    if ~isfile(delta)
        cmd=sprintf('"%s/nucmer" --maxmatch -t 4 -p "%s" "%s/reference.fasta" "%s" > "%s.nucmer.log" 2>&1',bin,prefix,work,query,prefix);
        [status,msg]=system(cmd);assert(status==0,'%s',msg);
        cmd=sprintf('"%s/delta-filter" -1 -i 95 -l 500 "%s.delta" > "%s"',bin,prefix,delta);
        [status,msg]=system(cmd);assert(status==0,'%s',msg);
    end
    q=fastaread(query);qids=string(regexp({q.Header},'^\S+','match','once'));
    lines=splitlines(string(fileread(delta)));n=3;count=zeros(1,numel(refseq),'uint16');
    observed=repmat('N',1,numel(refseq));mask=false(size(observed));
    queryIndex=[];recordCount=0;
    while n<=numel(lines)
        line=strtrim(lines(n));n=n+1;if strlength(line)==0,continue;end
        if startsWith(line,'>')
            fields=split(extractAfter(line,1));
            assert(fields(1)==string(regexp(ref.Header,'^\S+','match','once')));
            queryIndex=find(qids==fields(2));assert(isscalar(queryIndex));continue;
        end
        v=sscanf(line,'%d');assert(numel(v)==7 && ~isempty(queryIndex));
        assert(v(2)>=v(1));r=v(1);qp=v(3);direction=sign(v(4)-v(3));
        qseq=upper(q(queryIndex).Sequence);recordDiff=0;gapCount=0;ambiguousDiff=0;
        while true
            d=str2double(lines(n));n=n+1;
            if d==0,chunk=v(2)-r+1;else,chunk=abs(d)-1;end
            rr=r:r+chunk-1;qq=qp+direction*(0:chunk-1);
            bases=qseq(qq);
            if direction<0,bases=seqcomplement(bases);end
            observed(rr)=bases;count(rr)=count(rr)+1;
            recordDiff=recordDiff+sum(bases~=refseq(rr) & ismember(bases,'ACGT') & ismember(refseq(rr),'ACGT'));
            ambiguousDiff=ambiguousDiff+sum(~ismember(bases,'ACGT')|~ismember(refseq(rr),'ACGT'));
            r=r+chunk;qp=qp+direction*chunk;
            if d==0,break;end
            mask(max(1,r-5):min(numel(refseq),r+5))=true;
            gapCount=gapCount+1;
            if d>0,r=r+1;else,qp=qp+direction;end
        end
        assert(r==v(2)+1 && qp==v(4)+direction,'Delta endpoint mismatch.');
        % NUCmer does not count ambiguous reference bases as substitutions.
        assert(recordDiff+gapCount==v(5), ...
            'Decoded differences do not match delta errors: %s record %d; mismatch %d, gap %d, ambiguous %d; source %s',names(i),recordCount+1,recordDiff,gapCount,ambiguousDiff,mat2str(v'));
        recordCount=recordCount+1;
    end
    overlaps(i)=sum(count>1);indelmask(i)=sum(mask & count==1);
    % Thread-dependent raw-block order can change delta-filter's handling of
    % equally scored gap placements. Apply the unique-alignment rule explicitly.
    ambiguous=false(size(observed));
    spans=find_ambiguous_delta_spans(prefix+'.delta');
    for j=1:size(spans,1),ambiguous(spans(j,1):spans(j,2))=true;end
    ambiguousmask(i)=sum(ambiguous & count==1 & ~mask);
    keep=count==1 & ~mask & ~ambiguous & ismember(observed,'ACGT') & ismember(refseq,'ACGT');
    seq(i,keep)=observed(keep);aligned(i)=sum(keep);
    mismatches(i)=sum(observed(keep)~=refseq(keep));
    fprintf('%s: %d records, %d callable bases, %d reference SNPs\n',names(i),recordCount,aligned(i),mismatches(i));
    copyfile(delta,fullfile(root,names(i)+'.one.delta'));
end
st1=startsWith(names,'ST1_');
writeFasta(fullfile(root,'mummer_all23.fasta'),names,seq);
writeFasta(fullfile(work,'mummer_st1.fasta'),names(st1),seq(st1,:));
copyfile(fullfile(work,'mummer_st1.fasta'),fullfile(root,'mummer_st1.fasta'));
t=table(names,aligned,mismatches,overlaps,indelmask,ambiguousmask,'VariableNames',{'Strain','CallableBases','ReferenceSNPs','OverlapMaskedBases','IndelFlankMaskedBases','TiedAlignmentMaskedBases'});
writetable(t,fullfile(root,'mummer_callability.csv'));
ska=fastaread(fullfile(root,'st1_k17.aln'));sn=string({ska.Header});
disagreement=zeros(sum(st1),1);common=disagreement;
for i=1:sum(st1)
    a=upper(ska(sn==names(i)).Sequence);b=seq(i,:);
    keep=ismember(a,'ACGT')&ismember(b,'ACGT');
    common(i)=sum(keep);disagreement(i)=sum(a(keep)~=b(keep));
end
writetable(table(names(st1),common,disagreement,'VariableNames',{'Strain','JointCallableBases','DifferentBaseCalls'}),fullfile(root,'ska_mummer_base_concordance.csv'));
save(fullfile(root,'coordinate_alignment.mat'),'names','seq','refseq');
end

function writeFasta(path,names,seq)
fid=fopen(path,'w');assert(fid>=0);c=onCleanup(@()fclose(fid));
for i=1:numel(names),fprintf(fid,'>%s\n%s\n',names(i),seq(i,:));end
end
