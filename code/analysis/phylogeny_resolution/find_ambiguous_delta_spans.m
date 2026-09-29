function intervals=find_ambiguous_delta_spans(path)
% Identify equally best gap placements for the same alignment endpoints.
% delta-filter can keep or drop tied alternatives depending on input order.
% Exclude this ambiguity explicitly; do not exclude a uniquely better match.
lines=splitlines(string(fileread(path)));keys=strings(0,1);spans=zeros(0,2);
gaps=strings(0,1);errors=zeros(0,1);query="";i=3;
while i<=numel(lines)
    line=strtrim(lines(i));i=i+1;
    if strlength(line)==0,continue;end
    if startsWith(line,'>'),query=line;continue;end
    v=sscanf(line,'%d');assert(numel(v)==7);
    gap="";
    while str2double(lines(i))~=0,gap=gap+" "+lines(i);i=i+1;end
    i=i+1;
    if abs(v(2)-v(1))+1<500,continue;end
    keys(end+1,1)=query+" "+join(string(v(1:4)')," "); %#ok<AGROW>
    spans(end+1,:)=[min(v(1:2)),max(v(1:2))]; %#ok<AGROW>
    gaps(end+1,1)=gap; %#ok<AGROW>
    errors(end+1,1)=v(5); %#ok<AGROW>
end
[uniqueKeys,~,group]=unique(keys);ambiguous=false(size(keys));
for i=1:numel(uniqueKeys)
    rows=group==i;
    rows=rows & errors==min(errors(rows));
    if nnz(rows)>1 && numel(unique(gaps(rows)))>1,ambiguous(rows)=true;end
end
intervals=unique(spans(ambiguous,:),'rows');
end
