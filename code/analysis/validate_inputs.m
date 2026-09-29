function validate_inputs()
root=package_root();t=readtable(fullfile(root,'input_manifest.tsv'),'FileType','text','Delimiter','\t','TextType','string');
for i=1:height(t)
    path=fullfile(root,t.Proposed_package_path(i));
    fid=fopen(path,'rb');assert(fid>=0,'Missing input: %s',path);
    bytes=fread(fid,Inf,'*uint8');fclose(fid);
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
    hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
    assert(strcmp(hash,t.SHA256(i)),'Input checksum mismatch: %s',path);
end
fprintf('PASS: SHA-256 verified for all %d supplied files.\n',height(t));
end
