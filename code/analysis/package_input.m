function path=package_input(name)
% Resolve a unique supplied source, excluding validation-only workbooks.
items=dir(fullfile(package_root(),'input','**',char(name)));
assert(numel(items)==1,'Expected one input named %s; found %d.',name,numel(items));
path=fullfile(items.folder,items.name);
end
