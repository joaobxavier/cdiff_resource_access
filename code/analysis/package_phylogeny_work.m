function folder=package_phylogeny_work()
% Gubbins helper tools require a whitespace-free working directory.
folder=getenv('CDIFF_PHYLO_WORK');
if isempty(folder)
    folder=getappdata(0,'CDiffPhylogenyWork');
    owner=getappdata(0,'CDiffPhylogenyRoot');
    if isempty(folder) || ~strcmp(owner,package_root())
        folder=tempname; mkdir(folder);setappdata(0,'CDiffPhylogenyWork',folder);
        setappdata(0,'CDiffPhylogenyRoot',package_root());
    end
end
assert(isempty(regexp(folder,'\s','once')),'Set CDIFF_PHYLO_WORK to a new whitespace-free directory.');
if ~isfolder(folder),mkdir(folder);end
end
