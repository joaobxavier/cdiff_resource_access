function root=package_root()
% Resolve this release without relying on the current working directory.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
