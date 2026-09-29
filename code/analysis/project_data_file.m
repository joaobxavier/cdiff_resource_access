function path=project_data_file(~,varargin)
% Compatibility adapter: legacy study scripts resolve supplied basenames.
path=package_input(varargin{end});
end
