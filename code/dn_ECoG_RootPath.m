function rootPath = dn_ECoG_RootPath()
% function rootPath = temporalRootPath()

rootPath = which('dn_ECoG_RootPath');
rootPath = fileparts(fileparts(rootPath));

end