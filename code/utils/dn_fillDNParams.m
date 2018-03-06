function outputprm = dn_fillDNParams(inputprm, irfType)

outputprm = [];
%%
switch irfType
    case 'uniphasic'
        outputprm = [inputprm(1), 0, inputprm(2 : 5), 1];
    case 'biphasic'
        outputprm = [inputprm(1 : 3),  2, inputprm(4 : 5), 1];
end
end