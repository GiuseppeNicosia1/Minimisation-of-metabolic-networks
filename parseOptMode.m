function resObj = parseOptMode(resObj)
%parseOptMode method for the class optResult. It parses the string
%resObj.optMode in order to set the falgs that specify what optimization
%has been performed.
switch resObj.optMode
    case 'MinCell'
        resObj.flagMinCell = true;
        resObj.flagGene = false;
        resObj.flagRedirector = false;
        resObj.flagKO = false;
        resObj.flagWC = false;
    case 'MinCell_Gene'
        resObj.flagMinCell = true;
        resObj.flagGene = true;
        resObj.flagRedirector = false;
        resObj.flagKO = false;
        resObj.flagWC = false;
    case 'Redirector'
        resObj.flagMinCell = false;
        resObj.flagGene = false;
        resObj.flagRedirector = true;
        resObj.flagKO = false;
        resObj.flagWC = false;
    case 'KO'
        resObj.flagMinCell = false;
        resObj.flagGene = false;
        resObj.flagRedirector = false;
        resObj.flagKO = true;
        resObj.flagWC = false;
    case 'Whole-Cell'
        resObj.flagMinCell = false;
        resObj.flagGene = false;
        resObj.flagRedirector = false;
        resObj.flagKO = false;
        resObj.flagWC = true;
    otherwise
        warning('optMode not supported! Going for minimal cell optimization, using gene variables');
        resObj.flagMinCell = true;
        resObj.flagGene = true;
        resObj.flagRedirector = false;
        resObj.flagKO = false;
        resObj.flagWC = false;
end

end

