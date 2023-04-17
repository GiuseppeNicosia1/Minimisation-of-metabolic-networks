function GDMOobj = parseOptMode(GDMOobj)
%parseOptMode method of the class GDMOopt. It parses the string
%GDMOobj.optMode in order to set the optimization flags.
switch GDMOobj.optMode
    case 'FBA'
        GDMOobj.flagFBA = true;
        GDMOobj.flagpFBA = false;
        GDMOobj.flagMOMA = false;
    case 'pFBA'
        GDMOobj.flagFBA = false;
        GDMOobj.flagpFBA = true;
        GDMOobj.flagMOMA = false;
    case 'MOMA'
        GDMOobj.flagFBA = false;
        GDMOobj.flagpFBA = false;
        GDMOobj.flagMOMA = true;  
    otherwise
        warning('optMode not supported! Going for FBA');
        GDMOobj.flagFBA = true;
        GDMOobj.flagpFBA = false;
        GDMOobj.flagMOMA = false;
end

end

