function [vbiomass , v] = invokeFBAobject(GDMOobj)

fbamodel = GDMOobj.fbamodel;
%     if(length(bioIdx)==1)
%         if(flagMinBiomass)
%             fbamodel.vmin(bioIdx) =  GDMOobj.min_biomass(1);
%         end
%     end
    if GDMOobj.flagFBA
        [vbiomass , v] = fbamodel.fluxBalance();
    elseif GDMOobj.flagpFBA || (GDMOobj.flagMOMA && GDMOobj.currGen == 1)
        [vbiomass , ~] = fbamodel.fluxBalance();
        v = pFluxBalance(fbamodel, vbiomass);
    elseif GDMOobj.flagMOMA && GDMOobj.currGen > 1
        [vbiomass,v] = MOMAfun(fbamodel, GDMOobj.wildFluxes);
    end
end
