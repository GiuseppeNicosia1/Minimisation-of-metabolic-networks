function FBAobj = FBAobj_from_FBAstruct( FBAstruct )

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

FBAobj = FBAmodel(FBAstruct.S,FBAstruct.vmin,FBAstruct.vmax,FBAstruct.f,FBAstruct.genes,FBAstruct.optPts,FBAstruct.genesBU,FBAstruct.essentialGenes,FBAstruct.essentialCouples,FBAstruct.flagEssGenes,FBAstruct.flagEssCouples, FBAstruct.sense);

end

