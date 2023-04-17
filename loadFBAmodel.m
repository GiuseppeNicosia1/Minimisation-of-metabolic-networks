function FBAobj = loadFBAmodel( model, flagEssGenes, flagEssCouples, flagCOBRA )


model=[model '.mat'];
T = load(model);
%Ensuring that the FBA model loaded is called fbamodel
Ts = fieldnames(T);
FBAstruct = T.(Ts{1});
FBAstruct.ngenes = length(FBAstruct.genes);
if ~flagEssGenes
    FBAstruct.genesBU = FBAstruct.genes;
%     FBAstruct.ngenesBU = FBAstruct.ngenes;
    FBAstruct.essentialGenes = false(FBAstruct.ngenes,1);
elseif FBAstruct.ngenes == length(FBAstruct.essentialGenes)
    FBAstruct.genesBU = FBAstruct.genes;
    FBAstruct.genes = FBAstruct.genes(~FBAstruct.essentialGenes);
    FBAstruct.ngenes = length(FBAstruct.genes);
end
if ~flagEssCouples
    FBAstruct.essentialCouples = sparse(false(FBAstruct.ngenes));
end
if flagCOBRA
    FBAstruct.vmin = FBAstruct.lb;
    FBAstruct.vmax = FBAstruct.ub;
    FBAstruct.f = FBAstruct.c;
    FBAstruct.optPts = FBAstruct.rules;
end
FBAstruct.flagEssGenes = flagEssGenes;
FBAstruct.flagEssCouples = flagEssCouples;
FBAobj = FBAobj_from_FBAstruct( FBAstruct );
end

