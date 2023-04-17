function explStrain = exhaustive1LevelBreak(strain, fbamodel, utopianBio, maxred)
% disp('exh call')
explStrain = [];
if length(strain) == fbamodel.ngenes
    if fbamodel.flagEssCouples
        t = feasGenesSA(strain, fbamodel.essentialCouples, fbamodel.ngenes);
    else
        t = find(strain == 0);
    end
    strainBU = strain;
    for i = 1:length(t)
        if t(i) ~= 104
            strain(t(i)) = 1;
            f = fitnessFunctionBoundsSA(strain, fbamodel);
            if (1-f(1)/utopianBio)<maxred
                explStrain = [strain f];
                return
            end
            strain = strainBU;
        end
    end
else
    error('Strain - Genes dimensions mismatch')
end
% disp('exh return')
end