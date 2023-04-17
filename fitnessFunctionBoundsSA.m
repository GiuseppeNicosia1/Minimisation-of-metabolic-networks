function [f, v] = fitnessFunctionBoundsSA(x, fbamodel )
    fbamodel.present = true(fbamodel.nrxn,1);
    x = ~x;
    if fbamodel.flagEssGenes
        y = true(fbamodel.ngenesBU,1);
        y(~fbamodel.essentialGenes) = x;
        x = y;
    end
    for i = 1: fbamodel.nrxn
%         if any(fbamodel.rxnGeneMat(i,~x))
            if ~isempty(fbamodel.optPts{i})
                if ~eval(fbamodel.optPts{i})
                    fbamodel.vmin(i) = 0;
                    fbamodel.vmax(i) = 0;
                end
            end
%         end
    end
    %[vbiomass, ~] = fbamodel.fluxBalance();
    [vbiomass, v] = fluxBalanceBoundsSA(fbamodel);
    f(1) = vbiomass; %odd indexes are the various biomasses
    f(2) = sum(~x(1:fbamodel.ngenesBU)); %evene indexes are all the same in this case, i.e. the KOs perfermed
end
