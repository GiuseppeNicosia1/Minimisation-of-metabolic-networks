function [bestBioFit, bestKoFit, Close2UtopianFit, bestTradeFit_allCorners] = computeNotablePoints(resObj,pfNormFit,absPfPop)
%computeNotablePoints, method of the class resObj. It computes three
%"notable" solutions from the observed pareto front: 1) the one that has
%the highest biomass production.
%2) the one that has the highest synthetic objective. 3) The best trade-off
%among them.
%Usage: [bestBioFit, bestKoFit, Close2UtopianFit] = computeNotablePoints(resObj,pfNormFit,absPfPop)

[evenIdxs ,oddIdxs ] = resObj.getFFidxs();
kk = 1; %Corner index

bestBioFit = zeros(2, (resObj.M)/2);
bestKoFit = zeros(2,(resObj.M)/2);
Close2UtopianFit = zeros(2,(resObj.M)/2);

bestBioPop = zeros(resObj.V,(resObj.M)/2);
bestKoPop = zeros(resObj.V,(resObj.M)/2);
Close2UtopianPop = zeros(resObj.V,(resObj.M)/2);

for ii = 1:resObj.M
    %best biomass
    if(oddIdxs(ii)) %if odd then it is a biomass index
        [~,idx] = min(pfNormFit(:,ii) );
        bestBioFit(:,kk) = resObj.pfFit(idx,ii:(ii+1));
        %bestBioPop = absPfPop(idx,:);
        %bestBioPop(:,kk) = double(absPfPop{idx});
        bestBioPop(:,kk) = decodeStrain(absPfPop{idx} , resObj );
    end
    
    %best kos
    if(evenIdxs(ii)) %if even then it is a KO index
        [~,idx] = min(pfNormFit(:,ii) );
        bestKoFit(:,kk) = resObj.pfFit(idx,(ii-1):ii);
        %bestKoPop = absPfPop(idx,:);
        %bestKoPop(:,kk) = double(absPfPop{idx});
        bestKoPop(:,kk) = decodeStrain( absPfPop{idx} , resObj );
    end
    
    if(evenIdxs(ii)) %if even than I get ii and ii-1 and I have the fit associated to one corner
        %closest2utopian
        [Npf, ~] = size(pfNormFit);
        PfNormFitNorms = zeros(1,Npf);
        for jj = 1:Npf
            PfNormFitNorms(jj) = norm(pfNormFit(jj,(ii-1):ii)); 
        end
        [~,idx] = min(PfNormFitNorms);
        Close2UtopianFit(:,kk) = resObj.pfFit(idx,(ii-1):ii);
        %Close2UtopianPop = absPfPop(idx,:);
        %Close2UtopianPop(:,kk) = double(absPfPop{idx});
        Close2UtopianPop(:,kk) = decodeStrain( absPfPop{idx} , resObj );
        kk = kk +1; % I am done with this corner, I go for the next one...
    end
end


%Saving results to file

save([resObj.results_folder,'bestBiomassFit.txt'],'bestBioFit','-ascii');
save([resObj.results_folder,'bestBiomassTot.txt'],'bestBioPop','-ascii');


if(resObj.flagMinCell || resObj.flagWC )
    save([resObj.results_folder,'bestKoFit.txt'],'bestKoFit','-ascii');
    save([resObj.results_folder,'bestKoTot.txt'],'bestKoPop','-ascii');    
elseif(resObj.flagKO || resObj.flagRedirector)
    save([resObj.results_folder,'bestSynthFit.txt'],'bestKoFit','-ascii');
    save([resObj.results_folder,'bestSynthTot.txt'],'bestKoPop','-ascii');
end

save([resObj.results_folder,'Close2UtopianFit.txt'],'Close2UtopianFit','-ascii');
save([resObj.results_folder,'Close2UtopianPop.txt'],'Close2UtopianPop','-ascii');
%Finally I compute the best trade Off among all the objective functions.

if(resObj.flagCorner)
    [Npf, ~] = size(pfNormFit);
    PfNormFitNorms = zeros(1,Npf);
    for jj = 1:Npf
        PfNormFitNorms(jj) = norm(pfNormFit(jj,:)); 
    end
    [~,idx] = min(PfNormFitNorms);
    bestTradeFit_allCorners = resObj.pfFit(idx,:);
    %Close2UtopianPop = absPfPop(idx,:);
    bestTradePop_allCorners = double(absPfPop{idx});
    save([resObj.results_folder,'Close2UtopianFit_AllCorners_.txt'],'bestTradeFit_allCorners','-ascii');
    save([resObj.results_folder,'Close2UtopianPop_AllCorners_.txt'],'bestTradePop_allCorners','-ascii');    
else
    bestTradeFit_allCorners = [];
end


end