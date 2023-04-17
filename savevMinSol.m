function savevMinSol(model,ncores)

% Script per il calcolo dei flussi

fbamodel = load(['solutions' model filesep 'fbamodel.mat']);
fbamodel = fbamodel.fbamodel;
load(['solutions' model filesep 'minimalSolutions.mat']);
totalResults = minSol;
[m,~] = size(totalResults);
vMinSol = zeros(m,fbamodel.nrxn);
% vMinSolPars = zeros(m,fbamodel.nrxn);
% delete(gcp('nocreate'))
% parpool('local',ncores)
parfor i = 1:size(totalResults,1)
    [f, v] = fitnessFunctionBoundsSA(totalResults(i,1:fbamodel.ngenes), fbamodel);
    vMinSol(i,:) = v';
end
save(['solutions' model filesep 'vMinSol.mat'],'vMinSol')
% vMinSol = vMinSolPars;
% save('vMinSolPars.mat','vMinSol')
