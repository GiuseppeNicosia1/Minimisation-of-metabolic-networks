function matrixReacParsimonious(model,ncores)

% Script per il confronto tra le reazioni effettivamente attive usando la
% parsimonious flux balance analysis. L'obiettivo è quello di creare una
% campana come quella per le reazioni potenziali, ma stavolta uilizzano le
% reazioni effettivcamente attive

fbamodel = load(['solutions' model filesep 'fbamodel.mat']);
fbamodel = fbamodel.fbamodel;
load(['solutions' model filesep 'minimalSolutions.mat']);
totalResults = minSol;
[m,~] = size(totalResults);
matr = zeros(m);
delete(gcp('nocreate'))
parpool('local',ncores)
%     optPts = fbamodel.optPts;
    present = true(fbamodel.nrxn,m);
    vMinSol = zeros(m,fbamodel.nrxn);
    parfor i = 1:size(totalResults,1)
        [f, v] = fitnessFunctionParsBoundsSA(totalResults(i,1:fbamodel.ngenes), fbamodel);
        present(:,i) = v~=0;
        vMinSol(i,:) = v';
    end
    for i = 1:m
%         clc
%         fprintf('%.2f %% completed\n',100*i/m)
        present1 = present(:,i);
        parfor j = 1:i-1
            present2 = present(:,j);
            matr(i,j) = nnz(present1 & present2)/nnz(present1 | present2);
        end
    end
    save(['solutions' model filesep 'parsMatrReac.mat'],'matr')
    save(['solutions' model filesep 'vMinSolPars.mat'],'vMinSol')
%     figure
    % histogram(matr(matr~=0),32,'Normalization','Probability')
    % edges = 0.625:0.025:1;
    % histogram(matr(matr~=0),edges,'Normalization','Probability')
    figure
    histogram(matr(matr~=0),20,'Normalization','Probability')

% rmpath(folder)
