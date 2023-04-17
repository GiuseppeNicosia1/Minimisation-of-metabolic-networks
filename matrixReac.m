function matrixReac(model, ncores)

% load('totalResults.mat');
load(['solutions' model filesep 'minimalSolutions.mat']);
totalResults = minSol;
% addpath('classes')
% load('solutionsfbamodel_yeast7.6_SD_PLOS\MinCell_Gene\GDMOobj.mat')
fbamodel = load(['solutions' model filesep 'fbamodel.mat']);
fbamodel = fbamodel.fbamodel;
[m,~] = size(totalResults);
matr = zeros(m);
index = true(fbamodel.nrxn,1);
for j = 1: fbamodel.nrxn
    if isempty(fbamodel.optPts{j})
        index(j) = false;
    end
end
optPts = fbamodel.optPts;
present = true(fbamodel.nrxn,m);

delete(gcp('nocreate'))
parpool('local',ncores)
parfor i = 1:size(totalResults,1)
    x = ~logical(totalResults(i,1:fbamodel.ngenes));
    y = true(fbamodel.ngenesBU,1);
    y(~fbamodel.essentialGenes) = x;
    x = y;
    present(:,i) = evalSol(x,optPts);
end
% N = nnz(index);
for i = 1:m

    present1 = present(:,i) & index;
    parfor j = 1:i-1
        
        present2 = present(:,j) & index;
%         matr(i,j) = mean([nnz(present1 & present2)/nnz(present1) nnz(present1 & present2)/nnz(present2)]);
        matr(i,j) = nnz(present1 & present2)/nnz(present1 | present2);
    end
end
save(['solutions' model filesep 'matrReac.mat'],'matr')

% figure
% histogram(matr(matr~=0),20,'Normalization','Probability')
       
