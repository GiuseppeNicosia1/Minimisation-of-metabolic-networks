function findMinimal(model, ncores, maxred)

load(['solutions' model filesep 'totalChromosome.mat'])
totalChromosome = totalChromosome(:,1:end-2);
[N,M] = size(totalChromosome);
totalChromosome = sortrows(totalChromosome,-M);
index = true(N,1);
V = M - 2;
if ~exist(['solutions' model filesep 'minimalSolutionsOld.mat'],'file')
    for i = 1:N
%         clc
%         fprintf('%.2f %% completed\n',100*i/N);
        x = totalChromosome(i,:);
        nKO = x(M);
    %     t = totalChromosome(:,M) == nKO + 1;
    %     temp = totalChromosome(t,:);
    %     D = pdist2(x(1:V),temp(:,1:V),'cityblock');
    %     if any(D == 1)
    %         index(i) = false;
    %     end
        t = ((totalChromosome(:,M) == nKO - 1) & index);
        temp = totalChromosome(t,:);
        D = pdist2(x(1:V),temp(:,1:V),'cityblock');
        I = D == 1;
    %     keyboard
        if any(I)
            t = find(t);
            index(t(I)) = false;
        end
    end

    minSol = totalChromosome(index,:);
    save(['solutions' model filesep 'minimalSolutionsOld.mat'],'minSol')
else
    load(['solutions' model filesep 'minimalSolutionsOld.mat']);
end

delete(gcp('nocreate'))
parpool('local',ncores)

temp = load(['solutions' model filesep 'fbamodel.mat']);
fbamodel = temp.fbamodel;
f = fitnessFunctionBoundsSA(zeros(1,fbamodel.ngenes),fbamodel);
utopianBio = f(1);
% load('utopianBio.mat');
V = fbamodel.ngenes;
M = size(minSol,1);
I = true(M,1);
parfor i = 1:M
%     tic
%     if ~isempty(exhaustive1LevelBreak(minSol(i,1:V), fbamodel, utopianBio))
    if ~isempty(exhaustive1LevelBreak(minSol(i,1:V), fbamodel, utopianBio, maxred))
        I(i) = false;
    end
%     clc
%     fprintf('%.2f %% completed\n',100*i/M);
%     toc
end
minSol = minSol(I,:);
save(['solutions' model filesep 'minimalSolutions.mat'],'minSol')

% t = 1:fbamodel.ngenes;
% for i = 1:size(minSol,1)
%     t = intersect(t,find(minSol(i,1:fbamodel.ngenes)),'stable');
% end

