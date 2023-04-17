function saveTotalChromosome(model, optMode, gen, ncores)
% results_for_table

% gen = 5000;
% formatSpec = '%f %f\n';
count = 0;
% load(['solutionsfbamodel_yeast7.6NewSD' filesep 'MinCell_Gene' filesep 'solution1']);
load(['solutions' model filesep optMode filesep 'solution1']);
m = size(chromosome,2)-1;
totalChromosome1 = false(gen*100,m-4);
totalChromosome2 = zeros(gen*100,4);
for j = 1:gen
    load(['solutions' model filesep optMode filesep 'solution' num2str(j)]);
    totalChromosome1(count+1:count+size(chromosome,1),:) = logical(chromosome(:,1:m-4));
    totalChromosome2(count+1:count+size(chromosome,1),:) = [-chromosome(:,m-3:m-2) j*ones(size(chromosome,1),1) (1:size(chromosome,1))'];
    count = count + size(chromosome,1);
end
totalChromosome1 = totalChromosome1(1:count,:);
totalChromosome2 = totalChromosome2(1:count,:);
[totalChromosome1,I] = unique(totalChromosome1,'rows','stable');
totalChromosome2 = totalChromosome2(I,:);
[totalChromosome2,I] = sortrows(totalChromosome2,[1 2]);
totalChromosome1 = totalChromosome1(I,:);
totalChromosome = [totalChromosome1 totalChromosome2];
save(['solutions' model filesep 'totalChromosome.mat'],'totalChromosome','-v7.3');
index = true(size(totalChromosome2,1),1);

delete(gcp('nocreate'))
parpool('local',ncores)
parfor i = 1 :size(totalChromosome2,1)
    for j = i+1:size(totalChromosome2,1)
        if totalChromosome2(i,2)<=totalChromosome2(j,2)
            index(i) = false;
            break
        end
    end
end
totalFront = [totalChromosome1(index(:),:) totalChromosome2(index(:),:)];
save(['solutions' model filesep 'totalFront.mat'],'totalFront');
