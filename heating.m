function chromosome = heating(GDMOobj)
chromosome = GDMOobj.chromosome;
lastChromosome = GDMOobj.lastChromosome;
V = GDMOobj.V;
M = GDMOobj.M;
utopianBio = GDMOobj.utopianBio;
[n, ~] = size(chromosome);
% chromosome = [chromosome zeros(n,1)];
index = false(n,1);
R = rand(n,1);
P = zeros(n,1);
for i = 1:n
    [t, I] = ismember(chromosome(i,1:V),lastChromosome(:,1:V),'rows');
    if t
        chromosome(i,V+M+3) = lastChromosome(I,V+M+3)+1;
    end
    P(i) = 1/(1+exp(10-chromosome(i,V+M+3)/10));
    if R(i) < P(i)
        index(i) = true;
    end
end

if any(index)
    temp = chromosome(index,:);
    [m, ~] = size(temp);
    parfor j = 1:m
        temp(j,:) = heatingProcedure(temp(j,:),utopianBio,GDMOobj,V,M);
%         isFeasible = false;
%         backup = temp(j,:);
%         while ~isFeasible
%             to = sum(temp(j,1:V));
%             n_knockin=randi([0 to],1);
%             geni_knockout=find(temp(j,1:V)==1);
%             knockin = randi([1 length(geni_knockout)],n_knockin,1);
%             temp(j,geni_knockout(knockin))=0;
%             temp(j,(V+1):(V+M)) = fitnessFunction(temp(j,1:V),GDMOobj);
%             isFeasible = (1+temp(j,V+1)/utopianBio)<0.01;
%             if ~isFeasible
%                 temp(j,:) = backup;
%             end
%         end
    end
    chromosome(index,:) = temp;
    chromosome = non_domination_sort_mod(chromosome, M, V);
end
end

function x = heatingProcedure(x, utopianBio, GDMOobj, V, M)
isFeasible = false;
x(length(x)) = 0;
backup = x;
while ~isFeasible
    to = nnz(x(1:V));
    n_knockin=randi([0 to],1);
    geni_knockout=find(x(1:V)==1);
%     knockin = randi([1 length(geni_knockout)],n_knockin,1);
    knockin = randperm(length(geni_knockout),n_knockin);
    x(geni_knockout(knockin))=0;
    x((V+1):(V+M)) = fitnessFunction(x(1:V),GDMOobj);
    isFeasible = (1+x(V+1)/utopianBio)<GDMOobj.maxred;
    if ~isFeasible
        x = backup;
    end
end
return
end