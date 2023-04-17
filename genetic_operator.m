function f  = genetic_operator(parent_chromosome, GDMOobj)

%% function f  = genetic_operator(parent_chromosome, M, V)
% 
% This function is utilized to produce offsprings from parent chromosomes.
%
% parent_chromosome - the set of selected chromosomes.
% M - number of objective functions
% V - number of decision varaiables
%
% The genetic operation is performed only on the decision variables, that
% is the first V elements in the chromosome vector. 
[N,~] = size(parent_chromosome);
V = GDMOobj.V;
M = GDMOobj.M;
% for each parent are formed "value" offsprings, only the best is choosen
value = GDMOobj.value ; 

C = V-1; % maximum knockout selected
flagEssCouples = GDMOobj.fbamodel.flagEssCouples;
essentialCouples = GDMOobj.fbamodel.essentialCouples;
child_new = zeros(size(parent_chromosome,1),V+M);
if(GDMOobj.ncores>1)
    parent_chromosome = parent_chromosome(:,1:(V+M));
%     par
    for i=1:N
        child=parent_chromosome(i,:);
        backup = child;
        if flagEssCouples
            t = feasGenes(child, essentialCouples, V);
        else
            t = find(child == 0);
        end
%         h_j=randi([1 V],1);
        h_j=t(randi([1 length(t)],1));
        j=1;
        child(j,h_j)=not(child(j,h_j));
        while sum(child(j,1:V))>C
            from = sum(child(j,1:V))-C;
            to = sum(child(j,1:V));
            n_knockin=randi([from to],1);
            geni_knockout=find(child(j,1:V)==1);
            knockin = randi([1 length(geni_knockout)],n_knockin,1);
            child(j,geni_knockout(knockin))=0;  
        end
        child(j,(V+1):(V+M)) = fitnessFunction(child(1,1:V),GDMOobj);
        isFeasible = (1+child(j,V+1)/GDMOobj.utopianBio)<GDMOobj.maxred;
        if ~isFeasible
            child(j,:) = backup;
%             child(j,h_j) = not(child(j,h_j));
        end
        j=2;
        while j<=value && ~isFeasible      
            % random mutations in child
%             h_j=randi([1 V],1);
            h_j=t(randi([1 length(t)],1));
            child(j,:)=child(j-1,:);
            child(j,h_j)=not(child(j,h_j));

            % if the knockout number is greater than C, random genes are turned on
            while sum(child(j,1:V))>C
                from = sum(child(j,1:V))-C;
                to = sum(child(j,1:V));
                n_knockin=randi([from to],1);
                geni_knockout=find(child(j,1:V)==1);
                knockin = randi([1 length(geni_knockout)],n_knockin,1);
                child(j,geni_knockout(knockin))=0;  
            end
            child(j,V+1:V+M) = fitnessFunction(child(j,1:V),GDMOobj);
%             count = count + 1;
%             isFeasible = child(j,V+1);
            isFeasible = (1+child(j,V+1)/GDMOobj.utopianBio)<GDMOobj.maxred;
            if ~isFeasible
                child(j,:) = backup;
%                 child(j,h_j) = not(child(j,h_j));
            end
            j=j+1;
        end% while
%         child = non_domination_sort_mod(child, M, V);
%         child_new(i,:) = replace_chromosome(child, M, V, 1);
        child_new(i,:) = child(j-1,:);
    end
else
    for i=1:N
        child=parent_chromosome(i,:);
        backup = child;
        if flagEssCouples
            t = feasGenes(child, essentialCouples, V);
        else
            t = find(child == 0);
        end
%         h_j=randi([1 V],1);
        h_j=t(randi([1 length(t)],1));
        j=1;
        child(j,h_j)=not(child(j,h_j));
        while sum(child(j,1:V))>C
            from = sum(child(j,1:V))-C;
            to = sum(child(j,1:V));
            n_knockin=randi([from to],1);
            geni_knockout=find(child(j,1:V)==1);
            knockin = randi([1 length(geni_knockout)],n_knockin,1);
            child(j,geni_knockout(knockin))=0;  
        end
        child(j,(V+1):(V+M)) = fitnessFunction(child(1,1:V),GDMOobj);
        isFeasible = (1+child(j,V+1)/GDMOobj.utopianBio)<GDMOobj.maxred;
        if ~isFeasible
            child(j,:) = backup;
%             child(j,h_j) = not(child(j,h_j));
        end
        j=2;
        while j<=value && ~isFeasible         
            % random mutations in child
%             h_j=randi([1 V],1);
            h_j=t(randi([1 length(t)],1));
            child(j,:)=child(j-1,:);
            child(j,h_j)=not(child(j,h_j));

            % if the knockout number is greater than C, random genes are turned on
            while sum(child(j,1:V))>C
                from = sum(child(j,1:V))-C;
                to = sum(child(j,1:V));
                n_knockin=randi([from to],1);
                geni_knockout=find(child(j,1:V)==1);
                knockin = randi([1 length(geni_knockout)],n_knockin,1);
                child(j,geni_knockout(knockin))=0;  
            end
            child(j,V+1:V+M) = fitnessFunction(child(j,1:V),GDMOobj);
%             count = count + 1;
%             isFeasible = child(j,V+1);
            isFeasible = (1+child(j,V+1)/GDMOobj.utopianBio)<GDMOobj.maxred;
            if ~isFeasible
                child(j,:) = backup;
%                 child(j,h_j) = not(child(j,h_j));
            end
            j=j+1;
        end% while
%         child = non_domination_sort_mod(child, M, V);
%         child_new(i,:) = replace_chromosome(child, M, V, 1);
        child_new(i,:) = child(j-1,:);
    end
end
f = child_new(:,1:M+V);
