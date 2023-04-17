function f = non_domination_sort_mod(x, M, V)

%  Prende un cromosoma in input e restituisce un cromosoma in output, con
%  in più il numero di front e la crowding distance (the greater the better)
%% function f = non_domination_sort_mod(x, M, V)

[N, ~] = size(x);

% Initialize the front number to 1.
front = 1;

% There is nothing to this assignment, used only to manipulate easily in
% MATLAB.
F(front).f = [];
individual = struct('n',0,'p',cell(N,1));

%% Non-Dominated sort. 

KOs = sum(x(:,1:V),2)';
for i = 1 : N
    index = true(N,1);
    index(i) = false;
    D = pdist2(x(i,1:V),x(index,1:V),'cityblock');
    KOsT = KOs(index);
    tt = 1:N;
    tt = tt(index);
    % Individuals which this individual dominate
    individual(i).p = tt((KOsT == KOs(i) - 1) & (D == 1));
    % Number of individuals that dominate this individual
    individual(i).n = nnz((KOsT == KOs(i) + 1) & (D == 1));
    if individual(i).n == 0
        x(i,M + V + 1) = 1;
        F(front).f = [F(front).f i];
    end
    individual(i).d = D;
end
% Find the subsequent fronts
while ~isempty(F(front).f)
    Q = [];
    for i = 1 : length(F(front).f)
        if ~isempty(individual(F(front).f(i)).p)
            for j = 1 : length(individual(F(front).f(i)).p)
            	individual(individual(F(front).f(i)).p(j)).n = ...
                	individual(individual(F(front).f(i)).p(j)).n - 1;
                if individual(individual(F(front).f(i)).p(j)).n == 0
               		x(individual(F(front).f(i)).p(j),M + V + 1) = ...
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
        end
    end
    front =  front + 1;
    F(front).f = Q;
end

f = x;
for i = 1:size(f,1);
    % QUESTA è LA DEFINIZIONE DI DISTANZA, CAMBIARE t IN MODO OPPORTUNO
%     Prendo l'inverso del numero di elementi vicini 10 o meno
%     t = nnz(individual(i).d <= 10);
%     if t
%         f(i,V+M+2) = 1/t;
%     else
%         f(i,V+M+2) = Inf;
%     end

%     Prendo la media delle distanze del punto con tutti gli altri
%     t = mean(individual(i).d);
%     f(i,V+M+2) = t;

%     Prendo la media della distanza dei 10 elementi più vicini
%     [~, k] = minK(individual(i).d,10); % R2017b
%     t = mean(individual(i).d(k));
    nel = 10;
    val = individual(i).d;
    dm = zeros(nel,1);
    for j = 1:nel
        [dm(j), tempk] = min(val);
        val(tempk) = Inf;
    end
    t = mean(dm);

    f(i,V+M+2) = t;

end
    
end
