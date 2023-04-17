function present = evalSol(x,optPts)
n = length(optPts);
present= true(n,1);
for j = 1:n
    if ~isempty(optPts{j})
        present(j) = eval(optPts{j});
    end
end
end