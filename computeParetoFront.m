function pfIdxs = computeParetoFront(resObj)
%computeParetoFront, method of the class resObj. It computes the indexes associated to the pareto optimal points of the two
%features input vector fit. 

[N, ~] = size(resObj.fit);
pfIdxs = [];
for ii = 1:N
    candidate = resObj.fit(ii,:);
    notDominated = true;
    for jj = 1:N
        currPoint = resObj.fit(jj,:);
        if(isDominated(candidate , currPoint))
           notDominated = false; 
           break ;
        end
    end
    
    if(notDominated)
        pfIdxs = [pfIdxs , ii];
    end
    
end




end

function bool = isDominated(x,y) 
%is x dominated by y?

bool = all(y <= x) && any(y < x);

end