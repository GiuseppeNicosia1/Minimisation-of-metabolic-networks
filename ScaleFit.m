function normFit = ScaleFit(fit,A,B)
%Scales the two objective functions in the [0,1] interval, separatelly
%vec1 = fit(:,1);
%vec2 = fit(:,2);
%vec1 = normalize(vec1);
%vec2 = normalize(vec2);
%normFit = [vec1(:) , vec2(:)];
if nargin < 2
	A = [];
	B = [];
end 
normFit = zeros(size(fit,1),size(fit,2));
for ii=1:size(fit,2)
    fitVec = fit(:,ii);
    if(isempty(A) || isempty(B))
        fitVec = normalize(fitVec);
    else
        fitVec = normalize(fitVec,A(ii),B(ii));
    end
    
    normFit(:,ii) = fitVec(:);
end

end



function normVec = normalize(vec,a,b)
%given an input vector with min=a and max=b , a<b, returns (vec-a)/(b-a),
%i.e. it normalizes the vector entries in the [0,1] interval
if nargin < 2
	a = min(vec);
	b = max(vec);
end
normVec = (vec-a)/(b-a);
end
