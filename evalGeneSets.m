function KOGeneSets = evalGeneSets(fbamodel,x)
%evalGeneSets evaluates the field optPts of the FBAmodel object, using x as
%variable. The output, KOGeneSets, is a vector of the same length of fbamodel.pts
% in which 1 means that the corresponding gene set is KO, 0 means that it
% is not.

KOGeneSets = zeros(1,length(fbamodel.optPts));
x = ~ x;

for ii = 1:length(fbamodel.optPts)
    if prod(fbamodel.optPts{ii} == ' ')~=1
        KOGeneSets(ii) = ~eval(fbamodel.optPts{ii});
    end
end

end