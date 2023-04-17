function [a,b]= getAsBs(resObj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a = zeros(1,resObj.M);
b = zeros(1,resObj.M);
[evenIdx , oddIdx] = resObj.getFFidxs();
kk = 1;
for ii = 1:resObj.M
    if(oddIdx(ii))
        a(ii) = -resObj.utopianBio(kk);
        b(ii) = 0; %TODO MIN BIOMASS
    elseif(evenIdx(ii))
        a(ii) = -resObj.utopianSynth(kk);
        b(ii) = 0;
        kk = kk + 1;
    end

end
end

