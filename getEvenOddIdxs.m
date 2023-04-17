function [evenIdxs , oddIdxs] = getEvenOddIdxs(M)

idxs = 1:M;
evenIdxs = mod(idxs,2) == 0;
oddIdxs = mod(idxs,2) == 1;

end