function pEpsIdxs = epsilonRelax(resObj, pfNormFit, normFit,kk)
%epsilonRelax, method of the class optResult. Given an input pareto front,
%it "relax" it, in order to compute the various epsilon fronts required by
%the user. Return the index associated to the front.

	[N,~]=size(normFit);
	[Npf,~]=size(pfNormFit);
    pfNormFit = abs(pfNormFit);
    normFit = abs(normFit);
    %total_scaled=ScaleFit(total_population,V+1,M);
    %pf_scaled=ScaleFit(pf_fit,1,M);
    pEpsIdxs = [];

	for ii=1:N
		candidate = normFit(ii,:);
		for jj=1:Npf
            pfPoint = pfNormFit(jj,:);
			if(norm(candidate-pfPoint) < resObj.epsilon(kk))	%%se sono abbastanza vicini allora ok
                pEpsIdxs = [pEpsIdxs, ii];
                break;
            end
        end
    end
end