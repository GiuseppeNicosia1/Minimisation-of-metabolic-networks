function resObj = getStrictFronts( resObj,epsilonFrontPop )
%getStrictFronts method of the class optResult. It compute strict fronts,
%from the results of the epsilon dominance. I.e. it gets pf points out of
%the first epsilon front. It gets points of the first epsilon front out of
%the second epsilon front....

prevFront = resObj.pfFit;          
for ii = 1:length(resObj.epsilon) %compute stricts epsilon-fronts (i.e. such that it is not in the previous fronts
    currFront = resObj.epsilonFrontFit{:,:,ii};
    currPop = epsilonFrontPop{:,:,ii};
    [~, diffIdxs] = setdiff(currFront,prevFront,'rows');
    prevFront = currFront;
    currFront= currFront(diffIdxs,:);
    %currPop = currPop(diffIdxs,:);
    currPop = currPop(diffIdxs);
    resObj.epsilonFrontFit{:,:,ii} = currFront;
    epsilonFrontPop{:,:,ii} = currPop;
    %save([option.results_folder,'strictPop_epsilon_',num2str(option.epsilon(ii)) ,'.txt'],'currPop','-ascii');
    save([resObj.results_folder,'strictPop_epsilon_',num2str(resObj.epsilon(ii)) ,'.mat'],'currPop');
    save([resObj.results_folder,'strictFit_epsilon_',num2str(resObj.epsilon(ii)) ,'.txt'],'currFront','-ascii');
end



end

