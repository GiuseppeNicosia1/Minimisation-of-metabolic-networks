classdef optResult
    %Class for handling optimization results
    properties
        results_folder              %folder in which solutions are stored
        gen                         %number of generation to analyse
        fit                         %fit explored
        pfFit                       %observed pf
        epsilon                     %vector of epsilon to compute
        epsilonFrontFit             %observed epsilon front fits
        optMode                     %optimization mode
%        
        flagMinCell                 %=1 if min cell optimization
        flagGene                    %=1 if genes are to be considered (and not gene sets)
        flagRedirector              %=1 if redirector optimization
        flagKO                      %=1 if KO optimization
        flagCorner                  %=1 if Corner analysis is being performed
        flagWC              %=1 if Whole cell for minimal cell of M. genitalium
%        
        M                           %codomain dim
        V                           %domain dim
        plotFlag                    %to plot or not to plot.
        bestBioFit                  %best biological fit
        bestSynthFit                %best synthetic fit
        Close2UtopianFit            %closest to utopian strain
        postProcTime                %computational time of Post Processing
        bestTradeFit_allCorners     %for Corner Analysis. It is the best trade off among all Corners
        utopianBio                  %wild type biomass (it's the best possible)
        utopianSynth                %best possible synthetic objective
    end
    
    methods
        
        resObj = parseOptMode(resObj)
        [bestBioFit, bestKoFit, Close2UtopianFit,bestTradeFit_allCorners] = computeNotablePoints(resObj,pfNormFit,absPfPop)
        pfIdxs = computeParetoFront(resObj)
        pEpsIdxs = epsilonRelax(resObj, pfNormFit, normFit,kk)
        resObj = getStrictFronts( resObj,epsilonFrontPop )
        plotting_wrap(resObj)
        postProcess(resObj,GDMOobj)
        [a,b]= getAsBs(resObj)
        
        function resObj = optResult(results_folder,gen,optMode,epsilon,V,M,plotFlag,flagCorner,utopianBio,utopianSynth)
            %Constructor for the resObj class.
            resObj.results_folder = results_folder;
            resObj.gen = gen;
            resObj.epsilon = epsilon;
            resObj.optMode = optMode;
            resObj.fit = [];
            resObj.pfFit = [];
            resObj.epsilonFrontFit = cell(1,1,1);
            resObj.V = V;
            resObj.M = M;
            resObj.plotFlag = plotFlag;
            resObj.bestBioFit = [];
            resObj.bestSynthFit = [];     
            resObj.Close2UtopianFit = [];
            resObj.postProcTime = 0;
            resObj = resObj.parseOptMode();
            resObj.flagCorner = flagCorner;
            resObj.bestTradeFit_allCorners = [];
            resObj.utopianBio = utopianBio;
            resObj.utopianSynth = utopianSynth;
        end
    end


end
