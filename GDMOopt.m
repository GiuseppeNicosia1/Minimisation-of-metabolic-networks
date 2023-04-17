classdef GDMOopt
% Class describing the GDMO optimizer
    properties
        pop                 %population dimension
        gen                 %total number of generation
        currGen             %current generation
        M                   %number of objective functions to be considered
        V                   %dimension of the domain space
        results_folder      %folder in which to save optimization results
        fbamodel            %FBAmodel object to be optimized
        optMode             %how this optimization shall be performed      
        ncores              %number of cores for the optimization (default = 1)
        best_KO
        %Optimization flags, initilized by parsing the optMode property
        flagFBA
        flagpFBA
        flagMOMA
        chromosome          %current population
        initKOprob          %probability of KO in initial population (default = 0)
        plotFlag            %if true produces plots (default = true)
        feasibleFits        %feasible fit in curr pop
        numOfNotDominated   %num of non dom sol in curr pop
        maxred              %
        value               %cloning value (default = 10)
        utopianBio          %wild type biomass (it's the best possible)
        lastChromosome      %Last generation Chromosome
        wildFluxes          %flux vector for wild type
    end
    
    methods
       
        GDMOobj = runOptimiser(GDMOobj)
        GDMOobj = preOpt(GDMOobj)
        V = getV(GDMOobj)
        chromosome = initPop(GDMOobj)
        GDMOobj = evalFit(GDMOobj)
        GDMOobj = endGenRoutine(GDMOobj)
        GDMOobj = parseOptMode(GDMOobj)
        GDMOobj = computeMinBiomass(GDMOobj)
        
        %Constructor for the FBAmodel GDMOopt class
        function GDMOobj = GDMOopt(pop,gen,fbamodel,optMode,resdir,ncores,maxred,initKOprob,plotFlag,value)
            if nargin < 10
                value = 10;
                if nargin < 9
                    plotFlag = true;
                    if nargin < 8
                        initKOprob = 0;
                        if nargin < 7
                            maxred = 0.01;
                            if nargin < 6
                                ncores = 1;
                            end
                        end
                    end
                end
            end
            GDMOobj.fbamodel = fbamodel;
            GDMOobj.optMode = optMode;
            GDMOobj.gen = gen;
            GDMOobj.pop = pop;
            GDMOobj.ncores = ncores;
            GDMOobj.M = 2;
            GDMOobj.currGen = 1;
            GDMOobj.chromosome = [];
            GDMOobj.results_folder = resdir;
            GDMOobj.initKOprob = initKOprob;
            GDMOobj.plotFlag = plotFlag;
            GDMOobj.feasibleFits = [];       
            GDMOobj.numOfNotDominated=0;
            GDMOobj.best_KO=0;
            GDMOobj.value = value;           
            GDMOobj = GDMOobj.parseOptMode();
            GDMOobj.V = GDMOobj.getV();
            GDMOobj.utopianBio = 0;         %wild type biomass (it s the best possible)
            GDMOobj.maxred = maxred;
            GDMOobj.wildFluxes = [];
        end
        
    end







end