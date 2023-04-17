
rng('shuffle')
addpath('classes');
addpath(genpath('utils'));
addpath('/usr/local/Cluster-Apps/gurobi/8.1.1/linux64/matlab')

%%%%%%
restart = false;
pop = 100;                                  % population for the optimization
gen = 5000;                                 % maximum number of generation
ncores = 32;                                % number of concurrent threads to be used for optimization and post processing
model = 'iCEL1314_dual';                          % mat file of the FBA model to use
optMode = 'FBA';                            % one of FBA, (pFBA), MOMA
flagEssGenes = false;
flagEssCouples = false;
flagCOBRA = true;
maxred = 0.01; 
%%%%%%

fbamodel = loadFBAmodel(['model' filesep model], flagEssGenes, flagEssCouples, flagCOBRA);

if ~restart
    %Directory in which results will be stored
    results_folder = ['solutions', model, filesep, optMode, filesep];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Building the Optimizer Object%%%
    GDMOobj = GDMOopt(pop,gen,fbamodel,optMode,results_folder,ncores,maxred,0,true,10);

else
    load(['solutions' model filesep optMode filesep 'GDMOobj.mat']);
    GDMOobj.gen = gen;
    GDMOobj.ncores = ncores;
    GDMOobj.results_folder = ['solutions' model filesep optMode filesep];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Running the Optimization%%%%
GDMOobj = GDMOobj.runOptimiser();

rmpath('classes');
rmpath(genpath('utils'));
fbamodel = struct(fbamodel);
fbamodel.S = sparse(fbamodel.S);
save(['solutions' model filesep 'fbamodel'], 'fbamodel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Running Post Process Analysis%%%
addpath('post')

saveTotalChromosome(model, optMode, gen, ncores);
findMinimal(model, ncores, maxred)
matrixReac(model,ncores)
matrixReacParsimonious(model,ncores)
savevMinSol(model,ncores)
minSolsCompGenes(model)

rmpath('post')
