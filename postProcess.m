function postProcess(resObj,GDMOobj)
%postProcess, method of the class optResult. It applies the post processing routines to the results of GDMO
%optimization.
%Input: 
%      -resObj, i.e. an optResult object
%optimization
%computational time of the postProcess routine.
%Usage: 
%      postProcess(resObj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description: It applies the post processing routine to the results of the
%GDMO optimization. First off it loads all the feasible points explored
%during the optimization and stores them in two ways: fit_tot, i.e. the
%codomain points explored, and total_population, i.e. the domain points
%explored. As a second step we retrieve the unique points of fit_tot (up to
%a three digit approximation), and save the corrisponing domain points. The
%function hence compute the pareto front, and if specified in the input
%structure, the espilon-fronts. All the results are saved in .txt and
%matlab file. Finally the results are plotted using the plotting_wrap
%auxiliary function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parpool('local',2); %AGGIUNTO
initGen = 1; %AGGIUNTO, TOGLIERE E SOSTITUIRE CON 1 NEL FOR O ASSEGNARE A 1

if ~resObj.flagMinCell
    ss=[resObj.results_folder 'solution'];
    oldChromosome = [];
    oldChromosomeReduced = [];
    disp('Removing redundant manipulation from the whole explored feasible region.')
    disp('This can be very expensive, consider skipping this part if it takes too much...')
    for ii = initGen:resObj.gen
        disp(['Gen: ', int2str(ii)]);
        tic

        s = [ss,int2str(ii),'.mat'];   
        load(s);
        chromosomeIn = chromosome;
        chromosome = spotDumKO(chromosomeIn,GDMOobj,oldChromosome,oldChromosomeReduced);
        oldChromosome = chromosomeIn;
        oldChromosomeReduced = chromosome;
        save(s,'chromosome');
        toc;
        if ~mod(ii,100) % AGGIUNTO
            delete(gcp);
            parpool('local',2);
        end
    end
end




tic;
ss=[resObj.results_folder 'solution'];
total_population = {};
disp('Loading all the feasible points found during the optimization...')
for ii=1:(resObj.gen)
	s = [ss,int2str(ii),'.mat'];      
	load(s);
    %if 1
    %    chromosome = spotDumKO(chromosome,GDMOobj);
    %end
	%%loading population
    for jj=1:size(chromosome,1)

        total_population{end+1} = encodeStrain(chromosome( jj,1:resObj.V) , resObj   ) ;
    end
    resObj.fit = [resObj.fit ; chromosome(:,(resObj.V+1):(resObj.V+resObj.M))];        
end
    
clear chromosome
clear GDMOobj
if(isempty(resObj.fit))
    disp('No feasible points have been found during the optimization... Consider running the optimization again increasing the population size...')
    return
end
    
    
%salvo i punti esplorati al completo
fit = resObj.fit;
save([resObj.results_folder, 'FULLfit.mat'],'fit');
save([resObj.results_folder, 'FULLtotal_population'],'total_population');
clear fit;
%%%Qui dovrei eseguire il seguente statement:
%total_population=unique(total_population, 'rows');	%%elimino le ripetizioni 
%Ma essendo troppo costoso faccio il lavoro su fit_tot:
disp('Deleting repition from population. This may take a while....')
resObj.fit = 10^-3*round(10^3*resObj.fit); %rounding up small variations on the objective functions %SHOULD BE SCALE FREE!
[resObj.fit, uniqueIdxs] = unique(resObj.fit, 'rows');
disp('...Done.')
%Così facendo tengo UNA sola copia tra quelle che mi permettono di
%ottenere tale punto nel codominio
total_population = total_population(uniqueIdxs);
clear uniqueIdxs
absFitTot = abs(resObj.fit);
%absPopTot = abs(total_population);
save([resObj.results_folder 'absFitTot.mat'],'absFitTot');
%save([option.results_folder 'absPopTot.mat'],'absPopTot');
save([resObj.results_folder 'total_population.mat'],'total_population');
clear absFitTot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARETO AND EPSILON NON-DOM ANALYSIS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
disp('Computing Pareto Front...')
pfIdxs = resObj.computeParetoFront();
%absPfPop = absPopTot(pfIdxs,:);
totPfPop = total_population(pfIdxs);
resObj.pfFit = abs(resObj.fit(pfIdxs,:));
absPfFit = resObj.pfFit;
save([resObj.results_folder,'paretoFit.txt'],'absPfFit','-ascii');
%save([option.results_folder,'absPfPop.txt'],'absPfPop','-ascii');
save([resObj.results_folder,'totPfPop.mat'],'totPfPop');
disp('...Done')
%clear totPfPop
clear absPfFit

    
%normalizzo
[a,b]= resObj.getAsBs();
normFit = ScaleFit(resObj.fit,a,b);
pfNormFit = normFit(pfIdxs,:);
clear pfIdxs
resObj.epsilon = sort(abs(resObj.epsilon));
    
%Parallel computation of Epsilon Fronts...
disp('Computing epsilon-fronts...');
epsilonFrontPop = cell(1,1,1);
%if(option.parallel) %TODO!
%    parfor ii=1:length(option.epsilon)
%        pEpsIdxs = epsilon_relaxation(pfNormFit, normFit, option.epsilon(ii));
        %tempPop = absPopTot(pEpsIdxs, :);
%        tempFit = absFitTot(pEpsIdxs, :);
%        epsilonFrontPop{:,:,ii} = tempPop;
%        epsilonFrontFit{:,:,ii} = tempFit;        
%    end
%    for ii=1:length(option.epsilon)
%        tempPop = epsilonFrontPop{:,:,ii};
%        tempFit = epsilonFrontFit{:,:,ii};
%        save([option.results_folder,'pop_epsilon_',num2str(option.epsilon(ii)) ,'.txt'],'tempPop','-ascii');
%        save([option.results_folder,'fit_epsilon_',num2str(option.epsilon(ii)) ,'.txt'],'tempFit','-ascii');     
%    end
%else
for ii=1:length(resObj.epsilon)
    pEpsIdxs = resObj.epsilonRelax(pfNormFit, normFit,ii); 
   %tempPop = absPopTot(pEpsIdxs, :);
    tempPop = total_population(pEpsIdxs);
    tempFit = abs(resObj.fit(pEpsIdxs, :));
    epsilonFrontPop{:,:,ii} = tempPop;
    resObj.epsilonFrontFit{:,:,ii} = tempFit;        
   %save([option.results_folder,'pop_epsilon_',num2str(option.epsilon(ii)) ,'.txt'],'tempPop','-ascii');
    save([resObj.results_folder,'pop_epsilon_',num2str(resObj.epsilon(ii)) ,'.mat'],'tempPop');
    save([resObj.results_folder,'fit_epsilon_',num2str(resObj.epsilon(ii)) ,'.txt'],'tempFit','-ascii');     
end
        
%end
disp('...Done')
clear tempPop
clear tempFit
    
resObj = resObj.getStrictFronts(epsilonFrontPop);
%[bestBioFit, bestSynthFit, Close2UtopianFit] = computeNotablePoints(pfNormFit,absPfPop,absPfFit,option);
if(resObj.flagCorner)
    [resObj.bestBioFit, resObj.bestSynthFit, resObj.Close2UtopianFit, resObj.bestTradeFit_allCorners] = resObj.computeNotablePoints(pfNormFit,totPfPop); 
else
    [resObj.bestBioFit, resObj.bestSynthFit, resObj.Close2UtopianFit, ~] = resObj.computeNotablePoints(pfNormFit,totPfPop); 
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOTTING%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(resObj.plotFlag)
    resObj.plotting_wrap();
end

resObj.postProcTime = toc;
save([resObj.results_folder, 'resObj.mat'],'resObj');
disp(['Post processing took ',num2str(resObj.postProcTime)]);
end






