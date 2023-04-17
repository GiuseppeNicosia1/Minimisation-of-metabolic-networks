function GDMOobj = runOptimiser(GDMOobj)
%runOptimiser, method of the class GDMOopt. It run the optimization
%specified by the input object GDMOobj.

startTime = tic;

%Print out the kind of simulation run, and checking that input is
%consistent
if GDMOobj.currGen <= 1
    GDMOobj = GDMOobj.preOpt();
end

if (GDMOobj.ncores > 1)
    if(isempty(gcp('nocreate')))
        stops = 0;
        while true
            try
                parpool('local',GDMOobj.ncores);
                break
            catch ME
                stops = stops + 1;
                if stops > 3
                    error(ME.identifier,ME.message)
                end
            end
        end
    end
end

%making output directory
if(~exist(GDMOobj.results_folder,'dir'))
    mkdir(GDMOobj.results_folder);
end

%computing an initial chromosome
if GDMOobj.currGen <= 1
    GDMOobj.chromosome = GDMOobj.initPop();
end
%evaluating objective function in initial chromosome
if GDMOobj.currGen <= 1
    GDMOobj = GDMOobj.evalFit();%chromosome(i,1:V), M, V, fbamodel,option);
end  
%% Sort the initialized population
% Sort the population using non-domination-sort. This returns two columns
% for each individual which are the rank and the crowding distance
% corresponding to their position in the front they belong. At this stage
% the rank and the crowding distance for each chromosome is added to the
% chromosome vector for easy of computation.
%GDMOobj.chromosome = GDMOobj.nonDominationSort();
if GDMOobj.currGen <= 1
    GDMOobj.chromosome = non_domination_sort_mod(GDMOobj.chromosome, GDMOobj.M, GDMOobj.V);
    GDMOobj.chromosome = [GDMOobj.chromosome zeros(size(GDMOobj.chromosome,1),1)];
end

%Initializing statistic plots
if(GDMOobj.plotFlag)
    initPlot(GDMOobj);
end
%% Start the evolution process
% The following are performed in each generation
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.
initGen = GDMOobj.currGen; 
for ii = initGen : GDMOobj.gen
    GDMOobj.currGen = ii;
    startGen = tic;
    disp(['Current Generation: ', int2str(ii)]);

    % Perform Mutation operator
    %disp('Applying genetic operators...')
    %disp('Curr size of parent chromosome:')
    %size(parent_chromosome)
%     offspring_chromosome = genetic_operator(parent_chromosome, GDMOobj);
    tic
    offspring_chromosome = genetic_operator(GDMOobj.chromosome, GDMOobj);
    disp(['Genetic Operator took: ' num2str(toc) 'seconds'])
    % Intermediate population
    % Intermediate population is the combined population of parents and
    % offsprings of the current generation. The population size is two
    % times the initial population.
    
    [main_pop,~] = size(GDMOobj.chromosome);
    [offspring_pop,~] = size(offspring_chromosome);
    % intermediate_chromosome is a concatenation of current population and
    % the offspring population.

    intermediate_chromosome(1:main_pop,:) = GDMOobj.chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : GDMOobj.M+GDMOobj.V) = ...
        offspring_chromosome;

    % Non-domination-sort of intermediate population
    % The intermediate population is sorted again based on non-domination sort
    % before the replacement operator is performed on the intermediate
    % population.
    %disp('Non domination sorting...')
    tic
    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, GDMOobj.M, GDMOobj.V);
    disp(['Non-dom-sort took: ' num2str(toc) ' seconds'])
    % Perform Selection
    % Once the intermediate population is sorted only the best solution is
    % selected based on it rank and crowding distance. Each front is filled in
    % ascending order until the addition of population size is reached. The
    % last front is included in the population based on the individuals with
    % least crowding distance
    %disp('Raplicing chromosome...')
    GDMOobj.lastChromosome = GDMOobj.chromosome;
    GDMOobj.chromosome = replace_chromosome(intermediate_chromosome, GDMOobj.M, GDMOobj.V, GDMOobj.pop);
    % for each iteration you can eliminate the ridontand gene knockout or at the end of the optimization   
    
    % NOW A PROCEDURE FOR HEATING THE OBSOLETE
    % SOLUTIONS AND UPDATE THE AGING
    
    tic
    GDMOobj.chromosome = heating(GDMOobj);
    disp(['Heating took: ' num2str(toc) ' seconds'])
    
    %Current generation cycle performed. Perforfming some final routine for
    %statistical purposes
    %disp('end generation routine...')
    GDMOobj = GDMOobj.endGenRoutine();%(chromosome,option,ii); 
    calc = toc(startGen) ;
    calc = compTime(calc,GDMOobj.gen,ii);   
    if(GDMOobj.ncores > 1)
        if(~mod(ii,1000))%Every 1000 generation I restart the parallel pool to avoid memory waste...
            delete(gcp);
            stops = 0;
            while true
                try
                    parpool('local',GDMOobj.ncores);
                    break
                catch ME
                    stops = stops + 1;
                    if stops > 3
                        error(ME.identifier,ME.message)
                    end
                end
            end
        end
    end
    if ~mod(ii,100)
    %if(1)
        %clc
        disp([int2str(ii),' generations completed'])
        if(GDMOobj.plotFlag)
            EstimatedTimesH = calc/3600; %estimated remaining time in hours    	
            endGenPlots(EstimatedTimesH, GDMOobj.best_KO, GDMOobj.numOfNotDominated , GDMOobj.feasibleFits,ii, GDMOobj.results_folder );
        end
    end
end
%closing down the parallel pool
if(GDMOobj.ncores > 1)
    delete(gcp);
end

close all
save([GDMOobj.results_folder, 'GDMOobj.mat'],'GDMOobj');

end



function calc = compTime(calc,gen,ii)
    remain=gen-ii;
    calc=calc*remain;
    disp('remaining time, approximately:');
    disp(datestr(calc/86400, 'dd:HH:MM:SS.FFF'));
end

function endGenPlots(EstimatedTimes,best_synt,numOfNotDominated , feasibleFits,ii,results_folder)
    subplot(2,2,1);
    hold on
    %plot(EstimatedTimes,'r');
    plot(ii,EstimatedTimes,'or');
    if(~isempty(best_synt))
        subplot(2,2,2);
        hold on
        plot(ii,best_synt,'xr');
    end
    if(~isempty(numOfNotDominated))
        subplot(2,2,3);
        hold on
        plot(ii,numOfNotDominated,'xb');
    end
    if(~isempty(feasibleFits))
        subplot(2,2,4);
        hold on
        plot(abs(feasibleFits(:,1)), abs(feasibleFits(:,2)) ,'*g');
        title(['Feasible Points at generation ', int2str(ii)]);
    end
    drawnow
    savefig([results_folder, 'statistics.fig']);
end

function initPlot(GDMOobj)
figure
subplot(2,2,1);
xlabel('generation');
ylabel('Aprroxymate ramaining time [h]');
title('Generation vs. Remaining time')
grid on

subplot(2,2,2);
xlabel('generation');
ylabel('Best KO');
title('Generation vs. Greatest Number of KO')
grid on

subplot(2,2,3);
xlabel('generation');
ylabel('# of Non Dominated Points');
title('# of (Distinct) Non Dominated Points in current population')
grid on

subplot(2,2,4);
xlabel('Biomass [h^{-1}]');
ylabel('KO');
title('Feasible Points in current population (only at first Corner...)')
grid on

end
