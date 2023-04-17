function GDMOobj = endGenRoutine(GDMOobj)
%endGenRoutine, method of the class GDMOobj. It perform routine of end
%generation to the optimization. 1) It takes out of the population
%non-deasible strains. 2)It takes out of the population redundant strains.
%3) It cleans up each strain of the population from unrelevant genetic
%manipulation. 4) It saves the current population to file. 5) It save the
%GDMOopt object to file, which shall be used for restarting the
%optimization. 6) It computes statistics and print them to the standard
%output. 7) If redirector optimization, it updates the gamma value.
%

chromosome = GDMOobj.chromosome;
save([GDMOobj.results_folder ,'solution',int2str(GDMOobj.currGen),'.mat'],'chromosome');    
%%Faccio un backup in caso l'ottimizzazione si blocchi:
save([GDMOobj.results_folder,'GDMOobj.mat'],'GDMOobj')    
GDMOobj.best_KO = -min(chromosome(:,GDMOobj.V+GDMOobj.M));
disp(['Highest Number of KO so far: ', num2str(GDMOobj.best_KO)]);
GDMOobj.feasibleFits = chromosome(:,(GDMOobj.V+1):(GDMOobj.V+GDMOobj.M));
GDMOobj.numOfNotDominated = length(  find( chromosome(:,GDMOobj.V+GDMOobj.M+1)==1  )  );
end
