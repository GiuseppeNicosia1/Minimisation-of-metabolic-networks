function chromosome = spotDumKO(chromosome,GDMOobj,oldChromosome,oldChromosomeReduced)
%% delete genetic manipuation which do not actually affect the objective function.
V = GDMOobj.V;
tollerance = 0.0001;

if ~isempty(oldChromosome)
   oldChromosome = oldChromosome(:,1:V);
%    oldChromosomeReduced = oldChromosomeReduced(1,1:V); C'è un errore?
   oldChromosomeReduced = oldChromosomeReduced(:,1:V);
end

if GDMOobj.ncores > 1
    parfor ii=1:size(chromosome,1)
        x = chromosome(ii,:);
        if isempty(oldChromosome)
            boolValue = false;
        else
            [boolValue, idx] = ismember(x(1:V),oldChromosome,'rows');
        end

        if boolValue
            x(1:V) = oldChromosomeReduced(idx,:);

        else


            if(GDMOobj.flagRedirector)
                for jj=1:GDMOobj.fbamodel.nbin

                    if(x(jj)==1 && x(jj+GDMOobj.fbamodel.nbin)==1)
                        x(jj)=0;
                        x(jj+GDMOobj.fbamodel.nbin)=0;

                    end
                end
            end
            %% processing: metto i geni a 0 ad uno ad uno e controllo se cambia la funzione obiettivo o no. Se non cambia allora tengo quelli d prima
            trueFit = x((V+1):(V+GDMOobj.M));
            for jj=1:GDMOobj.V
                if(x(jj)==1)          %% se sta ad 1 controllo se effettivamente e' influente
                    z=x;
                    z(jj)=0;
                    %fx = fitnessFunction(x, GDMOobj);
                    fz = fitnessFunction(z(1:V), GDMOobj);
                    fitVariation = abs(trueFit-fz);
                    if(norm(fitVariation,1) < tollerance )
                        x(jj)=0;
                        trueFit = fz;
                        x((GDMOobj.V+1):(GDMOobj.V+GDMOobj.M)) = trueFit;
                    end
                end
            end
        end
    chromosome(ii,:) = x;
    %chromosome(ii,1:GDMOobj.V) = x;
    %chromosome(ii,(GDMOobj.V+1):(GDMOobj.V+GDMOobj.M) ) = fitnessFunction(x, GDMOobj );

    end
else
    for ii=1:size(chromosome,1)
        x = chromosome(ii,:);
        if isempty(oldChromosome)
            boolValue = false;
        else
            [boolValue, idx] = ismember(x(1:V),oldChromosome,'rows');
        end

        if boolValue
            x(1:V) = oldChromosomeReduced(idx,:);

        else


            if(GDMOobj.flagRedirector)
                for jj=1:GDMOobj.fbamodel.nbin

                    if(x(jj)==1 && x(jj+GDMOobj.fbamodel.nbin)==1)
                        x(jj)=0;
                        x(jj+GDMOobj.fbamodel.nbin)=0;

                    end
                end
            end
            %% processing: metto i geni a 0 ad uno ad uno e controllo se cambia la funzione obiettivo o no. Se non cambia allora tengo quelli d prima
            trueFit = x((V+1):(V+GDMOobj.M));
            for jj=1:GDMOobj.V
                if(x(jj)==1)          %% se sta ad 1 controllo se effettivamente e' influente
                    z=x;
                    z(jj)=0;
                    %fx = fitnessFunction(x, GDMOobj);
                    fz = fitnessFunction(z(1:V), GDMOobj);
                    fitVariation = abs(trueFit-fz);
                    if(norm(fitVariation,1) < tollerance )
                        x(jj)=0;
                        trueFit = fz;
                        x((GDMOobj.V+1):(GDMOobj.V+GDMOobj.M)) = trueFit;
                    end
                end
            end
        end
    chromosome(ii,:) = x;
    %chromosome(ii,1:GDMOobj.V) = x;
    %chromosome(ii,(GDMOobj.V+1):(GDMOobj.V+GDMOobj.M) ) = fitnessFunction(x, GDMOobj );

    end
    
end


end

