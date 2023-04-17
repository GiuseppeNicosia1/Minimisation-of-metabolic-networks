classdef FBAmodel
%Class describing an FBAmodel. I keep on the class only the properties
%which are strictly necessary for the optimization/analysis of the latter.
    properties
        nrxn            %number of reactions of the FBAmodel
        nmetab          %number of metabolites of the FBAmodel
        present         %vector of logical values. If the ith entry is set to 0 then it means that the ith reaction is KO
        S               %stoichiometric matrix of the FBAmodel
        vmin            %lower bounds on fluxes
        vmax            %upper bounds on fluxes
        f               %biological objective
        G               %matrix of GPR mapping
        pts             %gene sets of the model
        genes           %cell array of all the genes in the model
        ngenes          %number of genes of the model
        nbin            %number of gene sets in the model
        optPts          %An evaluable string associated to the FBAmodel genesets boolean expression, in which gene names are replaced by x(ii)
        g               %synthetic objective vector
        numGenesInvlvd  %number of genes in each gene sets
        genesBU
        ngenesBU
        essentialGenes
        essentialCouples
        flagEssGenes
        flagEssCouples
        sense
    end
    
    methods
        
        [vbiomass,v] = fluxBalance(fbamodel)
        [vbiomass,v] = pfluxBalance(fbamodel)
        [vbiomass,v] = MOMAfun(fbamodel,vWT)
        KOGeneSets = evalGeneSets(fbamodel,x)
        numGenesInvlvd = CardinalityGeneSets(fbamodel)
        fbamodel = setSynthObjective(fbamodel,idx,coeff)
        fbamodel = redirectFluxes(fbamodel , x, wildFluxes)
        
        function fbamodel = FBAmodel(S,vmin,vmax,f,genes,optPts,genesBU,essentialGenes,essentialCouples, flagEssGenes, flagEssCouples, sense)
            %Constructor for the FBAmodel class.
            %Input:
            %      -S, the stoichiometric matrix of the model
            %      -vmin, lower bound on fluxes
            %      -vmax, upper bound on fluxes
            %      -f, biological objective vector
            %      -optPts, logical genes rules of reactions
            %      -pts, cell array of gene sets of the model
            %Output:
            %       -fbamodel, an FBAmodel object.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Description: All the properties not directly provided as input
            %are extrapolated from the input features. The FBAmodel object
            %returned has a null synthetic objective vector. This has to be
            %set after the model is built.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fbamodel.S = S;
            [fbamodel.nmetab,fbamodel.nrxn] = size(fbamodel.S);
            fbamodel.present = ones(fbamodel.nrxn,1);
            fbamodel.vmin = vmin;
            fbamodel.vmax = vmax;
            fbamodel.f = f;
%             fbamodel.G = G;
%             fbamodel.pts = pts;
            fbamodel.nbin = length(fbamodel.pts);
            fbamodel.genes = genes;
            fbamodel.ngenes= length(fbamodel.genes);
            fbamodel.optPts = optPts;
            fbamodel.g = zeros(size(fbamodel.f,1),size(fbamodel.f,2));
            fbamodel.genesBU = genesBU;
            fbamodel.ngenesBU = length(fbamodel.genesBU);
            fbamodel.essentialGenes = essentialGenes;
            fbamodel.essentialCouples = essentialCouples;
            fbamodel.flagEssGenes = flagEssGenes;
            fbamodel.flagEssCouples = flagEssCouples;
            fbamodel.sense = sense;
        end
        
    end







end