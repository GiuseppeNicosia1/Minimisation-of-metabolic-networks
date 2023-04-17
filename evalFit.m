function GDMOobj = evalFit(GDMOobj)
%evalFit, method of the class GDMOopt. It evaluates the fitness function of
%the current population (i.e. GDMOobj.chromosome)
%tic;
for ii = 1:size(GDMOobj.chromosome,1)
    x = GDMOobj.chromosome(ii,1:GDMOobj.V);
    GDMOobj.chromosome(ii,(GDMOobj.V+1):(GDMOobj.V+GDMOobj.M)) = fitnessFunction(x,GDMOobj);
    %if ~mod(ii,100)
        %toc
        %tic;
    %end
end
%toc
end
