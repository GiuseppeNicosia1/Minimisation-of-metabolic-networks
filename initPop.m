function chromosome = initPop(GDMOobj)
%initPop, method of the class GDMOopt. It initilize the population of the
%first generation of the population.
    chromosome = (rand(GDMOobj.pop, GDMOobj.V) > (1-GDMOobj.initKOprob) );

end