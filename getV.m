function V = getV(GDMOobj)
%getV, method for the GDMOopt class. It extrapolates the dimension of the
%domain space of the optimization problem using the other property of the
%object.
  
    V = length(GDMOobj.fbamodel.genes);

end