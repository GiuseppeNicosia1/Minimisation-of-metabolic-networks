function GDMOobj = preOpt(GDMOobj)
%preOpt, method of the GDMOopt class. It performs initial checking on the
%GDMO obj, to see if the properties are consistent. Then it initialize
%fields associated to the wyld type flux distribution and to the utopian
%maximum yield of synthetic objective function.
%

disp('Optimization set for Minimal Cell GDMO'); 
if(any(GDMOobj.fbamodel.g))
    warning('...However a synthetic objective is provided... I will ignore it.'); 
end
GDMOobj.fbamodel.g = zeros( size(GDMOobj.fbamodel.g,1), size(GDMOobj.fbamodel.g,2)  );

GDMOobj = initUtopianFields(GDMOobj);

end

function GDMOobj = initUtopianFields(GDMOobj)

GDMOobj.fbamodel.g = zeros(size(GDMOobj.fbamodel.g,1),size(GDMOobj.fbamodel.g,2)) ;

[GDMOobj.utopianBio ,GDMOobj.wildFluxes ] =  GDMOobj.invokeFBAobject();

% GDMOobj = GDMOobj.computeMinBiomass();
    
end