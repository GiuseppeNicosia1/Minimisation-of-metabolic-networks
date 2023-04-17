function f = fitnessFunction(x, GDMOobj )

%Computing Objective Functions for Minimal Cell
GDMOobj.fbamodel.present = true(GDMOobj.fbamodel.nrxn,1);
x = ~x;
if GDMOobj.fbamodel.flagEssGenes
    y = true(GDMOobj.fbamodel.ngenesBU,1);
    y(~GDMOobj.fbamodel.essentialGenes) = x;
    x = y;
end
for i = 1: GDMOobj.fbamodel.nrxn
    if ~isempty(GDMOobj.fbamodel.optPts{i})
        if ~eval(GDMOobj.fbamodel.optPts{i})
            GDMOobj.fbamodel.vmin(i) = 0;
            GDMOobj.fbamodel.vmax(i) = 0;
        end
%         GDMOobj.fbamodel.present(i) = eval(GDMOobj.fbamodel.optPts{i});
    end
end
%[vbiomass, ~] = fbamodel.fluxBalance();
[vbiomass, ~] = GDMOobj.invokeFBAobject();
f(1) = -vbiomass; % predicted biomass
f(2) = -sum(~x(1:GDMOobj.fbamodel.ngenesBU)); % KOs performed

%% Check for error
if length(f) ~= GDMOobj.M
    error('The number of decision variables does not match you previous input. Kindly check your objective function');
end

end

