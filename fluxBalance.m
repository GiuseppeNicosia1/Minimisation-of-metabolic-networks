function [vbiomass,v] = fluxBalance(fbamodel)
% fluxBalance method for the flux balance analysis of the FBAmodel object
% fbamodel.
% Output: 
%         -vbiomass, i.e. optimal biomass value
%         -v, optimal flux vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description: the code just perform the FBA analysis as usual. KO reactions
%are taken off the Stoichiometric matrix and linear programms are solved.
%Two cases are considered: 1)The synthetic objective vector is identically
%null. The code thus assumes that there is no interest in synthetic obj and
%hence only perform the linear program optimization on the biological
%objective. If whereas the synthetic objective has non-null entries, the
%code solves solve two linear programs. The first one on the biological
%objective. The second one on the synthetic objective, contrained to the
%optimality of the biological fluxes, as dictated by the firt optimization
%performed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrxn   = fbamodel.nrxn;
nmetab = fbamodel.nmetab;

yt = fbamodel.present;

A = [ fbamodel.S; 
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn) ];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1) ];

model = struct();
model.A = A;
model.obj = fbamodel.f;
% model.sense = '=';
model.sense = fbamodel.sense;
model.rhs = b;
model.lb = fbamodel.vmin;
model.ub = fbamodel.vmax;
model.vtype = 'C';
model.modelsense = 'max';
params = struct();
params.TimeLimit = 1000;
params.OutputFlag = 0;
params.FeasibilityTol = 1.0e-9;
params.OptimalityTol = 1.0e-6;

stops = 0;
while true
    try
        result = gurobi(model, params);
        break
    catch ME
        stops = stops + 1;
        if stops > 1000
            keyboard
            error(ME.identifier,ME.message)
        end
    end
end

% result = gurobi(model, params);
if isfield(result, 'objval')
    vbiomass = result.objval;
    v = result.x;
else
    vbiomass = 0;
    v = zeros(nrxn,1);
end


end
