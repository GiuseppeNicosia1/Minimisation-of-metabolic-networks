 function [vbiomass,v,distance] = MOMA(fbamodel, vWT)

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
% model.obj = fbamodel.f;
model.obj = -vWT;
model.Q = 0.5*sparse(eye(nrxn));
model.sense = '=';
model.rhs = b;
model.lb = fbamodel.vmin;
model.ub = fbamodel.vmax;
model.vtype = 'C';
model.modelsense = 'min';
params = struct();
%     params.TimeLimit = 1000;
params.OutputFlag = 0;
params.FeasibilityTol = 1.0e-9;
params.OptimalityTol = 1.0e-6;keyboard
result = gurobi(model, params);
if isfield(result, 'objval')
    distance = result.objval;
    v = result.x;
    vbiomass = v(find(fbamodel.f));
else
    vbiomass = 0;
    distance = Inf;
    v = zeros(nrxn,1);
end

end
