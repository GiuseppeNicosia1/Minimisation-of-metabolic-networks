function v = pFluxBalance(fbamodel, vbiomass)

revRxns = find(fbamodel.vmin<0);
index = zeros(length(revRxns),1);

for i = 1:length(revRxns)
    fbamodel.S = [fbamodel.S -fbamodel.S(:,revRxns(i))];
    fbamodel.vmin = [fbamodel.vmin; 0];
    fbamodel.vmax = [fbamodel.vmax; -fbamodel.vmin(revRxns(i))];
    fbamodel.vmin(revRxns(i)) = 0;
    fbamodel.f = [fbamodel.f; 0];
    index(i) = length(fbamodel.vmin);
end

% yt = fbamodel.present;

A = [ fbamodel.S; 
%       sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), fbamodel.nrxn); 
      fbamodel.f' ];
b = [ zeros(fbamodel.nmetab, 1);
%       zeros(nnz(~yt), 1); 
      vbiomass ];

model = struct();
model.A = A;
% model.obj = fbamodel.f;
model.obj = ones(length(fbamodel.vmin),1);
% model.sense = '=';
model.sense = fbamodel.sense;
model.rhs = b;
model.lb = fbamodel.vmin;
model.ub = fbamodel.vmax;
model.vtype = 'C';
%%%%%%%%%%%%%%%
% model.modelsense = 'max';
model.modelsense = 'min';

params = struct();
%     params.TimeLimit = 1000;
params.OutputFlag = 0;
params.FeasibilityTol = 1.0e-9;
params.OptimalityTol = 1.0e-9;
result = gurobi(model, params);
if isfield(result, 'objval')
%     vbiomass = result.objval;
    v = result.x;
else
%     vbiomass = 0;
    v = zeros(fbamodel.nrxn+length(revRxns),1);
end

for i = 1:length(revRxns)
    v(revRxns(i)) = v(revRxns(i)) - v(index(i));
end
v = v(1:fbamodel.nrxn);
end
