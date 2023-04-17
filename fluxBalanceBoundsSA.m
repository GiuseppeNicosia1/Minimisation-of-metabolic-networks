function [vbiomass,v] = fluxBalanceBoundsSA(fbamodel)
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



%param.tmlim  = 300;
param.msglev = 1;
param.save   = 0;
%param.presol = 0;
param.itlim = -1;

nrxn   = fbamodel.nrxn;
nmetab = fbamodel.nmetab;

if(~any(fbamodel.g)) %if there is no synthetic objective, I work only on the biological one.
    yt = fbamodel.present;

    A = [ fbamodel.S; 
          sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn) ];
    b = [ zeros(nmetab, 1);
          zeros(nnz(~yt), 1) ];

    model = struct();
    model.A = A;
    model.obj = fbamodel.f;
%     if ~isfield(fbamodel,'sense')
%         model.sense = '=';
%     else
        model.sense = fbamodel.sense;
%     end
    model.rhs = b;
    model.lb = fbamodel.vmin;
    model.ub = fbamodel.vmax;
    model.vtype = 'C';
    model.modelsense = 'max';
    params = struct();
%     params.TimeLimit = 1000;
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
                error(ME.identifier,ME.message)
            end
        end
    end
    if isfield(result, 'objval')
        vbiomass = result.objval;
        v = result.x;
    else
        vbiomass = 0;
        v = zeros(nrxn,1);
    end

else % I work on both the objective
    yt = fbamodel.present;  
     
    A = [ fbamodel.S; sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn) ];
    
    b = [ zeros(nmetab, 1); zeros(nnz(~yt), 1) ];

    model = struct();
    model.A = A;
    model.obj = fbamodel.f;
    model.sense = '=';
    model.rhs = b;
    model.lb = fbamodel.vmin;
    model.ub = fbamodel.vmax;
    model.vtype = 'C';
    model.modelsense = 'max';
    params = struct();
    params.TimeLimit = 1000;
    params.OutputFlag = 0;

    result = gurobi(model, params);
    if isfield(result, 'objval')
        vbiomass = result.objval;
    else
        vbiomass = 0;
        v = zeros(nrxn,1);
        return
    end
%         v = result.x;
    %adding the constraint that require the optimality of the biological
    %fluxes
    A = [ fbamodel.S; sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn); fbamodel.f' ]; 
    b = [ zeros(nmetab, 1); zeros(nnz(~yt), 1); vbiomass ];

    model = struct();
    model.A = A;
    model.obj = fbamodel.g;
    model.sense = '=';
    model.rhs = b;
    model.lb = fbamodel.vmin;
    model.ub = fbamodel.vmax;
    model.vtype = 'C';
    model.modelsense = 'max';
    params = struct();
    params.TimeLimit = 1000;
    params.OutputFlag = 0;

    result = gurobi(model, params);
%         vbiomass = result.objval;
    if isfield(result, 'x')
        v = result.x;
    else
        vbiomass = 0;
        v = zeros(nrxn,1);
        return
    end

end
%caller = dbstack;
%disp(caller(4).name)
%fprintf('.')
end
