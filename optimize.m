function ideal = optimize(landings, conditions)
% Finds the conditions that produce the best launch
% Inputs:
%   landings: landing point matrix given by montecarlo
%   conditions: condition array given by montecarlo
% Outputs:
%   ideal: the condition set that are best (a cell array)
[~, where] = max(landings(:,2));
ideal = conditions{where};
end