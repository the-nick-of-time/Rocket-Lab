function [Cd, stdCd] = calculateCd(data, area, mode)
% Calculates the drag coefficient of the rocket.
% Inputs:
%   data: a windtunneldata object holding the data from our 
%   area: the cross-sectional area of the rocket
% Outputs:
%   Cd: the drag coefficient, as determined by the formula F = q C_d A
%   stdCd: the standard deviation of this C_d

rho = 1.04187;
airspeed = data.getAirspeed();
force = data.getSting('axial');
Cd_ = [];
stdCd_ = [];
for s = [5 10 15 20 25]
    i = data.findByAirspeed(s, 1);
    if numel(i)
        temp = 2 * force(i) ./ (airspeed(i) .^2 * rho * area);
        stdCd_ = [stdCd_;std(temp)];
        Cd_ = [Cd_;mean(temp)];
    end
end

[~, place] = min(stdCd_);

% Cd_ = 2 * force ./ (airspeed .^ 2 * rho * area);
if nargin == 2 || strcmpi(mode, 'average')
%     Cd = mean(Cd_);
%     stdCd = std(Cd_);
    Cd = Cd_(place);
    stdCd = stdCd_(place);
else
    Cd = Cd_;
end

end