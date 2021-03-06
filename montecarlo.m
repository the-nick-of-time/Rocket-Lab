function [landings, conditions] = montecarlo(type, many, static, varargin)
% Performs Monte Carlo simulation of rocket paths
% Inputs:
%   type: 1, 2, or 3, the type of rocket to use
%   many: how many simulations to perform
%   static: a statictest object for use with the Isp and thrust models
%   varargin: A number of variable name : mean value : standard variance
%       triplets, enough to define the rocket
% Outputs:
%   landings: matrix of landing points
%   conditions: the condition set associated with each landing point
if mod(length(varargin), 3)
    error('You must define var:mean:sigma triplets')
end
optsthermo = { 'position';
               'theta';
               'm_w';
               'V_a';
               'p_bottle';
               'm_bottle';
               'c_d';
               'A';
               'mu';
               'T_atm';
               'P_atm';
               'discharge';
               'wind'
             };
optsthrust = { 'null';
               'position';
               'theta';
               'm_bottle';
               'c_d';
               'A';
               'mu';
               'T_atm';
               'P_atm';
               'm_w';
               'p_bottle';
               'wind'
             };
optsisp = { 'null';
            'position';
            'theta';
            'm_bottle';
            'c_d';
            'A';
            'mu';
            'T_atm';
            'P_atm';
            'wind'
          };    

fig = figure('Position', [100 100 900 700]);
      
landings = zeros(many, 2);
conditions = cell(many, 1);
      
ispinit = {static};
thrustinit = {static};
thermoinit = {};

for j = 1:many
    for i = 1:3:length(varargin)
        opt = char(varargin{i});
        val = double(varargin{i+1});
        var = double(varargin{i+2});
        switch type
            case 1
                %Isp
                [loc, ~] = find(strcmpi(opt, optsisp));
                ispinit{loc} = val + randn(size(var)) .* var;
            case 2
                %thrust
                [loc, ~] = find(strcmpi(opt, optsthrust));
                thrustinit{loc} = val + randn(size(var)) .* var;
            case 3
                %thermo
                [loc, ~] = find(strcmpi(opt, optsthermo));
                thermoinit{loc} = val + randn(size(var)) .* var;
        end
    end
    switch type
        case 1
            %Isp
            rocket = IspModel(ispinit{:});
            TITLE = 'I_{sp} Model';
            conditions{j} = ispinit;
        case 2
            %thrust
            rocket = ThrustModel(thrustinit{:});
            TITLE = 'Thrust Model';
            conditions{j} = thrustinit;
        case 3
            %thermo
            rocket = ThermoModel(thermoinit{:});
            TITLE = 'Thermodynamic Model';
            conditions{j} = thermoinit;
    end
    try
        landing = integrate(rocket);
    catch
        warning('something bad happened, but whatever')
    end
    
    landings(j,:) = landing(1:2);
    
    rocket.makeplot3d('x', 'y', 'z', {}, {}, fig);
end

figure('Position', [100 100 900 700])
hold on

X = landings(:,1);
Y = landings(:,2);
scatter(X, Y)

m = [mean(X) mean(Y)];
P = cov(X, Y);
[eigvec, eigval] = eig(P);
colors = 'bgr';
pts = 0:pi/100:2*pi;
for i = 1:3
    xy = [cos(pts') sin(pts')] * sqrt(eigval) * eigvec';
    plot(i*xy(:,1) + m(1), i*xy(:,2) + m(2), colors(i))
end
daspect([1 1 1])
title(['Landing Points of the ' TITLE])
xlabel('Cross Range Distance')
ylabel('Range Distance')


print(gcf, '-dpng', ['./figures/Landings ' TITLE '.png'])
print(fig, '-dpng', ['./figures/' TITLE '.png'])
end