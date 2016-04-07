function landingpoint = integrate(rocket)
% Integrates the model
% Inputs:
%   rocket: model object
% Outputs:
%   landingpoint: [x y z] coordinates of impact (z should always be 0)
switch floor(rocket.type())
    case 1
        % Isp model
        opts = odeset('Events', @rocket.endcondition,...
            'MaxStep', .1);
    case 2
        % Thrust model
        opts = odeset('Events', @rocket.endcondition,...
            'MaxStep', .0025,...
            'NonNegative', [7]);
    case 3
        % Thermodynamic model
        opts = odeset('NonNegative', [7 8 9],...
            'Events', @rocket.endcondition,...
            'MaxStep', .1);
end
trange = [0 20];
[t, vars] = ode45(@rocket.derivatives, trange, ...
    rocket.initialconditions(), opts);
rocket.finalize(t, vars);

landingpoint = [vars(end,4) vars(end,5) vars(end,6)];
end