function landingpoint = integrate(rocket)

switch floor(rocket.type())
    case 1
        % Isp model
        opts = odeset('Events', @rocket.endcondition);
    case 2
        % Thrust model
        opts = odeset('Events', @rocket.endcondition,...
            'MaxStep', .0025);
    case 3
        % Thermodynamic model
        opts = odeset('NonNegative', [7 8 9],...
            'Events', @rocket.endcondition);
end
trange = [0 20];
[t, vars] = ode45(@rocket.derivatives, trange, ...
    rocket.initialconditions(), opts);
rocket.finalize(t, vars);

landingpoint = [vars(end,4) vars(end,5) vars(end,6)];
end