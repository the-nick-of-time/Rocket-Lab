function integrate(rocket)

switch floor(rocket.type())
    case 1
        % Isp model
        opts = odeset('Events', @rocket.endcondition);
    case 2
        % Thrust model
        opts = odeset('Events', @rocket.endcondition);
    case 3
        % Thermodynamic model
        opts = odeset('NonNegative',[0 0 0 0 0 0 1 1 1],...
            'Events', @rocket.endcondition);
end
trange = [0 20];
ode45(@rocket.derivatives, trange, rocket.initialconditions(), opts);
end