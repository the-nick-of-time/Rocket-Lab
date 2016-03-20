function integrate(init_cond, derivsfn, rocket)

switch length(init_cond)
    case 6
        % Isp model
        
    case 7
        % Thrust model
        
    case 9
        % Thermodynamic model
        odeset('NonNegative',[0 0 0 0 0 0 1 1 1],...
            'Events', @object.endcondition);
end
end