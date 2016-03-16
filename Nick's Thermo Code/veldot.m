function dvdt=veldot(t,vars,thrust,drag,gravTangent)
% gravTangent is the component of the gravitational acceleration acting along
% the velocity of the rocket.
% gravity is the acceleration of gravity, around 9.81 m/s^2.
% Each of the input functions are parametrized to take as arguments only
% the integration variables, $t$ and $vars$.
    dvdt = (thrust(t,vars) + drag(t,vars)) /vars(5) + ...
           gravTangent(t,vars);
end