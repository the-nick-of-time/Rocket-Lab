function D=drag(velocity,rho_air,drag_coefficient,A)
% This function calculates the aerodynamic drag force acting on the rocket
% during flight. It uses the same equation for the entire flight time.
    D= -.5 * velocity.^2 * rho_air * drag_coefficient * A;
end