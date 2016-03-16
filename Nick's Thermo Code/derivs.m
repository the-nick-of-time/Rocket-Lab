function d=derivs(t,vars,vdot,mdot,thetadot,Vdot)
% returns relevant derivatives
% elements of the vars vector: (d is the elementwise derivatives)
% 
% V=      velocity magnitude (speed)
% x=      horizontal distance
% y=      vertical distance
% theta=  angle from horizontal
% m=      mass of the whole rocket
% vol=    volume of air inside the rocket

global stage setstage;
setstage(t,vars)

d=[ vdot(t,vars);
    vars(1)*cos(vars(4));
    vars(1)*sin(vars(4));
    thetadot(t,vars);
    mdot(t,vars);
    Vdot(t,vars)
  ];
end