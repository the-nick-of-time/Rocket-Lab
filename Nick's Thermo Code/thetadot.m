function dadt=thetadot(theta,velocity,gravity)
% This function gives the rate of change in velocity with respect to time.
% The only reason this is difficult is if V=0, the function is undefined.
% This can be rectified by appealing to the physical analog to this
% problem, as near the start of launch the rocket will be supported by a
% stand and therefore cannot drop. Therefore $\dot{\theta}=0$ for low V.
if velocity < 1
    dadt=0;
else
    dadt=-gravity*cos(theta)/velocity;
end
end