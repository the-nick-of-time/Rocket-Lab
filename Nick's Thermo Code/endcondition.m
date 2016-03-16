function [value,isterminal,direction]=endcondition(~,vars)
% This function defining the end condition simply follows the format laid
% out on the corresponding help page.
% This will cause ode45 to terminate when vars(3), the height, becomes zero
% while falling.
value=[vars(3)];%;vars(5)-.07];
isterminal=[1];%;0];
direction=[-1];%-1]; %only care about the falling direction
end