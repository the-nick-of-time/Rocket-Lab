function [value,isterminal,direction]=endstage3(t,vars,pressurefn,p_atm)
[p,~]=pressurefn(t,vars);
value=p-p_atm;
isterminal=1;
direction=0;
end