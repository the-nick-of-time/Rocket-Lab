function [value,isterminal,direction]=endstage2(t,vars,pressurefn,p_atm)
[~,pstar]=pressurefn(t,vars);
value=pstar-p_atm;
isterminal=1;
direction=0;
end