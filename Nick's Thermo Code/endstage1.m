function [value,isterminal,direction]=endstage1(~,vars,Vbottle)
value=Vbottle-vars(6);
isterminal=1;
direction=0;
end