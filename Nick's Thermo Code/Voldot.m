function dVdt=Voldot(p_i,p_atm,Vi,V,Vbottle,gamma,c_d,A,rho_water)
% This function returns the time-rate of change of the volume of air within
% the bottle. As water is being expelled, it makes more room for the air
% and the air expands to fill that space. Hence, during stage 1, the volume
% is increasing.
% After the air expands to fill the whole bottle, there is obviously no
% way for the air to take up more of the bottle. Therefore, $\frac{dV}{dt}$
% drops to zero.
if V < Vbottle
    dVdt = c_d * A * sqrt((2/rho_water) * ((p_i * (Vi/V)^gamma) - p_atm));
else
    dVdt=0;
end
end