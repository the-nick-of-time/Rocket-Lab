function dmdt=mdot(t,vars,pressurefn,m_air, c_d,rho_water,p_atm,A,Vbottle,gamma,R)
% This function gives the rate of change of the rocket's mass over time.
% First section: depends on pressure differential?
% Second section: depends on whether the flow is choked or not. The flow is
%   choked if $p > p_*$ where $$ p_* = p \left ( \frac{2}{\gamma +1} \right
%   )^\frac{\gamma}{\gamma -1}$$
% Third section: is zero
global stage;

[p,pstar]=pressurefn(t,vars);

switch stage
    case 1
        dmdt = -c_d*rho_water*A*sqrt(2*(p-p_atm)/rho_water);
    case 2
        rho=m_air/Vbottle;
        T=p/(rho*R);
        T_t= 2*T/(gamma+1);
        rho_t= pstar/(R*T_t);
        v_t= sqrt(gamma*R*T_t);
        dmdt = -rho_t*v_t*A*c_d;
    case 3
        rho=m_air/Vbottle;
        T=p/(rho*R);
        T_t=T*(p/p_atm)^((gamma-1)/gamma);
        rho_t=p_atm/(R*T_t);
        M_t= sqrt( (T_t/T - 1) * (2/(gamma-1)) );
        v_t=M_t*sqrt(gamma*R*T_t);
        dmdt = -rho_t*v_t*A*c_d;
    case 4
        dmdt=0;
end

end