function [F,I]=thrust(t,vars,Vbottle,c_d,pressurefn,p_atm,A,gravity,rho_water,mdotfn,m_air,R,gamma)
% This function returns the force of thrust at any given time. 

global stage;
[p,pstar,p_t]=pressurefn(t,vars);
switch stage
    case 1
        F = 2 .* c_d .* (p-p_atm) .* A;
        I = sqrt(2.*(p-p_atm)./rho_water)./gravity;
    case 2
        rho=m_air./Vbottle;
        T=p./(rho.*R);
        T_t= 2.*T./(gamma+1);
        v_t= sqrt(gamma.*R.*T_t);
        F = (v_t.*mdotfn(t,vars) + (pstar-p_atm).*A);
        I = (F./(mdotfn(t,vars).*gravity));
    case 3
        rho=m_air./Vbottle;
        T=p./(rho.*R);
        T_t=T.*(p./p_atm).^((gamma-1)./gamma);
        M_t= sqrt((T_t./T - 1).*(2./gamma-1));
        v_t=M_t.*sqrt(gamma.*R.*T_t);
        F = mdotfn(t,vars).*v_t;
        I = v_t;
    case 4
        F=0;
        I=0;
end


end