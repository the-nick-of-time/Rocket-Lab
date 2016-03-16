function [p,pstar,p_t]=pressure(t,vars,p_i,p_end,p_atm,V_i,m_air_i,m_bottle,gamma)
% This function calculates the pressure inside the bottle at a given set of
% conditions.
global stage;

switch stage
    case 1
        p = p_i * (V_i/vars(6))^gamma;
        p_t = p_atm;
        pstar=Inf;
    case 2
        m_air=(vars(5)-m_bottle);
        if m_air < 0
            m_air=m_air_i;
        end
        p = (p_end * (m_air/m_air_i)^gamma);
        pstar = p * (2/(gamma+1))^(gamma/(gamma-1));
        p_t=pstar;
    case 3
        m_air=(vars(5)-m_bottle);
        if m_air < 0
            m_air=m_air_i;
        end
        p = (p_end * (m_air/m_air_i)^gamma);
        pstar = p * (2/(gamma+1))^(gamma/(gamma-1));
        p_t=p_atm;
    case 4
        p=p_atm;
        pstar=Inf;
        p_t=p_atm;
end

end