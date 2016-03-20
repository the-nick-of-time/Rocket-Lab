classdef ThermoModel < BottleRocket
    properties(Access = {?BottleRocket})
        m_air
        m_air_i
        m_water
        V_air
        V_air_i
        P_bottle
        P_bottle_i
        Pstar
        A_throat% = a number
        discharge% = a number
        stage = 1;
        type = 3 + rand(1);
    end
    methods
        function self = ThermoModel(vi, ri, theta, m_air, m_water, V_air,...
                P_bottle, m_bottle, c_d, A, mu, T_atm, P_atm)
            self.vx = vi(1);
            self.vy = vi(2);
            self.vz = vi(3);
            self.x = ri(1);
            self.y = ri(2);
            self.z = ri(3);
            self.initialheading = [0 cos(theta) sin(theta)];
            
            self.m_air = m_air;
            self.m_air_i = m_air;
            self.m_water = m_water;
            self.V_air = V_air;
            self.V_air_i = V_air;
            self.P_bottle = P_bottle;
            
            self.m_bottle = m_bottle;
            self.c_d = c_d;
            self.A = A;
            self.mu = mu;
            
            self.T_atm = T_atm + 273.15;
            self.P_atm = P_atm;
            self.rho_atm = P_atm / (self.R * T_atm);
        end
        function setstage(self)
            switch self.stage
                case 1
                    %water thrusting phase
                    if self.V_air >= self.volume
                        self.stage = self.stage+1;
                    end
                case 2
                    %air thrusting phase, choked flow
                    if self.Pstar < self.P_atm
                        self.stage = self.stage+1;
                    end
                case 3
                    %air thrusting phase, nonchoked flow
                    if self.P_bottle <= self.P_atm
                        self.stage = self.stage+1;
                    end
                case 4
                    % This is the ballistic phase and does not end until the
                    % integration ends.
            end
        end
        function [p, pstar, p_t] = pressure(self)
            switch self.stage
                case 1
                    p = p_i * (self.V_air_i/self.V_air) ^ self.gamma;
                    p_t = p_atm;
                    pstar = Inf;
                case 2
                    p_end = self.P_bottle_i * ...
                        (self.V_air_i / self.volume) ^ self.gamma;
                    p = (p_end * (self.m_air / self.m_air_i) ^ self.gamma);
                    pstar = p * (2 / (self.gamma + 1)) ^ ...
                        (self.gamma / (self.gamma - 1));
                    p_t = pstar;
                case 3
                    p_end = self.P_bottle_i * ...
                        (self.V_air_i / self.volume) ^ self.gamma;
                    p = (p_end * (self.m_air / self.m_air_i) ^ self.gamma);
                    pstar = p * (2 / (self.gamma + 1)) ^ ...
                        (self.gamma / (self.gamma - 1));
                    p_t = p_atm;
                case 4
                    p = p_atm;
                    pstar = Inf;
                    p_t = p_atm;
            end
        end
        function dmdt = madot(self)
            [p,pstar]=pressurefn(t,vars);
            switch self.stage
                case 1
                    dmdt = 0;
                case 2
                    rho = self.m_air/self.volume;
                    T = p/(rho*R);
                    T_t = 2*T/(gamma+1);
                    rho_t = pstar/(R*T_t);
                    v_t = sqrt(gamma*R*T_t);
                    dmdt = -rho_t * v_t * self.A_throat * c_d;
                case 3
                    rho=m_air/Vbottle;
                    T=p/(rho*R);
                    T_t=T*(p/p_atm)^((gamma-1)/gamma);
                    rho_t=p_atm/(R*T_t);
                    M_t= sqrt( (T_t/T - 1) * (2/(gamma-1)) );
                    v_t=M_t*sqrt(gamma*R*T_t);
                    dmdt = -rho_t*v_t*A*c_d;
                case 4
                    dmdt = 0;
            end
        end
        function dmdt = mwdot(self)
            switch self.stage
                case 1
                    dmdt = -self.discharge * self.rho_water * self.A * ...
                        sqrt(2*(self.P_bottle - self.P_atm) / ...
                        self.rho_water);
                case 2
                    dmdt = 0;
                case 3
                    dmdt = 0;
                case 4
                    dmdt = 0;
            end
        end
        function Vdot(self)
            
        end
        function DYDT = derivatives(self, t, vars)
            % output    input
            %   vx       ax
            %   vy       ay
            %   vz       az
            %   x        vx
            %   y        vy
            %   z        vz
            %   m_a      madot
            %   m_w      mwdot
            %   V_a      Vdot
            
            self.update(t, vars)
            self.setstage()
            
            a = self.vdot();
            
            DYDT = [ a(1);
                     a(2);
                     a(3);
                     vars(1);
                     vars(2);
                     vars(3);
                     self.madot();
                     self.mwdot();
                     self.Vdot()
                   ];
        end
        function update(self, t, vars)
            l = length(self.vx) + 1;
            self.vx(l) = vars(1);
            self.vy(l) = vars(2);
            self.vz(l) = vars(3);
            self.x(l) = vars(4);
            self.y(l) = vars(5);
            self.z(l) = vars(6);
            self.m_air(l) = vars(7);
            self.m_water(l) = vars(8);
            self.V_air(l) = vars(9);
            [self.P_bottle(l), self.Pstar(l)] = self.pressure();
            self.t(l) = t;
        end
        function T = thrust(self)
            [p, pstar, p_t] = self.pressure();
            dir = self.normalizeV();            
            switch self.stage
                case 1
                    Tmag = 2 * self.discharge * (p - P_atm) * self.A_throat;
                case 2
                    rho = self.m_air / self.volume;
                    t = p / (rho * R);
                    t_t = 2 * t / (self.gamma + 1);
                    v_t = sqrt(gamma * R * t_t);
                    Tmag = v_t * self.madot() + (pstar - p_atm) * self.A_throat;
                case 3
                    rho = self.m_air / self.volume;
                    t = p / (rho * R);
                    t_t = t * (p / p_atm) ^ ((gamma - 1) / gamma);
                    M_t = sqrt((t_t / t - 1) * (2 / gamma - 1));
                    v_t = M_t * sqrt(gamma * R * t_t);
                    Tmag = self.madot() * v_t;
                case 4
                    Tmag = 0;
            end
            T = Tmag * dir;
        end
        function m = mass(self)
            m = self.m_bottle + self.m_air + self.m_water;
        end
        function [opts, vars, titles, units] = maps(self)
            opts = {'x', 'crossrange';
                    'y', 'distance';
                    'z', 'height';
                    'vx', 'x velocity';
                    'vy', 'y velocity';
                    'vz', 'z velocity';
                    'v', 'speed';
                    'm', 'mass';
                    'vol', 'volume';
                    'm_a', 'air mass';
                    'm_w', 'water mass';
                    't', 'time'};
            vars = {self.x;
                    self.y;
                    self.z;
                    self.vx;
                    self.vy;
                    self.vz;
                    self.v;
                    self.mass();
                    self.V_air;
                    self.m_air;
                    self.m_water;
                    self.t};
            titles = {'Cross Range Distance';
                      'Range Distance';
                      'Height';
                      'Cross Range Velocity';
                      'Range Velocity';
                      'Vertical Velocity';
                      'Speed';
                      'Rocket Mass';
                      'Air Volume';
                      'Air Mass';
                      'Water Mass';
                      'Time'};
            units = {'m';
                     'm';
                     'm';
                     'm/s';
                     'm/s';
                     'm/s';
                     'm/s';
                     'kg';
                     'm^3'
                     'kg';
                     'kg';
                     's'};
        end
    end
end