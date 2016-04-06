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
        A_throat = pi * (21.67e-3 / 2) ^ 2;
        discharge;
        stage = 1;
        ID = 3 + rand(1);
        thrustvec = 0
    end
    methods
        function self = ThermoModel(ri, theta, m_water, V_air,...
                P_bottle, m_bottle, c_d, A, mu, T_atm, P_atm, ...
                c_discharge, wind)
            self.vx = 0;
            self.vy = 0;
            self.vz = 0;
            self.x = ri(1);
            self.y = ri(2);
            self.z = ri(3);
            self.initialheading = [0 cosd(theta) sind(theta)];
            
            self.m_water = m_water;
            self.V_air = V_air;
            self.V_air_i = V_air;
            self.P_bottle = P_bottle;
            self.P_bottle_i = P_bottle;
            %assumes that the bottle is in thermal equilibrium with outside
            self.m_air = (P_bottle * V_air)/(self.R * (T_atm + 273.15));
            self.m_air_i = self.m_air;
            
            self.m_bottle = m_bottle;
            self.c_d = c_d;
            self.A = A;
            self.mu = mu;
            self.discharge = c_discharge;
            
            self.T_atm = T_atm;
            self.P_atm = P_atm;
            self.rho_atm = P_atm / (self.R * T_atm);
            self.wind_data = wind;
        end
        function setstage(self)
            switch self.stage
                case 1
                    %water thrusting phase
                    if self.V_air(end) >= self.volume
                        self.stage = self.stage+1;
                    end
                case 2
                    %air thrusting phase, choked flow
                    if self.Pstar(end) < self.P_atm
                        self.stage = self.stage+1;
                    end
                case 3
                    %air thrusting phase, nonchoked flow
                    if self.P_bottle(end) <= self.P_atm
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
                    p = self.P_bottle_i * (self.V_air_i/self.V_air(end)) ^ self.gamma;
                    p_t = self.P_atm;
                    pstar = Inf;
                case 2
                    p_end = self.P_bottle_i * ...
                        (self.V_air_i / self.volume) ^ self.gamma;
                    p = (p_end * (self.m_air(end) / self.m_air_i) ^ self.gamma);
                    pstar = p * (2 / (self.gamma + 1)) ^ ...
                        (self.gamma / (self.gamma - 1));
                    p_t = pstar;
                case 3
                    p_end = self.P_bottle_i * ...
                        (self.V_air_i / self.volume) ^ self.gamma;
                    p = (p_end * (self.m_air(end) / self.m_air_i) ^ self.gamma);
                    pstar = p * (2 / (self.gamma + 1)) ^ ...
                        (self.gamma / (self.gamma - 1));
                    p_t = self.P_atm;
                case 4
                    p = self.P_atm;
                    pstar = Inf;
                    p_t = self.P_atm;
            end
        end
        function dmdt = madot(self)
            %[p,pstar]=pressurefn(t,vars);
            switch self.stage
                case 1
                    dmdt = 0;
                case 2
                    rho = self.m_air(end) / self.V_air(end);
                    T = self.P_bottle(end) / (rho * self.R);
                    T_t = 2 * T / (self.gamma + 1);
                    if self.Pstar(end) == Inf
                        rho_t = self.P_bottle(end) / (self.R * T_t);
                    else
                        rho_t = self.Pstar(end) / (self.R * T_t);
                    end
                    v_t = sqrt(self.gamma * self.R * T_t);
                    dmdt = -rho_t * v_t * self.A_throat * self.discharge;
                case 3
                    rho = self.m_air(end) / self.volume;
                    T = self.P_bottle(end) / (rho * self.R);
                    T_t = T * (self.P_bottle(end) / self.P_atm) ^ ...
                        ((self.gamma - 1) / self.gamma);
                    rho_t = self.P_atm / (self.R * T_t);
                    M_t = sqrt( (T_t/T - 1) * (2/(self.gamma-1)) );
                    v_t = M_t * sqrt(self.gamma * self.R * T_t);
                    dmdt = -rho_t * v_t * self.A_throat * self.discharge;
                case 4
                    dmdt = 0;
            end
        end
        function dmdt = mwdot(self)
            switch self.stage
                case 1
                    dmdt = -self.discharge * self.rho_water * ...
                        self.A_throat * sqrt(2*(self.P_bottle(end) - ...
                        self.P_atm) / self.rho_water);
                otherwise
                    dmdt = 0;
            end
        end
        function dVdt = Vdot(self)
            switch self.stage
                case 1
                    dVdt = self.discharge * self.A_throat * ...
                        sqrt((2/self.rho_water) * ((self.P_bottle_i * ...
                        (self.V_air_i/self.V_air(end)) ^ self.gamma) - ...
                        self.P_atm));
                otherwise
                    dVdt = 0;
            end
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
            self.vx(l,1) = vars(1);
            self.vy(l,1) = vars(2);
            self.vz(l,1) = vars(3);
            self.x(l,1) = vars(4);
            self.y(l,1) = vars(5);
            self.z(l,1) = vars(6);
            self.m_air(l,1) = vars(7);
            self.m_water(l,1) = vars(8);
            self.V_air(l,1) = vars(9);
            [self.P_bottle(l,1), self.Pstar(l,1)] = self.pressure();
            self.t(l,1) = t;
        end
        function T = thrust(self)
            %[p, pstar] = self.pressure();
            dir = self.relativeV();            
            switch self.stage
                case 1
                    Tmag = 2 * self.discharge * (self.P_bottle(end) - self.P_atm(end)) * self.A_throat;
                case 2
                    rho = self.m_air(end) / self.V_air(end);
                    t = self.P_bottle(end) / (rho * self.R);
                    t_t = 2 * t / (self.gamma + 1);
                    v_t = sqrt(self.gamma * self.R * t_t);
                    if self.Pstar(end) == Inf
                        Tmag = v_t * self.madot() + (self.P_bottle(end) - self.P_atm(end)) * self.A_throat;
                    else
                        Tmag = v_t * self.madot() + (self.Pstar(end) - self.P_atm(end)) * self.A_throat;
                    end
                case 3
                    rho = self.m_air(end) / self.volume;
                    t = self.P_bottle(end) / (rho * self.R);
                    t_t = t * (self.P_bottle(end) / self.P_atm(end)) ^ ((self.gamma - 1) / self.gamma);
                    M_t = sqrt((t_t / t - 1) * (2 / self.gamma - 1));
                    v_t = M_t * sqrt(self.gamma * self.R * t_t);
                    Tmag = self.madot() * v_t;
                case 4
                    Tmag = 0;
            end
            T = Tmag * dir;
            self.thrustvec = [self.thrustvec;Tmag];
        end
        function m = mass(self, all)
            if ~all
                m = self.m_bottle + self.m_air(end) + self.m_water(end);
            else
                m = self.m_bottle + self.m_air + self.m_water;
            end
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
            [~, v] = self.relativeV();
            vars = {self.x;
                    self.y;
                    self.z;
                    self.vx;
                    self.vy;
                    self.vz;
                    v;
                    self.mass(true);
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
        function rv = initialconditions(self)
            rv = [ self.vx;
                   self.vy;
                   self.vz;
                   self.x;
                   self.y;
                   self.z;
                   self.m_air;
                   self.m_water;
                   self.V_air
                 ];
        end
        function finalize(self, t, vars)
            self.vx = vars(:,1);
            self.vy = vars(:,2);
            self.vz = vars(:,3);
            self.x = vars(:,4);
            self.y = vars(:,5);
            self.z = vars(:,6);
            self.m_air = vars(:,7);
            self.m_water = vars(:,8);
            self.V_air = vars(:,9);
            self.t = t;
        end
    end
end