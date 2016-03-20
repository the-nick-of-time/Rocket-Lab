classdef ThrustModel < BottleRocket
    properties(Access = {?BottleRocket})
        m_fuel
        staticdata
        type = 2 + rand(1);
    end
    methods
        function self = Thrustmodel(staticdata, vi, ri, theta, m_bottle,...
                c_d, A, mu, T_atm, P_atm, m_water, P_bottle)
            self.vx = vi(1);
            self.vy = vi(2);
            self.vz = vi(3);
            self.x = ri(1);
            self.y = ri(2);
            self.z = ri(3);
            self.initialheading = [0 cos(theta) sin(theta)];
            
            self.staticdata = staticdata;
            self.m_fuel = m_water + (P_bottle / (T_atm * self.R)) * ...
                (self.volume - m_water / 1000);
            
            self.m_bottle = m_bottle;
            self.c_d = c_d;
            self.A = A;
            self.mu = mu;
            
            self.T_atm = T_atm;
            self.P_atm = P_atm;
            self.rho_atm = P_atm / (self.R * T_atm);
        end
        function m = mass(self)
            m = self.m_fuel + self.m_bottle;
        end
        function dmdt = mdot(self)
            % possibilities: 
            %   1. mdot is assumed constant throughout the duration of
            %   thrusting. this would need to verify the current value of T
            %   2. exhaust velocity is assumed constant. mdot can be
            %   extracted from thrust data.
        end
        function T = thrust(self, t, ~)
            dir = self.normalizeV(true);
            T = self.staticdata.interpolate(t) * dir;
        end
        function update(self, t, vars)
            l = length(self.vx) + 1;
            self.vx(l) = vars(1);
            self.vy(l) = vars(2);
            self.vz(l) = vars(3);
            self.x(l) = vars(4);
            self.y(l) = vars(5);
            self.z(l) = vars(6);
            self.m_fuel(l) = vars(7);
            self.t(l) = t;
        end
        function DYDT = derivatives(self, t, vars)
            % Elements of the vector:
            % Inputs        Outputs
            %   ax             vx
            %   ay             vy
            %   az             vz
            %   vx             x
            %   vy             y
            %   vz             z
            %  mdot           m_f
            
            self.update(t, vars);
            
            a = self.vdot();
            
            DYDT = [ a(1);
                     a(2);
                     a(3);
                     vars(1);
                     vars(2);
                     vars(3);
                     self.mdot();
                   ];
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
                    't', 'time'};
            vars = {self.x;
                    self.y;
                    self.z;
                    self.vx;
                    self.vy;
                    self.vz;
                    self.v;
                    self.mass();
                    self.t};
            titles = {'Cross Range Distance';
                      'Range Distance';
                      'Height';
                      'Cross Range Velocity';
                      'Range Velocity';
                      'Vertical Velocity';
                      'Speed';
                      'Rocket Mass'
                      'Time'};
            units = {'m';
                     'm';
                     'm';
                     'm/s';
                     'm/s';
                     'm/s';
                     'm/s';
                     'kg';
                     's'};
        end
    end
end