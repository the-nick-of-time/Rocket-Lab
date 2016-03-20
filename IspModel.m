classdef IspModel < BottleRocket
    properties(Access = {?BottleRocket})
        staticdata
        type = 1 + rand(1);
    end
    methods
        function self = IspModel(data, ri, theta, m_bottle, c_d, A, mu,...
                T_atm, P_atm)
            self.staticdata = data;
            self.initialheading = [0 cos(theta) sin(theta)];
            vi = self.staticdata.deltaV() * self.initialheading;
            self.vx = vi(1);
            self.vy = vi(2);
            self.vz = vi(3);
            self.x = ri(1);
            self.y = ri(2);
            self.z = ri(3);
            
            self.m_bottle = m_bottle;
            self.c_d = c_d;
            self.A = A;
            self.mu = mu;
            
            self.T_atm = T_atm;
            self.P_atm = P_atm;
            self.rho_atm = P_atm / (self.R * T_atm);
        end
        function m = mass(self)
            m = self.m_bottle + self.rho_atm * self.volume;
        end
        function T = thrust(~)
            T = [0 0 0];
        end
        function update(self, t, vars)
            l = length(self.vx) + 1;
            self.vx(l) = vars(1);
            self.vy(l) = vars(2);
            self.vz(l) = vars(3);
            self.x(l) = vars(4);
            self.y(l) = vars(5);
            self.z(l) = vars(6);
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
            %  mdot            m
            
            self.update(t, vars);
            
            a = self.vdot();
            
            DYDT = [ a(1);
                     a(2);
                     a(3);
                     vars(1);
                     vars(2);
                     vars(3);
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
                    't', 'time'};
            vars = {self.x;
                    self.y;
                    self.z;
                    self.vx;
                    self.vy;
                    self.vz;
                    self.v;
                    self.t};
            titles = {'Cross Range Distance';
                      'Range Distance';
                      'Height';
                      'Cross Range Velocity';
                      'Range Velocity';
                      'Vertical Velocity';
                      'Speed';
                      'Time'};
            units = {'m';
                     'm';
                     'm';
                     'm/s';
                     'm/s';
                     'm/s';
                     'm/s';
                     's'};
        end
    end
end