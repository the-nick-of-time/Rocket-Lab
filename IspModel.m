classdef IspModel < BottleRocket
    properties(Access = {?BottleRocket})
        staticdata
        ID = 1 + rand(1);
    end
    methods
        function self = IspModel(data, ri, theta, m_bottle, c_d, A, mu,...
                T_atm, P_atm, wind)
            self.staticdata = data;
            self.initialheading = [0 cosd(theta) sind(theta)];
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
            
            self.T_atm = T_atm + 273.15;
            self.P_atm = P_atm;
            self.rho_atm = self.P_atm / (self.R * self.T_atm);
            self.wind_data = wind;
        end
        function m = mass(self, ~)
            m = self.m_bottle + self.rho_atm * self.volume;
        end
        function T = thrust(~)
            T = [0 0 0];
        end
        function update(self, t, vars)
            l = length(self.vx) + 1;
            self.vx(l,1) = vars(1);
            self.vy(l,1) = vars(2);
            self.vz(l,1) = vars(3);
            self.x(l,1) = vars(4);
            self.y(l,1) = vars(5);
            self.z(l,1) = vars(6);
            self.t(l,1) = t;
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
            [~, v] = self.normalizeV(false);
            vars = {self.x;
                    self.y;
                    self.z;
                    self.vx;
                    self.vy;
                    self.vz;
                    v;
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
        function rv = initialconditions(self)
            rv = [ self.vx;
                   self.vy;
                   self.vz;
                   self.x;
                   self.y;
                   self.z
                 ];
        end
        function finalize(self, t, vars)
            self.vx = vars(:,1);
            self.vy = vars(:,2);
            self.vz = vars(:,3);
            self.x = vars(:,4);
            self.y = vars(:,5);
            self.z = vars(:,6);
            self.t = t;
        end
    end
end