classdef ThrustModel < BottleRocket
    properties(Access = {?BottleRocket})
        m_fuel
        staticdata
        ID = 2 + rand(1);
    end
    methods
        function self = ThrustModel(staticdata, ri, theta, m_bottle,...
                c_d, A, mu, T_atm, P_atm, m_water, P_bottle, wind)
            % Class constructor
            % Inputs:
            %   staticdata: the statictest object holding the data from the
            %       static test stand
            %   ri: initial position vector, should be [0 0 height]
            %       definitionally
            %   theta: launch angle, set by stand
            %   m_bottle: mass of the empty bottle
            %   c_d: drag coefficient
            %   A: cross-sectional area of bottle
            %   mu: coefficient of friction between bottle and launch stand
            %   T_atm: temperature of ambient air
            %   P_atm: atmospheric pressure
            %   m_water: initial mass of water
            %   P_bottle: initial pressure inside bottle
            %   wind: matrix giving wind data
            % Outputs:
            %   None explicit
            self.vx = 0;
            self.vy = 0;
            self.vz = 0;
            self.x = ri(1);
            self.y = ri(2);
            self.z = ri(3);
            self.initialheading = [0 cosd(theta) sind(theta)];
            
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
            self.wind_data = wind;
        end
        function m = mass(self, all)
            % Calculates mass of rocket
            % Inputs:
            %   none
            % Outputs:
            %   m: Total mass of rocket
            if ~all
                m = self.m_fuel(end) + self.m_bottle;
            else
                m = self.m_fuel + self.m_bottle;
            end
        end
        function dmdt = mdot(self)
            % Returns mdot for integration
            % Inputs:
            %   none
            % Outputs:
            %   dmdt: time rate of change of mass
            %
            % possibilities: 
            %   1. mdot is assumed constant throughout the duration of
            %   thrusting. this would need to verify the current value of T
            %   2. exhaust velocity is assumed constant. mdot can be
            %   extracted from thrust data.
            dmdt = self.staticdata.masschange();
        end
        function T = thrust(self)
            % Calculates thrust along the direction of the rocket's motion
            % relative to the wind
            % Inputs:
            %   none
            % Outputs:
            %   T: thrust vector
            dir = self.relativeV();
            T = self.staticdata.interpolate(self.t(end)) * dir;
%             disp([self.t(end) T])
        end
        function update(self, t, vars)
            % Updates internal variables after an integration step
            % Inputs:
            %   t: time of integration
            %   vars: variables of integration
            % Outputs:
            %   None explicit
            l = length(self.vx) + 1;
            self.vx(l,1) = vars(1);
            self.vy(l,1) = vars(2);
            self.vz(l,1) = vars(3);
            self.x(l,1) = vars(4);
            self.y(l,1) = vars(5);
            self.z(l,1) = vars(6);
            self.m_fuel(l,1) = vars(7);
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
            %  mdot           m_f
            %
            % Creates the vector of derivatives for ode45, as defined above
            % Inputs:
            %   t: time of integration
            %   vars: variables of integration
            % Outputs:
            %   DYDT: vector of derivatives
            
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
            % Creates name:value maps for usage by the makeplot functions
            % Inputs:
            %   none
            % Outputs:
            %   opts: string pairs identifying which variable is which
            %   vars: data vectors
            %   titles: official names associated with variables
            %   units: units of variable
            opts = {'x', 'crossrange';
                    'y', 'distance';
                    'z', 'height';
                    'vx', 'x velocity';
                    'vy', 'y velocity';
                    'vz', 'z velocity';
                    'v', 'speed';
                    'm', 'mass';
                    't', 'time'};
            [~, v] = self.normalizeV(false);
            vars = {self.x;
                    self.y;
                    self.z;
                    self.vx;
                    self.vy;
                    self.vz;
                    v;
                    self.mass(true);
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
        function rv = initialconditions(self)
            % Creates initial conditions vector
            % Inputs:
            %   none
            % Outputs:
            %   rv: initial conditions for integration
            rv = [ self.vx;
                   self.vy;
                   self.vz;
                   self.x;
                   self.y;
                   self.z;
                   self.m_fuel
                 ];
        end
        function finalize(self, t, vars)
            % Overwrites the raw values of variables from the integration
            % with cleaned, post-processed versions that are ode45's return
            % value.
            % Inputs:
            %   none
            % Outputs:
            %   None explicit
            self.vx = vars(:,1);
            self.vy = vars(:,2);
            self.vz = vars(:,3);
            self.x = vars(:,4);
            self.y = vars(:,5);
            self.z = vars(:,6);
            self.m_fuel = vars(:,7);
            self.t = t;
        end
    end
end