classdef (Abstract=true) BottleRocket
    properties(Access = 'protected')
        % position and velocity
        vx
        vy
        vz
        x
        y
        z
        initialheading
        % rocket parameters
        m_bottle
        c_d
        A
        mu
        volume = 2048e-6;
        % atmospheric parameters
        T_atm
        P_atm
        rho_atm
        wind_data
        % universal parameters
        g = [0 0 -9.81];
        R = 287;
        gamma = 1.4;
        rho_water = 1000;
        % time from integration
        t
    end
    methods
        function D = drag(self)
            [dir, mag] = self.normalizeV(true);
            D = .5 * self.rho_atm * mag^2 * self.c_d * self.A * (-dir);
        end
        function dvdt = vdot(self)
            dir = normalizeV(true);
            dvdt = (self.thrust() + self.drag() + ...
                dot(self.weight(), dir) * dir) / self.mass();
        end
        function [vhat, vmag] = normalizeV(self, onlyLatest)
            if onlyLatest
                v = [self.vx(end) self.vy(end) self.vz(end)];
            else
                v = [self.vx self.vy self.vz];
                vmag = zeros(length(self.vx),1);
                vhat = zeros(length(self.vx),3);
            end
            
            for i = 1:size(v, 1)
                vmag(i) = norm(v(i,:));
                if vmag(i) > 0
                    vhat(i,:) = v(i,:) / vmag(i);
                else
                    vhat(i,:) = [0 0 0];
                end
            end
        end
        function f = friction(self)
            railLength = 1; %measure and change later
            if norm([self.x self.y self.z]) < railLength
                f = self.mu * dot(self.weight(), self.initialheading);
            else
                f = [0 0 0];
            end
        end
        function w = windV(self)
            if self.z < 1
                % near the ground the wind speed is zero. This is primarily
                % to account for the impossibility of reorientation on the
                % stand.
                w = [0 0 0];
            else
                v1 = self.wind_data(1,1:3);
                h1 = self.wind_data(1,4);
                v2 = self.wind_data(2,1:3);
                h2 = self.wind_data(2,4);
                
                dv = v2-v1;
                dh = h2-h1;
                % assume that the wind speed varies linearly over the 
                % entire height of the flight (save that first meter)
                if dh ~= 0
                    w = v1 + self.z * (dv/dh);
                else
                    w = [0 0 0];
                end
            end
        end
        function [vhat, vmag] = relativeV(self)
            rv = [self.vx(end) self.vy(end) self.vz(end)];
            wv = self.windV();
            
            v_r = rv - wv;
            vmag = norm(v_r);
            vhat = v_r / vmag;
        end
        function W = weight(self)
            W = self.mass() * self.g;
        end
        function [value, isterminal, direction] = endcondition(~, ~, vars)
            value = vars(3);
            isterminal = 1;
            direction = -1;
        end
        function fig = makeplot(self, xvar, yvar, figureargs, plotargs, outfig)
            if strcmpi(xvar, yvar)
                error('Variables must be different')
            end
            
            [opts, vars, titles, units] = self.maps(xvar, yvar);
            
            tempx = find(sum(strcmpi(xvar, opts), 2));
            tempy = find(sum(strcmpi(yvar, opts), 2));
            if ~any(tempx) || ~any(tempy)
                error('Use valid variable identifiers')
            else
                X = vars{tempx};
                Y = vars{tempy};
                TITLE = sprintf('%s against %s', titles{tempx}, titles{tempy});
                XLABEL = sprintf('%s (%s)', titles{tempx}, units{tempx});
                YLABEL = sprintf('%s (%s)', titles{tempy}, units{tempy});
            end
            
            if nargin == 6
                fig = outfig;
            else
                fig = figure(figureargs{:});
            end
            figure(fig)
            plot(X, Y, plotargs{:})
            title(TITLE)
            xlabel(XLABEL)
            ylabel(YLABEL)
        end
        function rv = Type(self)
            rv = self.type;
        end
    end
    methods (Abstract)
        thrust(self, t, vars) %outputs thrust as a vector along the correct direction
        mass(self) %outputs mass of the whole system as kg
        maps(self) %does the plotting thing
        derivatives(self, t, vars) %returns vector of derivatives for ode45
        update(self, t, vars) %updates self variables after an ode45 step
%         Type(self) %returns a double identifying the object
    end
end