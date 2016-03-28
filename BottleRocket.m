classdef (Abstract=true) BottleRocket < handle
    properties%(Access = 'protected')
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
        t = 0;
    end
    methods
        function D = drag(self)
            [dir, mag] = self.normalizeV(true);
            D = .5 * self.rho_atm * mag^2 * self.c_d * self.A * (-dir);
        end
        function dvdt = vdot(self, t, vars)
            dir = self.relativeV();
            dvdt = (self.thrust(t, vars) + self.drag() + self.weight()) ...
                / self.mass(true);
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
                if norm([self.x self.y self.z]) < 1
                    vhat = self.initialheading;
                elseif vmag(i) > 0
                    vhat(i,:) = v(i,:) / vmag(i);
                else
                    vhat(i,:) = [0 0 0];
                end
            end
        end
        function f = friction(self)
            railLength = 1; %measure and change later
            if norm([self.x self.y self.z]) < railLength
                f = self.mu * dot(self.weight(), self.initialheading) * ...
                    self.initialheading;
            else
                f = [0 0 0];
            end
        end
        function w = windV(self)
            if self.z(end) < 1
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
                    w = v1 + self.z(end) * (dv/dh);
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
            W = self.mass(true) * self.g;
        end
        function [value, isterminal, direction] = endcondition(~, ~, vars)
            value = vars(6);
            isterminal = 1;
            direction = -1;
        end
        function fig = makeplot(self, xvar, yvar, figureargs, plotargs, outfig)
            if strcmpi(xvar, yvar)
                error('Variables must be different')
            end
            
            [opts, vars, titles, units] = self.maps();
            
            [tempx, ~] = find(strcmpi(xvar, opts));
            [tempy, ~] = find(strcmpi(yvar, opts));
            if ~numel(tempx) || ~numel(tempy)
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
            
            print(fig, '-dpng', ['./figures/' TITLE '.png'])
        end
        function rv = type(self)
            rv = self.ID;
        end
        function rv = get(self, which)
            [opts, vars] = self.maps();
            [where, ~] = find(strcmpi(which, opts));
            rv = vars{where};
        end
    end
    methods (Abstract)
        thrust(self, t, vars) %outputs thrust as a vector along the correct direction
        mass(self, onlyend) %outputs mass of the whole system as kg
        maps(self) %does the plotting thing
        derivatives(self, t, vars) %returns vector of derivatives for ode45
        update(self, t, vars) %updates self variables after an ode45 step
        initialconditions(self) %returns the initial condition vector
        finalize(self, t, vars) %performs final update of integration vars
    end
end