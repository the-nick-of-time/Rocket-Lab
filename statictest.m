classdef statictest
    properties(Access='protected')
        leftforce
        rightforce
        force
        freq
        emptymass
        watermass
        g = 9.81;
    end
    methods
        function self = statictest(filename, samplefreq, emptymass, ...
                fuelmass)
            % Class constructor
            % Inputs:
            %   filename: name of file to read data from
            %   samplefreq: sampling frequency for the data
            %   emptymass: mass of empty bottle
            %   fuelmass: mass of water and air within bottle
            % Outputs:
            %   
            if nargin > 0
                data = load(filename);
                self.leftforce = data(:,1) * 4.44822;
                self.rightforce = data(:,2) * 4.44822;
                self.force = data(:,3) * 4.44822;
                self.freq = samplefreq;
                self.emptymass = emptymass;
                self.watermass = fuelmass;
            end
        end
        function indices = isolate(self)
            % Finds the indices at which thrust is happening
            % Inputs:
            %   none
            % Outputs:
            %   indices: the vector of relevant indices
            startpoint = find((self.force) > 20, 1);
            [~, highpoint] = max(self.force);
            [~, relativeendpoint] = min(self.force(highpoint:end));
            indices = startpoint:(startpoint+relativeendpoint);
        end
        function rv = getForce(self, which)
            % Gets force data
            % Inputs:
            %   none
            % Outputs:
            %   which: left, right, or total
            if nargin == 1 || strcmpi(which, 't') || ...
                    strcmpi(which, 'total')
                rv = self.totalforce;
                return
            end
            if strcmpi(which, 'l') || strcmpi(which, 'left')
                rv = self.leftforce;
                return
            end
            if strcmpi(which, 'r') || strcmpi(which, 'right')
                rv = self.rightforce;
                return
            end
        end
        function Isp = getImpulse(self)
            % Integrates the thrust curve to get specific impulse
            % Inputs:
            %   none
            % Outputs:
            %   Isp: specific impulse imparted
            range = self.isolate();
            ind = 1:length(range);
            F = self.force(range);
            time = ind / self.freq;
            DF = F(end) - F(1);
            DT = length(range) / self.freq;
            base = (DF/DT) * time + F(1);
            
            %I = 0;
%             for x = ind
%                 % left Riemann sum 
%                 I = I + (F(x) - base(x)) * (1 / self.freq);
%             end
%             for x = ind(1:end-1)
%                 % trapezoidal Riemann sum
%                 I = I + .5 * (F(x) - base(x) + F(x+1) - ...
%                     base(x+1)) * (1 / self.freq);
%             end
            I = trapz(F - base') / self.freq;
            % Assumes that the mass of the air is negligible
            Isp = I / (self.watermass * self.g);
        end
        function F = interpolate(self, time)
            % Finds the thrust at any given time by linear interpolation
            % Inputs:
            %   none
            % Outputs:
            %   F: thrust force
            desiredindex = time * self.freq;
            indices = self.isolate();
            first = indices(1);
            DF = self.force(indices(end)) - self.force(1);
            DT = length(indices) / self.freq;
            base = (DF/DT) * time;
            if desiredindex > length(indices)
                F = 0;
            else
                fl = self.force(floor(desiredindex + first));
                fu = self.force(ceil(desiredindex + first));
                rm = mod(desiredindex, 1);
                F = fl + (fu - fl) * rm - base;
            end
        end
        function makeplot(self, everything)
            % Plots thrust against time
            % Inputs:
            %   none
            % Outputs:
            %   none
            figure

            if nargin == 1 || ~everything
                ind = self.isolate();
                time = (ind - ind(1))/self.freq;
                plot(time, self.force(ind))
            else
                time = (1:length(self.force))/self.freq;
                plot(time, self.force)
            end
            title('Time against Force')
            xlabel('Time (s)')
            ylabel('Force (N)')
        end
        function rv = change(self)
            % Finds the spike value of the thrust
            % Inputs:
            %   none
            % Outputs:
            %   rv
            rv = max(diff(self.force));
        end
        function DV = deltaV(self)
            % Finds the total velocity imparted to the rocket
            % Inputs:
            %   none
            % Outputs:
            %   DV: the delta V of the rocket by the ideal rocket equation
            DV = self.getImpulse() * self.g * ...
                log((self.emptymass + self.watermass) / self.emptymass);
        end
        function DmDt = masschange(self)
            % Finds the average rate of change of mass of the rocket
            % Inputs:
            %   none
            % Outputs:
            %   DmDt: average mass change rate 
            ind = self.isolate();
            Dt = length(ind) / self.freq;
            Dm = -self.watermass;
            
            DmDt = Dm / Dt;
        end
    end
end