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
            startpoint = find(diff(self.force) > 12, 1);
            endpoint = find(self.force == min(self.force(startpoint+1:end)));
            indices = startpoint:endpoint;
        end
        function rv = getForce(self, which)
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
            range = self.isolate();
            ind = 1:length(range);
            F = self.force(range);
            time = ind / self.freq;
            DF = F(end) - F(1);
            DT = length(range) / self.freq;
            base = (DF/DT) * time + F(1);
%             hold on
%             plot(time, base)
            
            I = 0;
%             for x = ind
%                 % left Riemann sum 
%                 I = I + (F(x) - base(x)) * (1 / self.freq);
%             end
            for x = ind(1:end-1)
                % trapezoidal Riemann sum
                I = I + .5 * (F(x) - base(x) + F(x+1) - ...
                    base(x+1)) * (1 / self.freq);
            end
            % Assumes that the mass of the air is negligible
            Isp = I / (self.watermass * 9.807);
        end
        function F = interpolate(self, time)
            desiredindex = time * self.freq;
            indices = self.isolate();
            if desiredindex > (indices(end) - indices(1))
                F = 0;
            else
                fl = self.force(floor(desiredindex));
                fu = self.force(ceil(desiredindex));
                rm = mod(desiredindex, 1);
                F = fl + (fu - fl) / rm;
            end
        end
        function makeplot(self, everything)
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
            rv = max(diff(self.force));
        end
        function DV = deltaV(self)
            DV = self.getImpulse() * self.g * ...
                log((self.emptymass + self.watermass) / self.emptymass);
        end
    end
end