classdef ThermoModel < BottleRocket
    properties(Access = 'protected')
        m_air
        m_water
        V_air
    end
    methods
        function fig = makeplot(self, xvar, yvar, figureargs, plotargs, outfig)
            opts = {'x', 'distance';
                    'y', 'height';
                    'v', 'velocity';
                    'theta', 'angle';
                    'm', 'mass';
                    'vol', 'volume';
                    'm_a', 'air mass';
                    'm_w', 'water mass'};
            titles = {'Distance';
                      'Height';
                      'Velocity';
                      'Angle';
                      'Rocket Mass';
                      'Air Volume';
                      'Air Mass';
                      'Water Mass'};
            units = {'m';
                     'm';
                     'm/s';
                     'rad';
                     'kg';
                     'm^3'
                     'kg';
                     'kg'};
            vars = {self.x;
                    self.y;
                    self.v;
                    self.getMass();
                    self.V_air;
                    self.m_water;
                    self.m_air};
            
            tempx = find(sum(strcmpi(xvar, opts), 2));
            tempy = find(sum(strcmpi(yvar, opts), 2));
            if strcmpi(xvar, yvar)
                error('Variables must be different')
            elseif ~any(tempx) || ~any(tempy)
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
    end
end