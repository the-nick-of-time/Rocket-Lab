function montecarlo(rocket, varargin)
if mod(nargin - 1, 3)
    error('You must define var:mean:sigma triplets')
end
opts = { 'p';
         'm_w';
         'T';
         'c_d';
         'theta';
         'mu';
         'p_atm';
         'wind';
         'mu';
         'A'
       };

for i = 2:3:nargin
    
end
end