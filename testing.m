clear
close all
clc
% objects = [statictest];
% for i = 4:-1:1
%     name = sprintf('./STATIC/LAtest%d',i);
%     objects(i) = statictest(name, 1.652e3, .2, 1);
%     objects(i).makeplot()
%     objects(i).getImpulse()
% end

static = statictest('./STATIC/LAtest3', 1.652e3, .2, 1);
initialpos = [0 0 .1]; %m
theta = 45; %degrees
m_bottle = .2; %kg
c_d = .4; %drag
A = .008659; %m^2
mu = .2; %friction with the launch rails
T_atm = 20; %degrees C
P_atm = 84; %kPa
wind = [0 0 0 0;0 0 0 1]; %3-vector velocity m/s then altitude m

rocket = IspModel(static, initialpos, theta, m_bottle, c_d, A, mu, ...
    T_atm, P_atm, wind);


m_water = 1; %kg
P_bottle = 40 * 6894.76; %input psi output Pa
% rocket = ThrustModel(static, initialpos, theta, m_bottle, c_d, A, mu, ...
%     T_atm, P_atm, m_water, P_bottle, wind);


integrate(rocket)

rocket.makeplot('t', 'z', {'Position', [100 100 900 700], ...
    'DefaultAxesFontSize', 16}, {});
rocket.makeplot('t', 'v', {'Position', [100 100 900 700], ...
    'DefaultAxesFontSize', 16}, {});
rocket.makeplot('y', 'z', {'Position', [100 100 900 700], ...
    'DefaultAxesFontSize', 16}, {});