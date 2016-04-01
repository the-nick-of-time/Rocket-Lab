clear
close all
clc
%%
% objects = [statictest];
% for i = 4:-1:1
%     name = sprintf('./STATIC/LAtest%d',i);
%     objects(i) = statictest(name, 1.652e3, .2, 1);
%     objects(i).makeplot()
%     objects(i).getImpulse()
%     objects(i).deltaV()
% end
%%
samplefreq = 1.652e3;
m_bottle = .085;
m_water_test = 1;
static = statictest('./STATIC/LAtest2', samplefreq, m_bottle, m_water_test);
% static.makeplot()

initialpos = [0 0 .1];
theta = 45;
c_d = .2;
A = .008659;
mu = .2;
T_atm = 20;
P_atm = 84000;
wind = [0 0 0 0;30 0 0 20];
isp = IspModel(static, initialpos, theta, m_bottle, c_d, A, mu,...
                T_atm, P_atm, wind);
integrate(isp)

% isp.makeplot('t', 'z', {'Position', [100 100 900 700]}, {});
% isp.makeplot('t', 'v', {'Position', [100 100 900 700]}, {});
isp.makeplot3d('x', 'y', 'z', {'Position', [100 100 900 700]}, {});

P_bottle = 40 * 6894.75729;
m_water = 1;
thrust = ThrustModel(static, initialpos, theta, m_bottle, c_d, A, mu,...
    T_atm, P_atm, m_water, P_bottle, wind);
integrate(thrust)

% thrust.makeplot('t', 'z', {'Position', [100 100 900 700]}, {});
% thrust.makeplot('t', 'v', {'Position', [100 100 900 700]}, {});
thrust.makeplot3d('x', 'y', 'z', {'Position', [100 100 900 700]}, {});

V_air = 1 * 1e-3;
c_discharge = .8;
thermo = ThermoModel(initialpos, theta, m_water, V_air, P_bottle, ...
    m_bottle, c_d, A, mu, T_atm, P_atm, c_discharge, wind);
integrate(thermo)

thermo.makeplot3d('x', 'y', 'z', {'Position', [100 100 900 700]}, {});
% thermo.makeplot('t', 'z', {'Position', [100 100 900 700]}, {});