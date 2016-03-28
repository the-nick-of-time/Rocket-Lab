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
% objects(1).isolate()

static = statictest('./STATIC/LAtest2', 1.652e3, .2, 1);
initialpos = [0 0 .1];
theta = 45;
m_bottle = .2;
c_d = .4;
A = .008659;
mu = .2;
T_atm = 20;
P_atm = 84000;
wind = [0 0 0 0;0 0 0 1];
rocket = IspModel(static, initialpos, theta, m_bottle, c_d, A, mu,...
                T_atm, P_atm, wind);

integrate(rocket)

rocket.makeplot('t', 'z', {'Position', [100 100 900 700]}, {});
rocket.makeplot('t', 'v', {'Position', [100 100 900 700]}, {});
rocket.makeplot('y', 'z', {'Position', [100 100 900 700]}, {});
% [~, vars] = rocket.maps();