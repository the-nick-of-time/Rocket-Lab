clear
close all
% objects = [statictest];
% for i = 4:-1:1
%     name = sprintf('./STATIC/LAtest%d',i);
%     objects(i) = statictest(name, 1.652e3, .2, 1);
%     objects(i).makeplot()
%     objects(i).getImpulse()
% end
% objects(1).isolate()

static = statictest('./STATIC/LAtest2', 1.652e3, .2, 1);
rocket = IspModel(static, [0 0 .1], 45, .2, .4, 1.6, .2, 20, 84000, ...
    [0 0 0 0;0 0 0 1]);

integrate(rocket)

rocket.makeplot('t', 'z', {'Position', [100 100 700 500]}, {});
rocket.makeplot('t', 'v', {'Position', [100 100 700 500]}, {});
rocket.makeplot('y', 'z', {'Position', [100 100 700 500]}, {});
% [~, vars] = rocket.maps();