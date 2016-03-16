clear
close all
% objects = [statictest];
for i = 4:-1:1
    name = sprintf('./STATIC/LAtest%d',i);
    objects(i) = statictest(name, 1.652e3, .2, 1);
    objects(i).makeplot()
    objects(i).getImpulse()
end
% objects(1).isolate()