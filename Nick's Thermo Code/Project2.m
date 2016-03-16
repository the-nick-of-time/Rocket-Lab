% ASEN 2012 Project 2
% Bottle Rocket Dynamics
% Author: 30a5bf672aff
% Written 2015-11-09 
% Last modified 2015-12-06
% This code models the motion of a pressure-based bottle rocket, solving
% the differential equations involved using the function ode45.
% The results of these integrations are then graphed.

% Known information is set out in the sections "Constants" and "Initial
% Conditions" of this script.

% Assumptions made:
% 1. The air expands within the bottle isentropically.
% 2. The water and the air both exit the bottle isentropically.
% Of these, the second is more likely to cause error. From viscous flow
% theory, we know that flow will be slowed near a solid surface. This means
% that, near the outside of the bottle mouth, the velocity of the fluid
% will not be quite as high as expected due to these frictional losses.

% This code performs its main purpose, calculating the values for the given
% equation, in a single call to ode45. This is possible as all of the
% functions are collected (by means of function handles) into a single
% function. Also, each function reacts to the current stage of the problem
% and returns the correct values (rather than needing chained
% integrations).



% p_atm was retrived from the NCAR website, 
% http://www.eol.ucar.edu/cgi-bin/weather.cgi?site=fl&units=si
% T_atm retrieved from NOAA using February data
% http://www.esrl.noaa.gov/psd/boulder/dailyrecords/

%%
clear all
close all
clc
%% Conversion Equations
% Degrees to radians
deg2rad=@(deg) deg*pi/180;
psi2pa=@(psi) psi*6894.76;

%% Constants
drag_coefficient=0.5*.9;
gravity=9.81; %m/s^2
cross_section=.008659; %m^2
rho_water=1000; %kg/m^3
gamma=1.4;
discharge_coefficient=0.8;
Vbottle=.002; %.002048; %m^3
R_air=287; %Pa*m^3 / kg*K
T_atm=280.4; %K
p_atm=psi2pa(12.03); %84389; %Pa
m_bottle=0.07; %kg
area_mouth=pi*.0105^2; %m^2 

%% Global Variables
global stage setstage;
stage=1;
%% Initial Conditions
v_initial=0;
x_initial=0;
y_initial=0.1;
theta_initial=deg2rad(45);
V_air_initial=Vbottle*0.5*1.1;

p_initial=(psi2pa(50)+p_atm);
T_initial=300;

m_air_initial= V_air_initial*(p_initial/(R_air*T_atm));
m_initial= m_bottle + rho_water*(Vbottle-V_air_initial) + m_air_initial;

%% Constant Calculated Values
% rho_atm=p_atm/(R_air*T_atm);
rho_atm=.961;
p_end=p_initial * (V_air_initial/Vbottle)^gamma;
T_end=T_initial * (V_air_initial/Vbottle)^gamma;

%% Parametrization Functions
% V=      velocity magnitude (speed)
% x=      horizontal distance
% y=      vertical distance
% theta=  angle from horizontal
% m=      mass of the whole rocket
% vol=    volume of air inside the rocket
pressureFn=@(t,vars) pressure(t,vars,p_initial,p_end,p_atm,V_air_initial,...
                              m_air_initial,m_bottle,gamma);

setstage=@(t,vars) findStage(t,vars,pressureFn,Vbottle,p_atm);

Thetadot=@(t,vars) thetadot(vars(4),vars(1),gravity);

Volumedot=@(t,vars) Voldot(p_initial,p_atm,V_air_initial,vars(6),Vbottle,...
                    gamma,discharge_coefficient,area_mouth,rho_water);
Massdot=@(t,vars) mdot(t,vars,pressureFn,vars(5)-m_bottle,...
                       discharge_coefficient,rho_water,...
                       p_atm,area_mouth,Vbottle,gamma,R_air);

Thrust=@(t,vars) thrust(t,vars,Vbottle,discharge_coefficient,pressureFn,...
                 p_atm,area_mouth,gravity,rho_water,Massdot,...
                 vars(5)-m_bottle,R_air,gamma);

Drag=@(t,vars) drag(vars(1),rho_atm,drag_coefficient,cross_section);
gravTangent=@(t,vars) -gravity*sin(vars(4));

Velocitydot=@(t,vars) veldot(t,vars,Thrust,Drag,gravTangent);

Derivatives=@(t,vars) derivs(t,vars,Velocitydot,Massdot,Thetadot,Volumedot);


%% Perform Integration
opts=odeset('Events',@endcondition,'RelTol',1e-9,'InitialStep',0.005);
trange=[0,20]; %nice and large to ensure that the end is captured
initialconditions=[v_initial;x_initial;y_initial;theta_initial;m_initial;...
    V_air_initial];
[t,vars,te,ye,ie]=ode45(Derivatives,trange,initialconditions,opts);


%% Display results
plot(vars(:,2),vars(:,3))
title('Trajectory','fontsize',16)
daspect([1 1 1])
xlabel('Horizontal Position (m)')
ylabel('VerticalPosition (m)')

%%
plot(vars(:,2),[vars(:,3) vars(:,1).*cos(vars(:,4)) vars(:,1).*sin(vars(:,4))],'-x')
h=legend('Position (m)','Horizontal velocity ($\frac{m}{s}$)',...
    'Vertical velocity ($\frac{m}{s}$)');
set(h,'Interpreter','LaTeX')
title('Trajectory Variables','fontsize',20)
xlabel('X position','fontsize',14)
daspect([1 1 1])
%%
inds=find(t < .3);
plot(t(inds),vars(inds,5)-m_bottle*ones(length(inds),1))
title('Mass over Time','fontsize',20)
xlabel('Time (s)','fontsize',14)
ylabel('Mass of Air in Rocket (kg)','fontsize',14)

%%
plot(t,vars(:,5))
title('Mass over Time','fontsize',20)
xlabel('Time (s)','fontsize',14)
ylabel('Mass of Rocket (kg)','fontsize',14)

%%
distance=vars(end,2)
%% Retrieve specific impulse
I=zeros(size(vars,1),1);
for i=1:size(vars,1)
    [~,I(i)]=Thrust(t(i),vars(i,:));
end
%% Variance of Parameters
x_initial_varied=[-.1 .1];
y_initial_varied=[1e-5 .2];
theta_initial_varied=[theta_initial*.9 theta_initial*1.1];
V_air_initial_varied=[V_air_initial*.9 V_air_initial*1.1];

drag_coefficient_varied=[.3 .5];
m_bottle_varied=[m_bottle*.9 m_bottle*1.1];

p_initial_varied=[p_initial*.9 p_initial*1.1];
%% Using varied parameters
opts=odeset('Events',@endcondition,'RelTol',1e-6,'InitialStep',0.005);
trange=[0 5];
%%
stage=1;
initialconditions_xlow=initialconditions;
initialconditions_xlow(2)=x_initial_varied(1);
[txlow,varsxlow]=ode45(Derivatives,trange,initialconditions_xlow,opts);

%%
stage=1;
initialconditions_xhigh=initialconditions;
initialconditions_xhigh(2)=x_initial_varied(2);
[txhigh,varsxhigh]=ode45(Derivatives,trange,initialconditions_xhigh,opts);

stage=1;
initialconditions_ylow=initialconditions;
initialconditions_ylow(3)=y_initial_varied(1);
[tylow,varsylow]=ode45(Derivatives,trange,initialconditions_ylow,opts);

stage=1;
initialconditions_yhigh=initialconditions;
initialconditions_yhigh(3)=y_initial_varied(2);
[tyhigh,varsyhigh]=ode45(Derivatives,trange,initialconditions_yhigh,opts);

stage=1;
initialconditions_thlow=initialconditions;
initialconditions_thlow(4)=theta_initial_varied(1);
[tthlow,varsthlow]=ode45(Derivatives,trange,initialconditions_thlow,opts);

stage=1;
initialconditions_thhigh=initialconditions;
initialconditions_thhigh(4)=theta_initial_varied(2);
[tthhigh,varsthhigh]=ode45(Derivatives,trange,initialconditions_thhigh,opts);

stage=1;
initialconditions_mlow=initialconditions;
m_initial_low= m_bottle_varied(1) + rho_water*(Vbottle-V_air_initial) + m_air_initial;
initialconditions_mlow(5)=m_initial_low;
[tmlow,varsmlow]=ode45(Derivatives,trange,initialconditions_mlow,opts);

stage=1;
initialconditions_mhigh=initialconditions;
m_initial_high= m_bottle_varied(2) + rho_water*(Vbottle-V_air_initial) + m_air_initial;
initialconditions_mhigh(5)=m_initial_high;
[tmhigh,varsmhigh]=ode45(Derivatives,trange,initialconditions_mhigh,opts);

stage=1;
initialconditions_Vlow=initialconditions;
initialconditions_Vlow(6)=V_air_initial_varied(1);
[tVlow,varsVlow]=ode45(Derivatives,trange,initialconditions_Vlow,opts);

stage=1;
initialconditions_Vhigh=initialconditions;
initialconditions_Vhigh(6)=V_air_initial_varied(2);
[tVhigh,varsVhigh]=ode45(Derivatives,trange,initialconditions_Vhigh,opts);

%%
stage=1;
initialconditions_Dlow=initialconditions;
Drag=@(t,vars) drag(vars(1),rho_atm,drag_coefficient_varied(1),cross_section);
Velocitydot=@(t,vars) veldot(t,vars,Thrust,Drag,gravTangent);
Derivatives=@(t,vars) derivs(t,vars,Velocitydot,Massdot,Thetadot,Volumedot);
[tDlow,varsDlow]=ode45(Derivatives,trange,initialconditions_Dlow,opts);

stage=1;
initialconditions_Dhigh=initialconditions;
Drag=@(t,vars) drag(vars(1),rho_atm,drag_coefficient_varied(2),cross_section);
Velocitydot=@(t,vars) veldot(t,vars,Thrust,Drag,gravTangent);
Derivatives=@(t,vars) derivs(t,vars,Velocitydot,Massdot,Thetadot,Volumedot);
[tDhigh,varsDhigh]=ode45(Derivatives,trange,initialconditions_Dhigh,opts);

stage=1;
initialconditions_plow=initialconditions;
pressureFn=@(t,vars) pressure(t,vars,p_initial_varied(1),p_end,p_atm,V_air_initial,...
                              m_air_initial,m_bottle,gamma);
setstage=@(t,vars) findStage(t,vars,pressureFn,Vbottle,p_atm);
Massdot=@(t,vars) mdot(t,vars,pressureFn,vars(5)-m_bottle,...
                       discharge_coefficient,rho_water,...
                       p_atm,area_mouth,Vbottle,gamma,R_air);
Thrust=@(t,vars) thrust(t,vars,Vbottle,discharge_coefficient,pressureFn,...
                 p_atm,area_mouth,gravity,rho_water,Massdot,...
                 vars(5)-m_bottle,R_air,gamma);
Velocitydot=@(t,vars) veldot(t,vars,Thrust,Drag,gravTangent);
Derivatives=@(t,vars) derivs(t,vars,Velocitydot,Massdot,Thetadot,Volumedot);
[tplow,varsplow]=ode45(Derivatives,trange,initialconditions_plow,opts);

stage=1;
initialconditions_phigh=initialconditions;
pressureFn=@(t,vars) pressure(t,vars,p_initial_varied(2),p_end,p_atm,V_air_initial,...
                              m_air_initial,m_bottle,gamma);
setstage=@(t,vars) findStage(t,vars,pressureFn,Vbottle,p_atm);
Massdot=@(t,vars) mdot(t,vars,pressureFn,vars(5)-m_bottle,...
                       discharge_coefficient,rho_water,...
                       p_atm,area_mouth,Vbottle,gamma,R_air);
Thrust=@(t,vars) thrust(t,vars,Vbottle,discharge_coefficient,pressureFn,...
                 p_atm,area_mouth,gravity,rho_water,Massdot,...
                 vars(5)-m_bottle,R_air,gamma);
Velocitydot=@(t,vars) veldot(t,vars,Thrust,Drag,gravTangent);
Derivatives=@(t,vars) derivs(t,vars,Velocitydot,Massdot,Thetadot,Volumedot);
[tphigh,varsphigh]=ode45(Derivatives,trange,initialconditions_phigh,opts);

%%
distances=[varsxlow(end,2);varsxhigh(end,2);varsylow(end,2);varsyhigh(end,2);
    varsthlow(end,2);varsthhigh(end,2);varsmlow(end,2);varsmhigh(end,2);
    varsVlow(end,2);varsVhigh(end,2);varsDlow(end,2);varsDhigh(end,2);
    varsplow(end,2);varsphigh(end,2)];
distances=[[txlow(end);txhigh(end);tylow(end);tyhigh(end);tthlow(end);
    tthhigh(end);tmlow(end);tmhigh(end);tVlow(end);tVhigh(end);tDlow(end);
    tDhigh(end);tplow(end);tphigh(end)] distances]