%% MBD_B: Assignment 3
% Question 3 - impact

%  Rick Staa (4511328)
%  Last edit: 05/03/2018
clear all; close all; % clc;

%% Parameters
% Segment 1
parms.L     = 0.55;                                             % [m]
parms.w     = 0.05;                                             % [m]
parms.t     = 0.004;                                            % [m]
parms.p     = 1180;                                             % [kg/m^3]
parms.m     = parms.p * parms.w * parms.t * parms.L;            % [kg]
parms.I     = (1/12) * parms.m * parms.L^2;                     % [kg*m^2]

% World parameters
parms.g      = 9.81;                     % [m/s^2]

% Motor constraint
omega = 120*(2*pi/60);
parms.e = 1; % Coefficient of restitution

%% Compute constraint matrices

% Use symbolic toolbox to calculate derivatives
syms x1 y1 phi1 x2 y2 phi2 t
syms dx1 dy1 dphi1 dx2 dy2 dphi2
x = [x1 y1 phi1 x2 y2 phi2];
dx = [dx1 dy1 dphi1 dx2 dy2 dphi2];
parms.x = x;
parms.dx = dx;

% The normal constrain equations
C = [x1-(parms.L/2)*cos(phi1); y1-(parms.L/2)*sin(phi1) ;   ...
    (x2-(parms.L/2)*cos(phi2))-(x1+(parms.L/2)*cos(phi1));  ...
    (y2-(parms.L/2)*sin(phi2))-(y1+(parms.L/2)*sin(phi1));  ...
    (x1 + (parms.L/2)*cos(phi1));                                 ...
    (x2 + (parms.L/2)*cos(phi2))]; % Extra constraint 

% Calculate the jacobian to create the constraint equations
Cx = jacobian(C,x);
Cx = simplify(Cx);

% Put Csx Cdp Cd cx in parms struct and feed tem into the solver
syst.Cx     = Cx;

%% C:
parms.phi1_0  = pi/2;   % angle of bar 1 (with x-axis)
parms.phi2_0  = pi/2;   % angle of bar 2 
dphi1_0 = 120*(2*pi/60);  % angular velocity of bar 1
dphi2_0 = 120*(2*pi/60);  % angular velocity of bar 2
dx1_0 = -(parms.L/2)*sin(parms.phi1_0)*dphi1_0;
dx2_0 = -(3*parms.L/2)*sin(parms.phi2_0)*dphi2_0;
dy1_0 = 0; 
dy2_0 = 0;

% Calculate impatct forces
x0      = [dx1_0;dy1_0;dphi1_0;dx2_0;dy2_0;dphi2_0];
[X]     = state_calc(x0,parms,syst);
X       = double(vpa(X));
M       = diag([parms.m,parms.m,parms.I,parms.m,parms.m,parms.I]);
K_be    = 0.5*x0'*M*x0;
K_af    = 0.5*X(1:6)'*M*X(1:6);
Work    = 0.5*X(7:end)'*x0*(1-parms.e);
display(X)
display(K_be)
display(K_af)
display(Work)

%% A: Create state space matrices for the case when we add a spring
function [X] = state_calc(x0,parms,syst)

% Get system matrices out
structname_fields = fields(parms);
for i = 1:size(fields(parms))
    eval_str = [structname_fields{i,:},'=','parms.',structname_fields{i,:},';'];
    eval(eval_str);
end
structname_fields = fields(syst);
for i = 1:size(fields(syst))
    eval_str = [structname_fields{i,:},'=','syst.',structname_fields{i,:},';'];
    eval(eval_str);
end

Cx = subs(Cx, {'phi1','phi2','dphi1','dphi2'},...
    [parms.phi1_0,parms.phi2_0,x0(3),x0(6)]);

% Create matrices

M = diag([parms.m,parms.m,parms.I,parms.m,parms.m,parms.I]);
A = [M Cx';Cx zeros(6,6)];
F = [M*x0];
b = [F;-parms.e*Cx*x0];
X = A\b;

end
