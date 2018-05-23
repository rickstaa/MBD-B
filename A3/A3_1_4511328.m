%% MBD_B: Assignment 3
% Question 1 - Passive Element

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

% Spring parameters
parms.k                      = (15/2)*parms.m*parms.g/parms.L;     % stiffness of spring

%% Compute constraint matrices

% Use symbolic toolbox to calculate derivatives
syms x1 y1 phi1 x2 y2 phi2
syms dx1 dy1 dphi1 dx2 dy2 dphi2
x = [x1 y1 phi1 x2 y2 phi2];
dx = [dx1 dy1 dphi1 dx2 dy2 dphi2];
parms.x = x;
parms.dx = dx;

% The normal constrain equations
C = [x1-(parms.L/2)*cos(phi1); y1-(parms.L/2)*sin(phi1) ;   ...
    (x2-(parms.L/2)*cos(phi2))-(x1+(parms.L/2)*cos(phi1));              ...
    (y2-(parms.L/2)*sin(phi2))-(y1+(parms.L/2)*sin(phi1))];

% Calculate the jacobian to create the constraint equations
Cx = jacobian(C,x);
Cx = simplify(Cx);

% Take the second derivative to create the gluing constraints
Cd = Cx*dx';
Cdp = jacobian(Cd,x)*dx';
Cdp = simplify(Cdp);

% Now also calculate the spring constraint terms
Cs = sqrt((x1 + (parms.L/6)*cos(phi1) + parms.L/2)^2 + (y1 + (parms.L/6)*sin(phi1))^2) - 2*(parms.L/3);
Csx = jacobian(Cs,x); % Jacobian of Cs
Csx = simplify(Csx)';

% Put Csx Cdp Cd cx in parms struct and feed tem into the solver
syst.Cx     = Cx;
syst.Cd     = Cd;
syst.Cdp    = Cdp;
syst.Cs     = Cs;
syst.Csx    = Csx;

%% A:
x0  = [0.5*pi 0.5*pi 0 0];
[X] = state_calc(x0,parms,syst);
X   = double(vpa(X));

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

% Calculate missing initial states
x1 = 0.5*parms.L*cos(x0(1));
y1 = 0.5*parms.L*sin(x0(2));
x1 = 0; 
y1 = parms.L/2;

Cx = subs(Cx, {'phi1','phi2','dphi1','dphi2'},...
    [x0(1),x0(2),x0(3),x0(4)]);

Cdp = subs(Cdp, {'phi1','phi2','dphi1','dphi2'},...
    [x0(1),x0(2),x0(3),x0(4)]);

Cs = subs(Cs, {'phi1','phi2','x1','y1'},...
    [x0(1),x0(2),x1,y1]);

Csx = subs(Csx, {'phi1','phi2','x1','y1'},...
    [x0(1),x0(2),x1,y1]);

sigma = parms.k*Cs; % Spring Force

% Create matrices
M = diag([parms.m,parms.m,parms.I,parms.m,parms.m,parms.I]);
A = [M Cx';Cx zeros(4,4)];
F = [parms.m*parms.g 0 0 parms.m*parms.g 0 0]' - Csx*sigma; % Updated F vector
b = [F;-Cdp];
X = A\b;

end
