%% MBD_B: Assignment 4 - Double pendulum systemetic approach
%  Rick Staa (4511328)
%  Last edit: 19/03/2018
% Question A: Redoo assignment 1
clear all; close all; clc;
fprintf('--- A4_a ---\n');
fprintf('First lets redo A3 - In this case there are no constraints\n')

%% Script settings and parameters
parms.accuracy_bool   = 0;                                          % If set to 1 A\b will be performed instead of inv(A)*B this is more accurate but slower
variables             = {'x1dd' 'y1dd' 'phi1dd' 'x2dd' 'y2dd' 'phi2dd'}';

%% Parameters
% Segment 1
parms.L               = 0.55;                                             % [parms.m]
parms.w               = 0.05;                                             % [parms.m]
parms.t               = 0.004;                                            % [parms.m]
parms.p               = 1180;                                             % [kg/parms.m^3]
parms.m               = parms.p * parms.w * parms.t * parms.L;            % [kg]
parms.I               = (1/12) * parms.m * parms.L^2;                     % [kg*parms.m^2]
parms.omega           = -120*(2*pi/60);

% World parameters
parms.g               = 9.81;                                             % [parms.m/s^2]
parms.omega           = -120*(2*pi/60);

%% Calculate the accelerations

x0(1) = 0.5*pi;
x0(2) = 0.5*pi;
x0(3) = parms.omega;
x0(4) = parms.omega;

[qdd_tmp, xdd_tmp, Motor_torque] = state_calc(x0,parms);
xdd.A3.b               = double(xdd_tmp);
qdd.A3.b               = double(qdd_tmp);
disp(table(variables,xdd.A3.b));

% Calculate additional info
fprintf('\nThe result for A3 - a is:\n');
disp(table(variables,xdd.A3.b));
disp(table({'phi1dd','phi2dd','T_motor'}',qdd.A3.b));

%% Compute the derivatives of constraint matrix
function [qdd, xdd, Motor_torque] = state_calc(x0,parms)
% Create symbolic variables
syms phi1 phi2 phi1d phi2d t

% Define genaralized coordinates
q = [phi1 phi2]';
qd = [phi1d phi2d]';

% Ex[ress coordinatates of CM in generalised coordinates
x1              = (parms.L/2)*cos(phi1);
y1              = (parms.L/2)*sin(phi1);
x2              = parms.L*cos(phi1) + (parms.L/2) * cos(phi2);
y2              = parms.L*sin(phi1) + (parms.L/2) * sin(phi2);
x               = [x1;y1;phi1;x2;y2;phi2];
Jx_q            = simplify(jacobian(x,q'));
xd              = Jx_q*qd;

%% Constraints
C = phi1 - parms.omega*t;
Cq = simplify(jacobian(C,q'));
Cd = Cq*qd;
Cdd = simplify(jacobian(Cd,q')*qd);

%% Calculate potential en kinetic energy
T = 0.5*xd'*diag([parms.m,parms.m,parms.I,parms.m,parms.m,parms.I])*xd;
V = -([parms.m*parms.g 0 0 parms.m*parms.g 0 0]*x);

%% Building the Lagrance equations of motion
% See page 50
Tqd             = simplify(jacobian(T,qd'))';
Tqdqd           = simplify(jacobian(Tqd,qd'));
Tqdq            = simplify(jacobian(Tqd,q'))*qd;
Tq              = simplify(jacobian(T,q'))';
Vq              = simplify(jacobian(V,q'))';
M               = [Tqdqd Cq';Cq 0];
F               = [-Tqdq + Tq - Vq; -Cdd];

% Solve Mqdp=F to get the accelerations
if parms.accuracy_bool == 0 
    qdd          = inv(M)*F;       % Less accurate but in our case faster
else
    qdd          = M\F;            % More accurate but it is slow
end

%% Express back in  CM coordinates
qdd = (subs(qdd, [phi1,phi2,phi1d,phi2d],[x0(1),x0(2),x0(3),x0(4)]));
xdd = simplify(jacobian(xd,qd'))*qdd(1:2) + simplify(jacobian(xd,q'))*qd;
Motor_torque = double(qdd(3));
xdd = double(subs(xdd, [phi1,phi2,phi1d,phi2d],[x0(1),x0(2),x0(3),x0(4)]));
end