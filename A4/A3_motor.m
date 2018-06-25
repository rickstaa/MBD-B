%% MBD_B: Assignment 4 - Double pendulum systemetic approach
%  Rick Staa (4511328)
%  Last edit: 19/03/2018
% Question A: Redoo assignment 1

% clear all; close all; clc;
fprintf('--- A4_a ---\n');
fprintf('First lets redo A3 - In this case there are no constraints\n')

%% Script settings and parameters
variable              = {'x1dp' 'y1dp' 'phi1dp' 'x2dp' 'y2dp' 'phi2dp'}';

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
parms.omega           =-120*(2*pi/60);

%% Calculate the accelerations

x0(1) = 0.5*pi;
x0(2) = 0.5*pi;
x0(3) = parms.omega;
x0(4) = parms.omega;

[dxx,Motor_torque] = state_calc(x0,parms);

%% Compute the derivatives of constraint matrix
function [xdp, Motor_torque] = state_calc(x0,parms)
% Create symbolic variables
syms phi1 phi2 phi1p phi2p t

% Define genaralized coordinates
q = [phi1 phi2]';
qp = [phi1p phi2p]';

% Ex[ress coordinatates of CM in generalised coordinates
x1              = (parms.L/2)*cos(phi1);
y1              = (parms.L/2)*sin(phi1);
x2              = x1 + (parms.L/2)*cos(phi1) + (parms.L/2) * cos(phi2);
y2              = y1 + (parms.L/2)*sin(phi1) + (parms.L/2) * sin(phi2);
x               = [x1;y1;phi1;x2;y2;phi2];
Xj_q = simplify(jacobian(x,q'));
xp = Xj_q*qp;

%% Constraints
C = phi1 - parms.omega*t;
Cq = simplify(jacobian(C,q'));
Cp = Cq*qp;
Cdp = simplify(jacobian(Cp,q')*qp);

%% Calculate potential en kinetic energy
T = 0.5*xp'*diag([parms.m,parms.m,parms.I,parms.m,parms.m,parms.I])*xp;
V = -([parms.m*parms.g 0 0 parms.m*parms.g 0 0]*xp);

%% Building the Lagrance equations of motion
Tqp = simplify(jacobian(T,qp'))';
Tqpqp = simplify(jacobian(Tqp,qp'));
Tqpq = simplify(jacobian(Tqp,q'))*qp;
Tq = simplify(jacobian(T,q'))';
Vq = simplify(jacobian(V,q'))';
M = [Tqpqp Cq';Cq 0];
F = [-Tqpq + Tq - Vq; -Cdp];

% Calculate acceleration
qdp = inv(M)*(F);
qdp = (subs(qdp, [phi1,phi1p,phi2,phi2p],[x0(1),x0(2),x0(3),x0(4)]));

%% Express back in  CM coordinates
xdp = simplify(jacobian(xp,qp'))*qdp(1:2) + simplify(jacobian(xp,q'))*qp;
Motor_torque = double(qdp(3));
xdp = double(subs(xdp, [phi1,phi2,phi1p,phi2p],[x0(1),x0(2),x0(3),x0(4)]));
end