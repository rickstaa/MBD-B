%% MBD_B: Assignment 4 - Double pendulum systemetic approach
%  Rick Staa (4511328)
%  Last edit: 19/03/2018
% Question A: Redoo assignment 3
clear all; close all; clc;
fprintf('--- A4_a ---\n');
fprintf('Now lets redo A3 (b) - In this we have a motor constraint\n')
 
%% Script settings and parameters
variable              = {'x1dd' 'y1dd' 'phi1dd' 'x2dd' 'y2dd' 'phi2dd'}';
 
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
 
% Calculate the second derivative of the state
[dxx_tmp,Motor_torque] = state_calc(x0,parms);
dxx.A3.b = dxx_tmp;

%% Compute the derivatives of constraint matrix
function [xdp, Motor_torque] = state_calc(x0,parms)
% Create symbolic variables
syms phi1 phi2 phi1p phi2p t
 
% Define genaralized coordinates
q               = [phi1 phi2]';
qd              = [phi1p phi2p]';
 
% Express coordinatates of CM in generalised coordinates
x1              = (parms.L/2)*cos(phi1);
y1              = (parms.L/2)*sin(phi1);
x2              = x1 + (parms.L/2)*cos(phi1) + (parms.L/2) * cos(phi2);
y2              = y1 + (parms.L/2)*sin(phi1) + (parms.L/2) * sin(phi2);
x               = [x1;y1;phi1;x2;y2;phi2];

% Calculate jacobian of the state to the generalised coordinates
Jx_q            = simplify(jacobian(x,q'));
xd              = Jx_q*qd;                % First derivative of the state to the generalised coordinates
 
%% Constraints
C = phi1 - parms.omega*t;
Cq = simplify(jacobian(C,q'));
Cd = Cq*qd;                               % First derivative of the motor constraint to the generalised coordinates
Cdd = simplify(jacobian(Cd,q')*qd);       % Convective term
 
%% Calculate potential en kinetic energy
T = 0.5*xd'*diag([parms.m,parms.m,parms.I,parms.m,parms.m,parms.I])*xd;
V = -([parms.m*parms.g 0 0 parms.m*parms.g 0 0]*xd);
 
%% Building the Lagrance equations of motion
Tqd = simplify(jacobian(T,qd'))';
Tqdqd = simplify(jacobian(Tqd,qd'));
Tqdq = simplify(jacobian(Tqd,q'))*qd;     % Convective terms of the state
Tq = simplify(jacobian(T,q'))';
Vq = simplify(jacobian(V,q'))';

% Create systen of equations
M = [Tqdqd Cq';Cq 0];
F = [-Tqdq + Tq - Vq; -Cdd];
 
% Calculate acceleration
qdd = M\F;
qdd = (subs(qdd, [phi1,phi1p,phi2,phi2p],[x0(1),x0(2),x0(3),x0(4)]));
 
%% Express back in  CM coordinates
xdp = simplify(jacobian(xd,qd'))*qdd(1:2) + simplify(jacobian(xd,q'))*qd;
Motor_torque = double(qdd(3));
xdp = double(subs(xdp, [phi1,phi2,phi1p,phi2p],[x0(1),x0(2),x0(3),x0(4)]));
end
