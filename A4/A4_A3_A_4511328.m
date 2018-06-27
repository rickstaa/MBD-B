%% MBD_B: Assignment 4 - Double pendulum systemetic approach
%  Rick Staa (4511328)
%  Last edit: 19/03/2018
% Question A: Redoo assignment 3
clear all; close all; clc;
fprintf('--- A4_a ---\n');
fprintf('Now lets redo A3 (a) - In this case we have a passive element\n')

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

% World parameters
parms.g               = 9.81;                                             % [parms.m/s^2]
parms.k               = (15/2)*parms.m*parms.g/parms.L;                   % stiffness of spring

%% Calculate state_dp for initial states
% A3 - a
x0                 = [0.5*pi 0.5*pi 0 0];
[qdd_tmp, xdd_tmp] = state_calc(x0,parms);
xdd.A3.a           = double(xdd_tmp);
qdd.A3.a           = double(qdd_tmp);
fprintf('\nThe result for A3 - a is:\n');
disp(table(variables,xdd.A3.a));
disp(table({'phi1dd','phi2dd'}',qdd.A3.a));

%% Express COM in generalised coordinates
function [qdd, xdd] = state_calc(x0,parms)
syms phi1 phi2 phi1p phi2p 

% Create generalized coordinate vectors
q               = [phi1; phi2];
qd              = [phi1p; phi2p];

% COM of the bodies expressed in generalised coordinates
x1              = (parms.L/2)*cos(phi1);
y1              = (parms.L/2)*sin(phi1);
x2              = parms.L*cos(phi1) + (parms.L/2) * cos(phi2);
y2              = parms.L*sin(phi1) + (parms.L/2) * sin(phi2);

% Calculate derivative of COM expressed in generalised coordinates (We need this for the energy equation)
x               = [x1;y1;phi1;x2;y2;phi2];
Jx_q            = simplify(jacobian(x,q));
xd              = Jx_q*qd;

%% Compute energies
T               = 0.5*xd.'*diag([parms.m;parms.m;parms.I;parms.m;parms.m;parms.I])*xd;          % Kinetic energy

% Now calculate spring potential energy
Cs              = sqrt((x1 + (parms.L/6)*cos(phi1) + parms.L/2)^2 + (y1 + (parms.L/6)*sin(phi1))^2) - 2*(parms.L/3);
Cs              = 0.5*parms.k*(Cs)^2;                                                           % The potential energy due to the spring

% Add to gravity potential energy
V               = -([parms.m*parms.g 0 0 parms.m*parms.g 0 0]*x)+Cs;                            % Potential energy

%% Calculate the terms of the jacobian

% Partial derivatives of Kinetic energy
T_q             = simplify(jacobian(T,q));
T_qd            = simplify(jacobian(T,qd));
T_qdqd          = simplify(jacobian(T_qd,qd));
T_qdq           = simplify(jacobian(T_qd,q));

% Partial derivatives of Potential energy
V_q             = simplify(jacobian(V,q));
V_qd            = simplify(jacobian(V,qd));
V_qdqd          = simplify(jacobian(V_qd,qd));

% Non-conservative forces
% Q                = Jx_q.'*[0 10 10*(parms.L/2)*cos(phi1) 0 0 0].';     
Q                = 0;

% Make matrix vector product
M                = T_qdqd;
F                = Q + T_q' - V_q' - T_qdq*qd;

% Solve Mqdp=F to get the accelerations
if parms.accuracy_bool == 0 
    qdd          = inv(M)*F;        % Less accurate but in our case faster
else
    qdd          = M\F;            % More accurate but it is slow
end

%% Get back to COM coordinates
qdd              = double(subs(qdd, [phi1,phi2,phi1p,phi2p],[x0(1),x0(2),x0(3),x0(4)]));
xdd              = simplify(jacobian(xd,qd))*qdd+simplify(jacobian(xd,q))*qd;
xdd              = double(subs(xdd,{phi1 phi2 phi1p phi2p},{x0(1) x0(2) x0(3) x0(4)}));

end