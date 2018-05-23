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

% World parameters
parms.g               = 9.81;                                             % [parms.m/s^2]
parms.k                      = (15/2)*parms.m*parms.g/parms.L;     % stiffness of spring

%% Calculate state_dp for initial states
% A2 - a:
% A1 - g
x0              = [0.5*pi 0.5*pi 0 0];
xdp_A1_e        = double(state_calc(x0,parms));
fprintf('\nThe result for A3 - a is:\n');
disp(table(variable,xdp_A1_e));


%% Express COM in generalised coordinates
function xdp = state_calc(x0,parms)
syms phi1 phi2 phi1p phi2p 

% Create generalized coordinate vectors
q               = [phi1; phi2];
qp              = [phi1p; phi2p];

% COM of the bodies expressed in generalised coordinates
x1              = (parms.L/2)*cos(phi1);
y1              = (parms.L/2)*sin(phi1);
x2              = x1 + (parms.L/2)*cos(phi1) + (parms.L/2) * cos(phi2);
y2              = y1 + (parms.L/2)*sin(phi1) + (parms.L/2) * sin(phi2);

% Calculate derivative of COM expressed in generalised coordinates (We need this for the energy equation)
x               = [x1;y1;phi1;x2;y2;phi2];
Jx_q            = simplify(jacobian(x,q));
xp              = Jx_q*qp;

%% Compute energies
T               = 0.5*xp.'*diag([parms.m;parms.m;parms.I;parms.m;parms.m;parms.I])*xp;           % Kinetic energy

% Now calculate spring potential energy
Ls = sqrt((x1 + (parms.L/6)*cos(phi1) + parms.L/2)^2 + (y1 + (parms.L/6)*sin(phi1))^2) - 2*(parms.L/3);
Vs = 0.5*parms.k*(Ls)^2; % The potential energy due to the spring

% Add to gravity potential energy
V               = -([parms.m*parms.g 0 0 parms.m*parms.g 0 0]*x)+Vs;                                % Potential energy

%% Calculate the terms of the jacobian

% Partial derivatives of Kinetic energy
T_q             = simplify(jacobian(T,q));
T_qp            = simplify(jacobian(T,qp));
T_qpqp          = simplify(jacobian(T_qp,qp));
T_qpq           = simplify(jacobian(T_qp,q));

% Partial derivatives of Potential energy
V_q             = simplify(jacobian(V,q));
V_qp            = simplify(jacobian(V,qp));
V_qpqp          = simplify(jacobian(V_qp,qp));

% Non-conservative forces
Q                = Jx_q.'*[0 10 10*(parms.L/2)*cos(phi1) 0 0 0].';     

% Make matrix vector product
M                = T_qpqp;
F                = Q + T_q' - V_q' - T_qpq*qp;

% Solve Mqdp=F to get the accelerations
qdp              = inv(M)*F;

%% Get back to COM coordinates
xdp              = simplify(jacobian(xp,qp))*qdp+simplify(jacobian(xp,q))*qp;
xdp              = subs(xdp,{phi1 phi2 phi1p phi2p},{x0(1) x0(2) x0(3) x0(4)});

end