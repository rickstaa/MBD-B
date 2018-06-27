%% MBD_B: Assignment 4 - Double pendulum systemetic approach
%  Rick Staa (4511328)
%  Last edit: 19/03/2018
% Question A: Redoo assignment 2
clear all; close all; %clc;
fprintf('--- A4_a ---\n');
fprintf('Now lets redo A2 (e,f) - In this case we have a slider constraint\n')

%% Script settings and parameters
parms.accuracy_bool   = 0;                                          % If set to 1 A\b will be performed instead of inv(A)*B this is more accurate but slower
variables             = {'x1dd' 'y1dd' 'phi1dd' 'x2dd' 'y2dd' 'phi2dd'}';

%% Parameters
% Segment 1
parms.L         = 0.55;                                             % [parms.m]
parms.w         = 0.05;                                             % [parms.m]
parms.t         = 0.004;                                            % [parms.m]
parms.p         = 1180;                                             % [kg/parms.m^3]
parms.m         = parms.p * parms.w * parms.t * parms.L;            % [kg]
parms.I         = (1/12) * parms.m * parms.L^2;                     % [kg*parms.m^2]

% World parameters
parms.g         = 9.81;                                             % [parms.m/s^2]

%% Calculate state_dp for initial states
% A2 - e
x0              = [0.5*pi 0.5*pi 0 0];
xdd.A2.e        = double(state_calc(x0,parms));
fprintf('\nThe result for A2 - e is:\n');
disp(table(variables,xdd.A2.e));

% A2 - f
w               = (60/60)*2*pi;             % Convert to rad/s
x0              = [0.5*pi 0.5*pi w 0];
xdd.A2.f        = double(state_calc(x0,parms));
fprintf('\nThe result for A2 - f is:\n');
disp(table(variables,xdd.A2.f));

%% Express COM in generalised coordinates
function xdd    = state_calc(x0,parms)
syms phi1 phi2 phi1p phi2p 

% Create generalized coordinate vectors
phi2            = pi - phi1;                                % Add extra constraint
q               = [phi1];
qp              = [phi1p];

% COM of the bodies expressed in generalised coordinates
x1              = (parms.L/2)*cos(phi1);
y1              = (parms.L/2)*sin(phi1);
x2              = x1 + (parms.L/2)*cos(phi1) + (parms.L/2) * cos(phi2);
y2              = y1 + (parms.L/2)*sin(phi1) + (parms.L/2) * sin(phi2);

% Calculate derivative of COM expressed in generalised coordinates (We need this for the energy equation)
x               = [x1;y1;phi1;x2;y2;phi2];
Jx_q            = simplify(jacobian(x,q));
xd              = Jx_q*qp;

%% Compute energies
T               = 0.5*xd.'*diag([parms.m;parms.m;parms.I;parms.m;parms.m;parms.I])*xd;           % Kinetic energy
V               = -([parms.m*parms.g 0 0 parms.m*parms.g 0 0]*x);                    % Potential energy

%% Calculate the terms of the jacobian
Q                = 0;                           % Non-conservative forces

% Partial derivatives of Kinetic energy
T_q             = simplify(jacobian(T,q));
T_qd            = simplify(jacobian(T,qp));
T_qdqd          = simplify(jacobian(T_qd,qp));
T_qdq           = simplify(jacobian(T_qd,q));

% Partial derivatives of Potential energy
V_q             = simplify(jacobian(V,q));
V_qd            = simplify(jacobian(V,qp));
V_qdqd          = simplify(jacobian(V_qd,qp));

% Make matrix vector product
M                = T_qdqd;
F                = Q + T_q' - V_q' - T_qdq*qp;

% Solve Mqdp=F to get the accelerations
if parms.accuracy_bool == 0 
    qdd          = inv(M)*F;        % Less accurate but in our case faster
else
    qdd          = M\F;            % More accurate but it is slow
end

%% Get back to COM coordinates
xdd              = simplify(jacobian(xd,qp))*qdd+simplify(jacobian(xd,q))*qp;
xdd              = subs(xdd,{phi1 phi2 phi1p phi2p},{x0(1) x0(2) x0(3) x0(4)});

end