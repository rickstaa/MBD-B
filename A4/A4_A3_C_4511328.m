%% MBD_B: Assignment 4 - Double pendulum systemetic approach
%  Rick Staa (4511328)
%  Last edit: 19/03/2018
% Question A: Redoo assignment 1
clear all; close all; clc;
fprintf('--- A4_a ---\n');
fprintf('Now lets redo A3 - In this case there is an impact constrant\n')

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
parms.omega           = 120*(2*pi/60);
parms.e               = 0;

%% Calculate the initial velocities
x0(1)                 = 0.5*pi;
x0(2)                 = 0.5*pi;
x0(3)                 = parms.omega;
x0(4)                 = parms.omega;

% Calculate the second derivative of the state
[qdd_tmp,xdd_tmp]     = state_calc(x0,parms);
qdd.A3.c              = qdd_tmp;
xdd.A3.c              = xdd_tmp;

% Calculate additional info
fprintf('\nThe result for A3 - (e = %0.1f) are:\n\n',parms.e);
disp(table(variables,xdd.A3.c));
disp(table({'phi1dd','phi2dd','lambda1','lambda2'}',qdd.A3.c));

%% Compute the derivatives of constraint matrix
function [qdd, xdd] = state_calc(x0,parms)

% Create symbolic variables
syms phi1 phi2 phi1p phi2p

% Define genaralized coordinates
q = [phi1 phi2]';
qd = [phi1p phi2p]';

% Ex[ress coordinatates of CM in generalised coordinates
x1              = (parms.L/2)*cos(phi1);
y1              = (parms.L/2)*sin(phi1);
x2              = parms.L*cos(phi1) + (parms.L/2) * cos(phi2);
y2              = parms.L*sin(phi1) + (parms.L/2) * sin(phi2);
x               = [x1;y1;phi1;x2;y2;phi2];
Jx_q            = simplify(jacobian(x,q'));
xd              = Jx_q*qd;

%% Set constraints
C               = [parms.L*cos(phi1); ...
                   parms.L*cos(phi1)+parms.L*cos(phi2)];
Cq              = simplify(jacobian(C,q'));

%% Define Matrices
M               = Jx_q.'*diag([parms.m,parms.m,parms.I,parms.m,parms.m,parms.I])*Jx_q;                                                        % Reduced generalised M matrix
M_big           = [M Cq.';Cq zeros(2,2)];                                                     % Full M matrix
F               = [M*x0(3:4).'; -parms.e*Cq*x0(3:4).'];                                       % We need the initial velocities

% Solve Mqdp=F to get the accelerations
if parms.accuracy_bool == 0 
    qdd         = inv(M_big)*F;        % Less accurate but in our case faster
else
    qdd         = M_big\F;            % More accurate but it is slow
end

%% Get back to COM coordinates
qdd             = double(subs(qdd, [phi1,phi2,phi1p,phi2p],[x0(1),x0(2),x0(3),x0(4)]));
xdd             = simplify(jacobian(x,q'))*qdd(1:2);
xdd             = double(subs(xdd, [phi1,phi2,phi1p,phi2p],[x0(1),x0(2),x0(3),x0(4)]));       % Impact velocities

end
