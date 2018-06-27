%% MBD_B: Assignment 5 - Passive dynamic walker (Lagrange and TMT aproach)
%  Rick Staa (4511328)
%  Last edit: 27/03/2018
% - Question A: Lagrange method -
clear all; close all; %clc;
fprintf('--- A5_a ---\n');
tic 

%% Script settings and parameters
parms.accuracy_bool   = 0;                                  % If set to 1 A\b will be performed instead of inv(A)*B this is more accurate but slower
e = 0;

%% Parameters and variables
% Initialise parameters and variables
syms l m I M                                                % Initialize model Parameters
syms x1 y1 phi1 x2 y2 phi2 x3 y3 phi3                       % Initialize variables
syms x1p y1p phi1p x2p y2p phi2p x3p y3p phi3p              % Initialize derivatives
syms alpha beta alphap betap                                % Initialize generalised coordinates
syms gamma g                                                % Initialize enviroment variables

% Calculate additional parameters
I       = (1/12)*m*(l^2);                                   % Mass moment of inertia of each leg about the COM

% Devine state and its derivative
q       = [alpha;beta];
qd      = [alphap;betap];

%% Generate equations of motion

% Express coordinates of the COM in terms of generalised coordinates
x1      = 0.5*l*sin(alpha);
y1      = 0.5*l*cos(alpha);
phi1    = 0.5*pi - alpha;
x2      = l*sin(alpha);
y2      = l*cos(alpha);
phi2    = 0;
x3      = l*sin(alpha) - 0.5*l*sin(beta);
y3      = l*cos(alpha) - 0.5*l*cos(beta);
phi3    = -(0.5*pi + beta);

% Put in one state vector
x       = [x1;y1;phi1;x2;y2;phi2;x3;y3;phi3];

% Compute the jacobian of these expressions
Jx_q    = simplify(jacobian(x,q)); 

% Calculate derivative of the state vector
xd      = Jx_q*qd;

% Compute energies
T       = 0.5*xd.'*diag([m;m;I;M;M;0;m;m;I])*xd;          % Kinetic energy
V       = -([m*g*sin(gamma) -m*g*cos(gamma) 0 M*g*sin(gamma) -M*g*cos(gamma) 0 m*g*sin(gamma) -m*g*cos(gamma) 0]*x);                    % Potential energy

%% Now compute the lagrangian constraint equations of motion

% Compute partial derivatives w.r.t q
T_q     = simplify(jacobian(T,q.')).';                          % Term 2 in Lagrange equation reader
V_q     = simplify(jacobian(V,q.')).';                          % Term 3 ...

% Compute partial derivatives w.r.t. qp
T_qd    = simplify(jacobian(T,qd.')).';

% Compute terms of total derivative
T_qd_qd  = simplify(jacobian(T_qd,qd.'));
T_qd_q   = simplify(jacobian(T_qd,q.'));

% Impact constraint
C       = l*cos(alpha) - l*cos(beta);
Jc_q    = simplify(jacobian(C,q.'));

% Calculate M_big
M_big = Jx_q.'*diag([m;m;I;M;M;0;m;m;I])*Jx_q;

% Combine matrixes to get the lagrangian euqations of motion in matrix
% vector form
Qj       = 0;                                                  % Non-conservative forces
Q        = [M_big Jc_q.';Jc_q zeros(1,1)];
F        = [M_big*qd;-e*Jc_q*qd];

% Calculate result expressed in generalized coordinates
if parms.accuracy_bool == 0
    qdd     = inv(Q)*F;     % Less accurate but in our case faster
else
    qdd     = Q\F;          % More accurate but it is slow
end

% Get result back in COM coordinates
xdd     = simplify(jacobian(xd,qd.'))*qdd(1:2) + simplify(jacobian(xd,q.'))*qd;

%% Calculate velocites of COM after impact
x0    = [pi/6,-pi/6,1.2,0.2];
parms = [0.75,3,9.81,5,0];                            % parms = [l,m,g,M,gamma];
xdd   = subs(xdd, [l,m,g,M,gamma],parms);
xdd   = double(subs(xdd, [alpha,beta,alphap,betap],[x0]));
disp(xdd)
toc 

%% Calculate velocities in C
qdd_value = double(subs(qdd, [l,m,g,M,gamma,alpha,beta,alphap,betap],[parms,x0]));% Calculate numeric value of qdd = qd_after
C_x       = l*sin(alpha)-l*sin(beta);
Jx_q      = jacobian(C_x,q);
x_cd      = Jx_q*qd;
x_cd      = double(subs(x_cd,[l,m,g,M,gamma,alpha beta alphap betap],[parms,pi/6 -pi/6  qdd_value(1) qdd_value(2)]));