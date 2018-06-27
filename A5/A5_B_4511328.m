%% MBD_B: Assignment 5 - Passive dynamic walker (Lagrange and TMT aproach)
%  Rick Staa (4511328)
%  Last edit: 27/03/2018
% - Question A: Lagrange method -
clear all; close all; %clc;
fprintf('--- A5_b ---\n');
tic

%% Script settings and parameters
parms.accuracy_bool   = 0;                                  % If set to 1 A\b will be performed instead of inv(A)*B this is more accurate but slower

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
Jx_q    = simplify(jacobian(x,q.')); 

% Calculate derivative of the state vector
xp      = Jx_q*qd;

% Solve with virtual power
M_bar   = Jx_q.'*diag([m,m,I,M,M,0,m,m,I])*Jx_q;
F       = [m*g*sin(gamma), -m*g*cos(gamma), 0, M*g*sin(gamma), -M*g*cos(gamma), 0, m*g*sin(gamma), -m*g*cos(gamma), 0];
T_qq    = simplify(jacobian(Jx_q*qd,q)*qd);
Q       = simplify(Jx_q.'*(F.' - diag([m,m,I,M,M,0,m,m,I])*T_qq));

% Calculate result expressed in generalized coordinates
if parms.accuracy_bool == 0
    qdd     = inv(M_bar)*Q;     % Less accurate but in our case faster
else
    qdd     = M_bar\Q;          % More accurate but it is slow
end

% Get result back in COM coordinates
xdd     = simplify(jacobian(xp,qd.'))*qdd + simplify(jacobian(xp,q.'))*qd;

%% Calculate for a initial state
x0    = [deg2rad(30),deg2rad(45),-pi,0];
parms = [0.8,12,9.81,36,pi/12];                            % parms = [l,m,g,M,gamma];
xdd = subs(xdd, [l,m,g,M,gamma],parms);
xdd = double(subs(xdd, [alpha,beta,alphap,betap],[x0]));
disp(xdd)
toc