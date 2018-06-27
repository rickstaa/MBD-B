%% MBD_B: Assignment 1 - Double pendulum EOM using FBD and newton-euler
%  Rick Staa (4511328)
%  Last edit: 24/02/2018
clear all; close all; tic; %clc;
fprintf('--- A1 ---\n');

%% Script parameters
parms.accuracy_bool = 1;                                     % If set to 1 A\b will be performed instead of inv(A)*B this is more accurate but slower

%% Parameters
% Segment 1
parms.L1     = 0.75;                                          % [m]
% parms.w1     = 0.05;                                        % [m]
% parms.t1     = 0.004;                                       % [m]
% parms.p1     = 1180;                                        % [kg/m^3]
% parms.m1     = parms.p1 * parms.w1 * parms.t1 * parms.L1;   % [kg]
parms.m1     = 3; %[kg]
parms.I1     = (1/12) * parms.m1 * parms.L1^2;              % [kg*m^2]

% Segment 2
parms.L2     = 0.75;                                        % [m]
% parms.w2     = 0.05;                                        % [m]
% parms.t2     = 0.004;                                       % [m]
% parms.p2     = 1180;                                        % [m]
% parms.m2     = parms.p2 * parms.w2 * parms.t2 * parms.L2;   % [kg]
parms.m2     = 3; %[kg]
parms.I2     = (1/12) * parms.m2 * parms.L2^2;              % [kg*m^2]

% World parameters
parms.g      = -9.81;                                        % [m/s^2] 

%% Initial states

% b).
x0           = [0 0.5*pi 0 0];                         % [phi_1 phi_2 phi_1_p phi_2_p];
x_dd.b       = state_calc(x0,parms);

% % c).
% x0           = [0 0 0 0];                                   % [phi_1 phi_2 phi_1_p phi_2_p];
% x_dd.c       = state_calc(x0,parms);
% 
% % d).
% w            = (60/60)*2*pi;                                % Convert to rad/s
% x0           = [0 0 w w];
% x_dd.d       = state_calc(x0,parms);

% Calculate velocities
x1_d         = -(parms.L1/2)*sin(x0(1))*x0(3);
y1_d         =  (parms.L1/2)*cos(x0(1))*x0(3);
x2_d         = x1_d - (parms.L1/2)*sin(x0(1))*x0(3) - (parms.L2/2)*sin(x0(2))*x0(4);
y2_d         = y1_d + (parms.L1/2)*cos(x0(1))*x0(3) + (parms.L2/2)*cos(x0(2))*x0(4);

% % Put in table
% results      = [x_dd.b x_dd.c x_dd.d];
% x_d.d        = [x1_d y1_d x2_d y2_d];
% toc;

%% Functions

function [x_dd] = state_calc(x0,parms)
% Equations of motions + Contraints in vector matrix form

% Get variables out of struct and initial state
phi_1    = x0(1);
phi_2    = x0(2);
phi_1_p  = x0(3);
phi_2_p  = x0(4);

names   = fieldnames(parms);
for i=1:length(names)
evalc([names{i} '=parms.' names{i} ]);
end

% Create matrices
M           = diag([m1 m1 I1 m2 m2 I2]);
Fg          = [0 m1*g 0 0 m2*g 0]';
a           = [ -(L1/2)*cos(phi_1)*phi_1_p^2; ...
                -(L1/2)*sin(phi_1)*phi_1_p^2; ...
                -(L1/2)*cos(phi_1)*phi_1_p^2 - (L2/2)*cos(phi_2)*phi_2_p^2; ...
                -(L1/2)*sin(phi_1)*phi_1_p^2 - (L2/2)*sin(phi_2)*phi_2_p^2];
b           = [Fg;a];
A           = [         -1                  0                  1                  0;          ...
                         0                 -1                  0                  1;          ...
                -(L1/2)*sin(phi_1)  (L1/2)*cos(phi_1)  -(L1/2)*sin(phi_1)  (L1/2)*cos(phi_1);  ...
                         0                  0                 -1                  0;          ...
                         0                  0                  0                 -1;          ...
                         0                  0          -(L2/2)*sin(phi_2)  (L2/2)*cos(phi_2)];
                     
B           = [          1                  0           (L1/2)*sin(phi_1)         0         0           0;          ...
                         0                  1          -(L1/2)*cos(phi_1)         0         0           0;          ...
                        -1                  0           (L1/2)*sin(phi_1)         1         0   (L2/2)*sin(phi_2);  ...
                         0                 -1          -(L1/2)*cos(phi_1)         0         1  -(L2/2)*cos(phi_2)];

% Calculate state out of Ax=b and the initial state
if parms.accuracy_bool == 0
    x_dd = inv([M A;B zeros(size(B,1),size(A,2))])*b;       % Less accurate but in our case faster
else
    x_dd = [M A;B zeros(size(B,1),size(A,2))]\b;            % More accurate but it is slow
end
end