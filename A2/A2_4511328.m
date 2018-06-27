%% MBD_B: Assignment 2 - Double pendulum systemetic approach
%  Rick Staa (4511328)
%  Last edit: 05/03/2018
clear all; close all; clc;
fprintf('--- A2 ---\n');

%% Script settings
parms.accuracy_bool = 0;                                        % If set to 1 A\b will be performed instead of inv(A)*B this is more accurate but slower
parms.booleans.ex_constr = 0;                                   % Put on 1 if you want to enable the extra constraint

%% Parameters
% Segment 1
parms.L1     = 0.55;                                            % [m]
parms.w1     = 0.05;                                            % [m]
parms.t1     = 0.004;                                           % [m]
parms.p1     = 1180;                                            % [kg/m^3]
parms.m1     = parms.p1 * parms.w1 * parms.t1 * parms.L1;       % [kg]
parms.I1     = (1/12) * parms.m1 * parms.L1^2;                  % [kg*m^2]

% Segment 2
parms.L2     = 0.55;                                            % [m]
parms.w2     = 0.05;                                            % [m]
parms.t2     = 0.004;                                           % [m]
parms.p2     = 1180;                                            % [m]
parms.m2     = parms.p2 * parms.w2 * parms.t2 * parms.L2;       % [kg]
parms.I2     = (1/12) * parms.m2 * parms.L2^2;                  % [kg*m^2]

% World parameters
parms.g      = 9.81;                                            % [m/s^2]

% Create state space matrices
[parms]      = create_state(parms);

%% Question 1: Initial states
% Set Force vector
parms.Fg     = [parms.m1*parms.g;0;0;parms.m2*parms.g;0;0];

% b) Both bars vertical up and zero speed
x0           = [0.5*pi 0.5*pi 0 0];
x_dd.b       = state_calc(x0,parms);

% c) Both bars horizontal to the right and zero speed
x0           = [0 0 0 0];
x_dd.c       = state_calc(x0,parms);

% d) Both bars horizontal and with an initial angular speed on both bars of 60
% rpm
x0           = [0 0 2*pi 2*pi];
x_dd.d       = state_calc(x0,parms);

% Calculate velocities
x1_d         = -(parms.L1/2)*sin(x0(1))*x0(3);
y1_d         =  (parms.L1/2)*cos(x0(1))*x0(3);
x2_d         = x1_d - (parms.L1/2)*sin(x0(1))*x0(3) - (parms.L2/2)*sin(x0(2))*x0(4);
y2_d         = y1_d + (parms.L1/2)*cos(x0(1))*x0(3) + (parms.L2/2)*cos(x0(2))*x0(4);

% Questions e-g Extra constraint
parms.booleans.ex_constr = 1;                    % Put on 1 if you want to enable the extra constraint
[parms] = create_state(parms);                   % Calculate new state matrixes

% e - Both bars vertical up and zero speed
x0           = [0.5*pi 0.5*pi 0 0];
x_dd.e       = state_calc(x0,parms);

% f - Both bars vertical up and angular speed of omega = 60 rmp on bar 1
% Calculate initial states
phi1_0  = 0.5*pi; 
phi2_0  = 0.5*pi; 
phi1d_0 = -2*pi;
phi2d_0 = (-parms.L1*sin(phi1_0)*phi1d_0)/(parms.L2*sin(phi2_0));
x0      = [phi1_0 phi2_0 phi1d_0 phi2d_0];
x_dd.f  = state_calc(x0,parms);

% f - Both bars horizontal with an initial angular speed of omega = 60
% rpm on bar 1
parms.Fg    = [parms.m1*parms.g;0;-10*(parms.L1/2);parms.m2*parms.g;0;0];        % Add a 10 N force in the y direction in B
x0          = [0 pi 2*pi 0];
x_dd.g      = state_calc(x0,parms);

% Test rank of Cx matrix
Cx        = double(vpa(subs(parms.Cx,{'phi1','phi2','phid1','phid2' 'm1' 'm2' 'I1' 'I2' 'L1' 'L2'},[x0(1) x0(2) x0(3) x0(4) parms.m1 parms.m2 parms.I1 parms.I2 parms.L1 parms.L2])));
rank_cx   = rank(Cx);
null_cx   = null(Cx);

%% Functions
% -- Create_state space --
% This function calculates the state matrixes for our double pendulum
% simulation

function [parms] = create_state(parms)

% Create symbolic variables
syms x1 y1 phi1 x2 y2 phi2                      % States
syms xd1 yd1 phid1 xd2 yd2 phid2                % State derivatives (dx/dt)
syms L1 L2 m1 m2 I1 I2 g                        % Parameters

% put the cm cogordinates and their time derivatives in a column vector x
x       = [x1; y1; phi1; x2; y2; phi2];
xd      = [xd1; yd1; phid1; xd2; yd2; phid2];

% Create constraints
ck_x1     = x1-(L1/2)*cos(phi1);                                  % X constraint on first body
ck_y1     = y1-(L1/2)*sin(phi1);                                  % Y constraint on second body
ck_x2     = (x2-L2/2*cos(phi2))-(x1+L1/2*cos(phi1));              % X constraint on second body
ck_y2     = (y2-L2/2*sin(phi2))-(y1+L1/2*sin(phi1));              % Y constraint on second body
if parms.booleans.ex_constr == 1
    ck_x2c    = x2+(L2/2)*cos(phi2);                              % Extra constraint on the right end of bar 2 (vertical path following at orgin)
end

% Calculate the jacobian (This can be used for the constraint equations of
% motion
if parms.booleans.ex_constr == 1
    C         = [ck_x1;ck_y1;ck_x2;ck_y2;ck_x2c];
else
    C         = [ck_x1;ck_y1;ck_x2;ck_y2];
end
Cx        = jacobian(C,x);                                   % Calculate partial derivative
Cx        = simplify(Cx);

% Calculate second derivative this is done with the jacobian and hessian
% (Chain rule)
Cd        = Cx*xd;
Cdp       = simplify(jacobian(Cd,x)*xd);                    % Convective term

% Save constraint jacobian and convective terms
parms.Cx = Cx;
parms.Cdp = Cdp;

end

% -- State Calc --
% This function calculates the second derivative of the double pendulum
% using the initial states.
function [xdd] = state_calc(x0,parms)

% Now build the system
M         = diag([parms.m1 parms.m1 parms.I1 parms.m2 parms.m2 parms.I2]);
Fg        = parms.Fg;
A         = [M,parms.Cx';parms.Cx,zeros(size(parms.Cx,1),size(parms.Cx',2))];
b         = [parms.Fg;-parms.Cdp];

% Substitude initial states in derived matrices
A        = double(vpa(subs(A,{'phi1','phi2','phid1','phid2' 'm1' 'm2' 'I1' 'I2' 'L1' 'L2'},[x0(1) x0(2) x0(3) x0(4) parms.m1 parms.m2 parms.I1 parms.I2 parms.L1 parms.L2])));
b        = double(vpa(subs(b,{'phi1','phi2','phid1','phid2' 'm1' 'm2' 'L1' 'L2' 'g'},[x0(1) x0(2) x0(3) x0(4) parms.m1 parms.m2 parms.L1 parms.L2 parms.g])));

% Calculate second derivative
if parms.accuracy_bool == 0
    xdd = inv(A)*b;             % Less accurate but in our case faster
else
    xdd = A\b;                  % More accurate but it is slow
end
end