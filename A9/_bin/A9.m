%% MBD_B: Assignment 9 - Euler angles and the human arm
%  Rick Staa (4511328)
clear all; close all; clc;
fprintf('--- A9 ---\n');

%% Set up needed symbolic parameters
% Create needed symbolic variables
syms alpha beta gamma alpha_d beta_d gamma_d alpha_dd beta_dd gamma_dd

% Put in parms struct for easy function handling
parms.syms.alpha             = alpha;
parms.syms.beta              = beta;
parms.syms.gamma             = gamma;
parms.syms.alpha_d           = alpha_d;
parms.syms.beta_d            = beta_d;
parms.syms.gamma_d           = gamma_d;
parms.syms.alpha_dd          = alpha_dd;
parms.syms.beta_dd           = beta_dd;
parms.syms.gamma_dd          = gamma_dd;

%% Intergration parameters
time                        = 5;                                        % Intergration time
parms.h                     = 1e-3;                                     % Intergration step size

%% Model Parameters
% Lengths and distances
parms.L1                     = 0.3;                                     % Length upper arm [m]
parms.L2                     = 0.4;                                     % Length lower arm [m]
parms.m1                     = 3;                                       % Mass upper arm [kg] 
parms.m2                     = 3;                                       % Mass lower arm [kg]

%% World parameters
% Gravity
parms.g                     = 9.81;                                     % [parms.m/s^2]

%% Set Initial states
alpha_0                     = deg2rad(30);
beta_0                      = deg2rad(-20);
gamma_0                     = deg2rad(-20);
alpha_d_0                   = 0;
beta_d_0                    = 0;
gamma_d_0                   = 0;
q0                          = [alpha_0;beta_0;gamma_0;alpha_d_0;beta_d_0;gamma_d_0];    % Put in initial state vector

%% Calculate equilibrium torques
torque_calc(parms);
q0_tmp                      = num2cell(q0',1);
parms.Q                     = subs_torque(q0_tmp{:});

%% Derive equation of motion
EOM_calc(parms);                                                                         % Calculate symbolic equations of motion and put in parms struct

%% Calculate movement with ode
opt                         = odeset('AbsTol',1e-6,'RelTol',1e-6,'Stats','on');
[t,q]                       = ode113(@(t,q) ODE_func(t,q), [0 time], q0',opt);

%% Animate movement
% animator(t,q,parms);

%% Plots

% Create plot label
title_str = strcat("\alpha = ", num2str(rad2deg(q0(1))),", \beta = ", num2str(rad2deg(q0(2))),"& \gamma = ", num2str(rad2deg(q0(3))));

% Plot angles vs time
figure;
subplot(3,1,1);
f1 = plot(t,q(:,1),'Linewidth',1.5);
ylim([-2*pi 2*pi])
x_lim = xlim;
x_point = (x_lim(2)-x_lim(1))/2+x_lim(1);
text(x_point,mean(q(:,1))+0.4*pi,(strcat('var =',num2str(var(q(:,1))))));
xlabel('t [s]');
ylabel('\alpha ]rad]');
title(title_str);
subplot(3,1,2)
f2 = plot(t,q(:,2),'Linewidth',1.5);
text(x_point,mean(q(:,2))+0.4*pi,(strcat('var =',num2str(var(q(:,2))))));
ylim([-2*pi 2*pi])
xlabel('t [s]');
ylabel('\beta ]rad]');
subplot(3,1,3)
f3 = plot(t,q(:,3),'Linewidth',1.5);
text(x_point,mean(q(:,3))+0.4*pi,(strcat('var =',num2str(var(q(:,3))))));
ylim([-2*pi 2*pi])
xlabel('t [s]');
ylabel('\gamma ]rad]');

%% FUNCTIONS

%% ANIMATION FUNCTION
function animator(t,q,parms)
%% Animation
% Adapted from A. Schwab's animation code

% Calculate elbow and wrist
x_full              = zeros(size(q,1),12);
for ii = 1:size(q,1)
    q_now            = num2cell(q(ii,1:3),1);
    x_full(ii,:)     = subs_x_full(q_now{:})';
end

% Create figure
figure
set(gca,'fontsize',16)
title('Animation arm')
shoulder          = plot3(0,0,0,'*g');
hold on
upper_arm         = plot3([0 x_full(1,7)],[0 x_full(1,8)],[0 x_full(1,9)],'-b');
elbow             = plot3(x_full(1,7),x_full(1,8),x_full(1,9),'*c'); 
lower_arm         = plot3([x_full(1,7) x_full(1,10)],[x_full(1,8) x_full(1,11)],[x_full(1,9) x_full(1,12)],'-r');
wrist             = plot3(x_full(1,10),x_full(1,11),x_full(1,12),'*m'); 
set(upper_arm,'LineWidth',5);
set(upper_arm,'Color','b')
set(lower_arm,'LineWidth',5);
set(lower_arm,'Color','r')
set(shoulder,'Linewidth',10);
set(elbow,'Linewidth',10);
set(wrist,'Linewidth',10);
axis([-(parms.L1+parms.L2) (parms.L1+parms.L2) -(parms.L1+parms.L2) (parms.L1+parms.L2) -(parms.L1+parms.L2) (parms.L1+parms.L2)]);
legend('shoulder','upper arm','elbow','lower arm','wrist');
nstep = length(t);
nskip = 10;
for istep = 2:nskip:nstep
    set(upper_arm,'XData',[0 x_full(istep,7)])
    set(upper_arm,'YData',[0 x_full(istep,8)])
    set(upper_arm,'ZData',[0 x_full(istep,9)])
    set(elbow,'XData',x_full(istep,7))
    set(elbow,'YData',x_full(istep,8))
    set(elbow,'ZData',x_full(istep,9))
    set(lower_arm,'XData',[x_full(istep,7) x_full(istep,10)])
    set(lower_arm,'YData',[x_full(istep,8) x_full(istep,11)])
    set(lower_arm,'ZData',[x_full(istep,9) x_full(istep,12)])    
    set(wrist,'XData',x_full(istep,10))
    set(wrist,'YData',x_full(istep,11))
    set(wrist,'ZData',x_full(istep,12))
    axis([-(parms.L1+parms.L2) (parms.L1+parms.L2) -(parms.L1+parms.L2) (parms.L1+parms.L2) -(parms.L1+parms.L2) (parms.L1+parms.L2)]);
    drawnow
    pause(1e-10)
end
end

%% ODE function handle
function [qdp] = ODE_func(t,q)
    q_now = num2cell(q',1);
    qdp   = subs_qdp(q_now{1:end});
    qdp   = [q(4);q(5);q(6);qdp];
end

%% Calculate COM torques
function torque_calc(parms)

% Unpack symbolic variables from varargin
alpha           = parms.syms.alpha;
beta            = parms.syms.beta;
gamma           = parms.syms.gamma;
alpha_d         = parms.syms.alpha_d;
beta_d          = parms.syms.beta_d;
gamma_d         = parms.syms.gamma_d;

% Create generalized coordinate vectors
q               = [alpha;beta;gamma];
qd              = [alpha_d;beta_d;gamma_d];

% Create rotation matrices
R_alpha         = rot_x(alpha);
R_beta          = rot_y(beta);
R_gamma         = rot_x(gamma);

% express COM positions in inertial frame and put them in a vector
r1_I            = R_alpha*R_beta*[0;0;-parms.L1/3];
r2_I            = R_alpha*R_beta*[0;0;-parms.L1]+R_alpha*R_beta*R_gamma*[0;0.5*parms.L2;0];

% Create mass matrix
parms.M         = diag([parms.m1,parms.m1,parms.m1,parms.m2,parms.m2,parms.m2]);

% Put in one state vector
x               = [r1_I;r2_I];

% Compute the jacobian of state and constraints
Jx_q            = simplify(jacobian(x,q.'));

%% Calculate convective component
Jx_dq           = jacobian(Jx_q*qd,q);

% Add forces F=[M2,F3_x,F3_y,M3,F4_x,F4_y,M4,F5_x,F5_y,M5,F6_x];
F_B             = [0;0;-parms.m1*parms.g;0;0;-parms.m2*parms.g];

% Calculate result expressed in generalized coordinates
Q_0 = (Jx_q.'*(F_B-parms.M*Jx_dq*qd));

% Get x_full to animate the body
matlabFunction(simplify(Q_0),'vars',[alpha,beta,gamma,alpha_d,beta_d,gamma_d],'file','subs_torque');                           % Create function handle of EOM in terms of generalised coordinates

end

%% Calculate (symbolic) Equations of Motion four our setup
function EOM_calc(parms)

% Unpack symbolic variables from varargin
alpha           = parms.syms.alpha;
beta            = parms.syms.beta;
gamma           = parms.syms.gamma;
alpha_d         = parms.syms.alpha_d;
beta_d          = parms.syms.beta_d;
gamma_d         = parms.syms.gamma_d;
alpha_dd        = parms.syms.alpha_dd;
beta_dd         = parms.syms.beta_dd;
gamma_dd        = parms.syms.gamma_dd;

% Create rotation matrices
R_alpha         = rot_x(alpha);
R_beta          = rot_y(beta);
R_gamma         = rot_x(gamma);

% Create generalized coordinate vectors
q               = [alpha;beta;gamma];
qd              = [alpha_d;beta_d;gamma_d];
qdd             = [alpha_dd;beta_dd;gamma_dd];

% express COM positions in inertial frame and put them in a vector
r1_I            = R_alpha*R_beta*[0;0;-parms.L1/3];
r2_I            = R_alpha*R_beta*[0;0;-parms.L1]+R_alpha*R_beta*R_gamma*[0;0.5*parms.L2;0];

% Create mass matrix
parms.M         = diag([parms.m1,parms.m1,parms.m1,parms.m2,parms.m2,parms.m2]);

% Put in one state vector
x               = [r1_I;r2_I];

% Create full state for animation
x_elbow         = R_alpha*R_beta*[0;0;parms.L1];
x_wrist         = R_alpha*R_beta*[0;0;parms.L1]+R_alpha*R_beta*R_gamma*[0;parms.L2;0];
x_full          = [x;x_elbow;x_wrist];

% Compute the jacobian of state and constraints
Jx_q            = simplify(jacobian(x,q.'));

%% Calculate convective component
Jx_dq           = jacobian(Jx_q*qd,q);

% Solve with virtual power
M_bar           = Jx_q.'*parms.M*Jx_q;

% Add forces F=[M2,F3_x,F3_y,M3,F4_x,F4_y,M4,F5_x,F5_y,M5,F6_x];
F_B             = [0;0;-parms.m1*parms.g;0;0;-parms.m2*parms.g];
F               = Jx_q.'*(F_B-parms.M*Jx_dq*qd)-parms.Q;                 % Forces expressed in the inertial frame in generalised coordinates

% Calculate result expressed in generalized coordinates
qdp             = M_bar\F;

% Get result back in COM coordinates
xd              = Jx_q*qd;
xdd             = simplify(jacobian(xd,qd.'))*qdd + simplify(jacobian(xd,q.'))*qdd;

%% Convert to function handles
matlabFunction(simplify(qdp),'vars',[alpha,beta,gamma,alpha_d,beta_d,gamma_d],'File','subs_qdp');                                 % Create function handle of EOM in terms of generalised coordinates

% Get back to COM coordinates
matlabFunction(simplify(x),'File','subs_x');                                     % Create function handle of EOM in terms of generalised coordinates

% Get xdp COM coordinates
matlabFunction(simplify(xdd),'vars',[alpha,beta,gamma,alpha_d,beta_d,gamma_d,alpha_dd,beta_dd,gamma_dd],'File','subs_xdd');                                 % Create function handle of EOM in terms of generalised coordinates

% Get x_full to animate the body
matlabFunction(simplify(x_full),'file','subs_x_full');                           % Create function handle of EOM in terms of generalised coordinates
end

%% Rotation matrices
% x rotation matrix
function R_alpha = rot_x(angle)
    R_alpha = [   1           0           0          ;...
                  0           cos(angle)    -sin(angle);...
                  0           sin(angle)    cos(angle)];
end

% y rotation matrix
function R_beta = rot_y(angle)
    R_beta = [ cos(angle)            0           sin(angle);...
                 0                   1           0;...
               -sin(angle)           0           cos(angle)];
end

% z rotation matrix
function R_gamma = rot_z(angle)
    R_gamma = [cos(angle)            -sin(angle)           0;...
               sin(angle)            cos(angle)            0;...
                   0                   0                   1];
end