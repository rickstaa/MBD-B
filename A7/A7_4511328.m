%% MBD_B: Assignment 7 - Quick return mechanism
%  Rick Staa (4511328)
%  Last edit: 09/05/2018
clear all;tic; close all; % clc;
fprintf('--- A7 ---\n');

%% Script parameters
parms.accuracy_bool = 0;                    % If set to 1 A\b will be performed instead of inv(A)*B this is more accurate but slower
animate_bool = 0;                           % Set on 1 if you want to see an animation

%% Set up needed symbolic parameters
% Create needed symbolic variables
syms phi2 phi4 phi5 phi2d phi4d phi5d

% Put in parms struct for easy function handling
parms.syms.phi2             = phi2;
parms.syms.phi4             = phi4;
parms.syms.phi5             = phi5;
parms.syms.phi2d            = phi2d;
parms.syms.phi4d            = phi4d;
parms.syms.phi5d            = phi5d;

%% Intergration parameters
time                        = 3;                                        % Intergration time
parms.h                     = 1e-3;                                     % Intergration step size
parms.tol                   = 1e-12;                                    % Intergration constraint error tolerance
parms.nmax                  = 10;                                       % Maximum number of Gauss-Newton drift correction iterations
parms.rot_sim               = 4*pi;                                     % How many rotations you want to simulate

%% Model Parameters
% Lengths and distances
parms.O2A                   = 0.2;                                      % Length segment 2 [m]
parms.O4B                   = 0.7;                                      % Length segment 4 [m]
parms.BC                    = 0.6;                                      % Length segment 5 [m]
parms.O4O2                  = 0.3;                                      % Distance between joint 4 and joint 2 [m]
parms.O4G4                  = 0.4;                                      % Distance bewteen COM4 and joint 4 [m]
parms.BG5                   = 0.3;                                      % Distance joint 5 and COM 5 [m]
parms.Yc                    = 0.9;                                      % Height joint C (COM body 6) [m]

% Masses and inertias
parms.m3                    = 0.5;                                      % Body 3 weight [kg]
parms.m4                    = 6;                                        % Body 4 weight [kg]
parms.m5                    = 4;                                        % Body 5 weight [kg]
parms.m6                    = 2;                                        % Body 6 weight [kg]
parms.J2                    = 100;                                      % Moment of inertia body 2 [kgm^2]
parms.J3                    = 0;                                        % Moment of inertia body 3 [kgm^2] - Put on 0 because no moment possible
parms.J4                    = 10;                                       % Moment of inertia body 4 [kgm^2]
parms.J5                    = 6;                                        % Moment of inertia body 5 [kgm^2]

%% World parameters
% Gravity
parms.g                     = 0;                                        % [parms.m/s^2]

% Forces
% Add forces F=[M2,F3_x,F3_y,M3,F4_x,F4_y,M4,F5_x,F5_y,M5,F6_x];
parms.F6_x                  = 1000;                                     % x force on body 6 [N]%
parms.T2                    = 0;                                        % Torque around joint 6 [Nm]
parms.F                     = [parms.T2, 0, -parms.m3*parms.g, 0, 0, -parms.m4*parms.g, 0, 0, -parms.m5*parms.g, 0, parms.F6_x];
parms.Q                     = [0;0;0];                                  % The generalised forces and torques

%% Calculate Initial states
% Set the independent states
phi2_init                   = 0;
phi2d_init                  = (150*pi)/60;

% Calculate the dependen states
x3_init                     = parms.O2A*cos(phi2_init);                                               % Calculate the x coordinate of slider 3
y3_init                     = parms.O4O2+parms.O2A*sin(phi2_init);                                    % Calculate the y coordinate of slider 3
phi4_init                   = atan2(y3_init,x3_init);
phi5_init                   = pi-asin((parms.Yc-parms.O4B*sin(phi4_init))/parms.BC);
phi4d_init                  = cos(phi4_init)^2*phi2d_init;                                       % Not real value but failed to calculate
phi5d_init                  = (parms.O4B*cos(phi4_init)*phi4d_init)/(-parms.BC*cos(phi5_init));  % Not real value but failed to calculate
q0                          = [phi2_init phi4_init phi5_init phi2d_init phi4d_init phi5d_init];

%% Derive equation of motion
EOM_calc(parms);                                                        % Calculate symbolic equations of motion and put in parms struct

%% Calculate movement by mean sof a Runge-Kuta 4th order intergration method
[t_RK4,q_RK4,x_RK4,xdp_RK4]                         = RK4_custom(q0,parms);

%% Calculate com velocities
xp = diff(x_RK4)/parms.h;

%% Create plots

%% Plot Angular speed crank as a function of time
figure;
plot(t_RK4,q_RK4(:,4:6),'linewidth',1.5);
set(gca,'fontsize',18);
title('Angular speed of the crank 2, rocker 4 and connecting bar 5');
xlabel('Time [s]');
ylabel('Angular speed [rad/s]');
legend('Crank 2 (\phi_2)','Rocker 4 (\phi_4)','Connecting bar 5 (\phi_5)','Location', 'Best');

%% Plot the sliding speedof slider 3 with respect to rocker 4
v_slider_rel = xp(2:end,4).*cos(q_RK4(3:end,3)) + xp(2:end,5).*sin(q_RK4(3:end,3));

figure;
plot(t_RK4,xdp_RK4(:,end),'linewidth',1.5);
set(gca,'fontsize',18);
xlabel('Time [s]');
ylabel('Acceleration slider 6 [m/s^2]');
title('Acceleration slider 6');
legend('acceleration slider 6','Location', 'Best');

figure;
plot(t_RK4,x_RK4(:,end),'linewidth',1.5);
set(gca,'fontsize',18);
xlabel('Time [s]');
ylabel('position slider 6 [m]');
title('position slider 6');
legend('position slider 6','Location', 'Best');

figure;
plot(t_RK4(4:end),xp(3:end,end),'linewidth',1.5);
set(gca,'fontsize',18);
xlabel('Time [s]');
ylabel('velocity slider 6 [m]');
title('velocity slider 6');
legend('velocity slider 6','Location', 'Best');

plot(t_RK4(4:end),v_slider_rel(2:end),'linewidth',1.5);
set(gca,'fontsize',18);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Speed of slider 3 with respect to rocker 4');
legend('relative velocity','Location', 'Best');

%% Normal forces
figure;
plot(t_RK4,q_RK4(:,end-1:end),'linewidth',1.5);
set(gca,'fontsize',18);
xlabel('Time [s]');
ylabel('Force [N]');
title('Reaction Forces [N]');
legend('slider 6 on ground','slider 3 on rocker 4','Location', 'Best')
toc

%% -- ANIMATE --
if animate_bool == 1
    % Adapted from A. Schwab's animation code and a animation script posted
    % on github by BitNide
    y = q_RK4.';
    g = figure;
    filename = 'testAnimated.gif';
    
    %% Loop through all time steps
    l = line([0, 0.2*cos(y(1,1))],[0.3, 0.3+0.2*sin(y(1,1))]);
    k = line([0, 0.7*cos(y(2,1))],[0, 0.7*sin(y(2,1))]);
    j = line([0.7*cos(y(2,1)), 0.7*cos(y(2,1)) + cos(y(3,1))*0.6],[0.7*sin(y(2,1)), 0.7*sin(y(2,1))+sin(y(3,1))*0.6]);
    xlim([-1.5 1])
    ylim([0 1])
    
    nstep = length(y);
    nskip = 10;
    
    for istep = 2:nskip:nstep
        set(l,'XData',[0, 0.2*cos(y(1,istep))])
        set(l,'YData',[0.3, 0.3+0.2*sin(y(1,istep))])
        set(k,'XData',[0, 0.7*cos(y(2,istep))])
        set(k,'YData',[0, 0.7*sin(y(2,istep))])
        set(j,'XData',[0.7*cos(y(2,istep)), 0.7*cos(y(2,istep)) + cos(y(3,istep))*0.6])
        set(j,'YData',[0.7*sin(y(2,istep)), 0.7*sin(y(2,istep))+sin(y(3,istep))*0.6])
        drawnow
        pause(parms.h)
    end
end

%% FUNCTIONS

%% Runge-Kuta numerical intergration function
% This function calculates the motion of the system by means of a
% Runge-Kuta numerical intergration. This function takes as inputs the
% parameters of the system (parms), the EOM of the system (parms.EOM)
% and the initial state.
function [t,q,x,xdd] = RK4_custom(q0,parms)

% Calculate x0
q_new_tmp               = num2cell(q0,1);
x0                      = subs_x(q_new_tmp{1:3}).';
xdd0                    = subs_xdp(q_new_tmp{:}).';

% Initialise variables
t                       = 0;                                                 % Initiate time
q                       = [q0 0 0 0 0 0];                                    % Put initial state in array
x                       = x0;
xdd                     = xdd0;
% Caculate the motion for the full simulation time by means of a
% Runge-Kutta4 method

% Perform intergration till two full rotations of the crank
ii = 1;                                                                      % Create counter
while abs(q(ii,1)) < parms.rot_sim
    
    % Calculate the next state by means of a RK4 method
    q_now_tmp           = num2cell(q(ii,1:end-5),1);                                                % Create cell for subs function function
    K1                  = [cell2mat(q_now_tmp(1,end-2:end)),subs_qdp(q_now_tmp{:}).'];              % Calculate the second derivative at the start of the step
    q1_tmp              = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K1(1,1:end-2));              % Create cell for subs function function
    K2                  = [cell2mat(q1_tmp(1,end-2:end)),subs_qdp(q1_tmp{:}).'];                    % Calculate the second derivative halfway the step
    q2_tmp              = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K2(1:end-2));                % Refine value calculation with new found derivative
    K3                  = [cell2mat(q2_tmp(1,end-2:end)),subs_qdp(q2_tmp{:}).'];                    % Calculate new derivative at the new refined location
    q3_tmp              = num2cell(cell2mat(q_now_tmp) + (parms.h)*K3(1:end-2));                    % Calculate state at end step with refined derivative
    K4                  = [cell2mat(q3_tmp(1,end-2:end)),subs_qdp(q3_tmp{:}).'];                    % Calculate last second derivative
    q_now_d             = (1/6)*(K1(end-4:end)+2*K2(end-4:end)+2*K3(end-4:end)+K4(end-4:end));      % Estimated current derivative
    q_next              = cell2mat(q_now_tmp) + (parms.h/6)*(K1(1:6)+2*K2(1:6)+2*K3(1:6)+K4(1:6));  % Perform euler intergration step
    
    % Save reaction forces and current derivative in state
    q(ii,end-4:end)     = q_now_d;
    
    % Save full state back in q array
    q                   = [q;[q_next 0 0 0 0 0]];
    
    % Correct for intergration drift
    q_now_tmp = q(ii+1,:);
    [q_new,~] = gauss_newton(q_now_tmp,parms);
    
    % Update the second derivative and the constraint forces
    q_new_tmp           = num2cell(q(ii,1:end-5),1);
    q_update            = subs_qdp(q_new_tmp{:}).';
    
    % Overwrite position coordinates
    q(ii+1,:)           = [q_new(1:6) q_update];
    
    % Create time array
    t                   = [t;t(ii)+parms.h];                                   % Perform Gauss-Newton drift correction
    ii                  = ii + 1;                                              % Append counter
    
    % Calculate COM coordinates
    % Calculate COM coordinates
    x_tmp               = subs_x(q_new_tmp{1:3}).';
    xdd_tmp             = subs_xdp(q_new_tmp{:}).';
    
    % Save x in state
    x                   = [x;x_tmp];
    xdd                 = [xdd;xdd_tmp];
    
end
end

%% Constraint calculation function
function [C,Cd] = constraint_calc(q)

% Get needed angles out
q_now_tmp       = num2cell(q,1);

% Calculate the two needed constraints
C               = subs_C(q_now_tmp{1:3});

% Calculate constraint derivative
Cd              = subs_Cd(q_now_tmp{1:3}).';

end

%% Speed correct function
function [q,error] = gauss_newton(q,parms)

% Get rid of the drift by solving a non-linear least square problem by
% means of the Gaus-Newton method
% Calculate the two needed constraints
[C,Cd] = constraint_calc(q);

%% Guass-Newton position correction
n_iter          = 0;                                                                        % Set iteration counter

% Solve non-linear constraint least-square problem
while (max(abs(C)) > parms.tol)&& (n_iter < parms.nmax)
    q_tmp           = q(1:3);
    n_iter          = n_iter + 1;
    q_del           = Cd*inv(Cd.'*Cd)*-C;
    q(1:3)          = q_tmp+ q_del.';
    
    % Recalculate constraint
    [C,Cd]          = constraint_calc(q);
end

% Calculate the corresponding speeds
q_tmp_vel          = q(4:6);
Dqd_n1             = -Cd*inv(Cd.'*Cd)*Cd.'*q_tmp_vel.';
q(4:6)             = q_tmp_vel + Dqd_n1.';

error = C;
end

%% Calculate (symbolic) Equations of Motion four our setup
function EOM_calc(parms)

% Unpack symbolic variables from varargin
phi2            = parms.syms.phi2;
phi4            = parms.syms.phi4;
phi5            = parms.syms.phi5;
phi2d           = parms.syms.phi2d;
phi4d           = parms.syms.phi4d;
phi5d           = parms.syms.phi5d;

% Create generalized coordinate vectors
q               = [phi2;phi4;phi5];
qd              = [phi2d;phi4d;phi5d];

% COM of the bodies expressed in generalised coordinates
x2              = 0;
y2              = parms.O4O2;
x3              = parms.O2A*cos(phi2);
y3              = parms.O4O2+parms.O2A*sin(phi2);
x4              = parms.O4G4*cos(phi4);
y4              = parms.O4G4*sin(phi4);
x5              = parms.O4B*cos(phi4)+parms.BG5*cos(phi5);
y5              = parms.O4B*sin(phi4)+parms.BG5*sin(phi5);
x6              = parms.O4B*cos(phi4)+parms.BC*cos(phi5);
y6              = parms.O4B*sin(phi4)+parms.BC*sin(phi5);

% Create mass matrix
% x2 = 0, y2 = 0 and y5 =0 also no moments around slider 3 and 6
parms.M         = diag([parms.J2,parms.m3,parms.m3,parms.J3,parms.m4,parms.m4,parms.J4,parms.m5,parms.m5,parms.J5,parms.m6]);

% Put in one state vector
x               = [phi2;x3;y3;phi4;x4;y4;phi4;x5;y5;phi5;x6];

% Calculate the two needed constraints
C               = [x3-sqrt(x3^2+y3^2)*cos(phi4);            ...
    y6-parms.Yc];

% Compute the jacobian of state and constraints
Jx_q            = simplify(jacobian(x,q.'));
JC_q            = simplify(jacobian(C,q));

%% Calculate convective component
Jx_dq           = jacobian(Jx_q*qd,q);
JC_dq           = jacobian(JC_q*qd,q);

% Solve with virtual power
M_bar           = simplify(Jx_q.'*parms.M*Jx_q);

% Create system of DAE
A               = [M_bar JC_q.'; JC_q zeros(size(JC_q,1))];
B               = [parms.Q + Jx_q.'*(parms.F.'-parms.M*Jx_dq*qd); ...
    -JC_dq*qd];

% Calculate result expressed in generalized coordinates
if parms.accuracy_bool == 0
    qdd             = inv(A)*B;        % Less accurate but in our case faster
else
    qdd             = A\B;             % More accurate but it is slow
end

% % Get result back in COM coordinates
% xdp     = simplify(jacobian(xp,qp.'))*qdp + simplify(jacobian(xp,q.'))*qp;

%% Convert to function handles
% xdp_handle       = matlabFunction(xdp);                                        % Create function handle of EOM in terms of COM positions
matlabFunction(simplify(qdd),'File','subs_qdp');                                 % Create function handle of EOM in terms of generalised coordinates

% Constraint function handle
matlabFunction(simplify(C),'File','subs_C');                                     % Create function handle of EOM in terms of generalised coordinates

% Constraint derivative function handle
matlabFunction(simplify(JC_q),'File','subs_Cd');                                 % Create function handle of EOM in terms of generalised coordinates

% Get back to COM coordinates
matlabFunction(simplify(x),'File','subs_x');                                     % Create function handle of EOM in terms of generalised coordinates

% Get xdp COM coordinates
xd              = Jx_q*qd;
xdd             = simplify(jacobian(xd,qd))*qdd(1:3)+simplify(jacobian(xd,q))*qd(1:3);
matlabFunction(simplify(xdd),'File','subs_xdd');                                 % Create function handle of EOM in terms of generalised coordinates

end