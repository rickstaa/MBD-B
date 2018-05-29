%% MBD_B: Assignment 8 - EzyRoller
%  Rick Staa (4511328)
%  Last edit: 29/05/2018
clear all; % close all; clc;
fprintf('--- A8 ---\n');

%% Set up needed symbolic parameters
% Create needed symbolic variables
syms x1 y1 phi1 phi2 x1d y1d phi1d phi2d

% Generalised state
parms.syms.x1               = x1;
parms.syms.y1               = y1;
parms.syms.phi1             = phi1;
parms.syms.phi2             = phi2;

% Generalised state derivative
parms.syms.x1d              = x1d;
parms.syms.y1d              = y1d;
parms.syms.phi1d            = phi1d;
parms.syms.phi2d            = phi2d;

%% Intergration parameters
time                        = 3;                                        % Intergration time
parms.h                     = 1e-3;                                     % Intergration step size
parms.tol                   = 1e-12;                                    % Intergration constraint error tolerance
parms.nmax                  = 10;                                       % Maximum number of Gauss-Newton drift correction iterations

%% Model Parameters
% Lengths and distances
parms.a                     = 0.5;                                      % Length wheel first segment to COM segment 1
parms.b                     = 0.5;                                      % Length COM to revolute joint B
parms.c                     = 0.125;                                    % Length revolute jonit B to COM segment 2
parms.d                     = 0.125;                                    % Length COM segment 2 to wheel 2

% Masses and inertias
parms.m1                    = 1;                                        % Body 1 weight [kg]
parms.m2                    = 0;                                        % Body 2 weight [kg]
parms.J1                    = 0.1;                                      % Moment of inertia body 1 [kgm^2]
parms.J2                    = 0;                                        % Moment of inertia body 2 [kgm^2]

%% World parameters
% Gravity
parms.g                     = 9.81;                                     % [parms.m/s^2]

%% Calculate Initial states
x1_init                     = 0.5;
y1_init                     = 0.5;
phi1_init                   = 0.5;
phi2_init                   = 0.5;
x1d_init                    = 0;
y1d_init                    = 0;
phi1d_init                  = 0;
phi2d_init                  = 0;
q0                          = [x1_init y1_init phi1_init phi2_init x1d_init y1d_init phi1d_init phi2d_init];

%% Derive equation of motion
[EOM_qdp,C_handle,Cd_handle,X_handle,Xp_handle] = EOM_calc(parms);                         % Calculate symbolic equations of motion and put in parms struct
parms.C_handle               = C_handle;
parms.Cd_handle              = Cd_handle;
parms.X_handle               = X_handle;
parms.Xp_handle              = Xp_handle;

%% Calculate movement by mean sof a Runge-Kuta 4th order intergration method
tic
[t_RK4,q_RK4,x_RK4,xdp_RK4]                         = RK4_custom(EOM_qdp,q0,parms);
toc

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

%% FUNCTIONS

%% Runge-Kuta numerical intergration function
% This function calculates the motion of the system by means of a
% Runge-Kuta numerical intergration. This function takes as inputs the 
% parameters of the system (parms), the EOM of the system (parms.EOM) 
% and the initial state.
function [t,q,x,xdp] = RK4_custom(EOM,q0,parms)

% Calculate x0
q_new_tmp        = num2cell(q0,1);
x0   = feval(parms.X_handle,q_new_tmp{1:3}).';
xdp0 = feval(parms.Xp_handle,q_new_tmp{:}).';

% Initialise variables
t                     = 0;                                                 % Initiate time
q                     = [q0 0 0 0 0 0];                                    % Put initial state in array
x                     = x0;
xdp                   = xdp0;
% Caculate the motion for the full simulation time by means of a
% Runge-Kutta4 method

% Perform intergration till two full rotations of the crank
ii = 1;                                                                    % Create counter
while abs(q(ii,1)) < (4*pi)
    
    % Calculate the next state by means of a RK4 method
    q_now_tmp         = num2cell(q(ii,1:end-5),1);                                                % Create cell for feval function
    K1                = [cell2mat(q_now_tmp(1,end-2:end)),feval(EOM,q_now_tmp{:}).'];             % Calculate the second derivative at the start of the step
    q1_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K1(1,1:end-2));              % Create cell for feval function
    K2                = [cell2mat(q1_tmp(1,end-2:end)),feval(EOM,q1_tmp{:}).'];                   % Calculate the second derivative halfway the step
    q2_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K2(1:end-2));                % Refine value calculation with new found derivative
    K3                = [cell2mat(q2_tmp(1,end-2:end)),feval(EOM,q2_tmp{:}).'];                   % Calculate new derivative at the new refined location
    q3_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h)*K3(1:end-2));                    % Calculate state at end step with refined derivative
    K4                = [cell2mat(q3_tmp(1,end-2:end)),feval(EOM,q3_tmp{:}).'];                   % Calculate last second derivative                         % Take weighted sum of K1, K2, K3
    q_now_p           = (1/6)*(K1(end-4:end)+2*K2(end-4:end)+2*K3(end-4:end)+K4(end-4:end));      % Estimated current derivative
    q_next            = cell2mat(q_now_tmp) + (parms.h/6)*(K1(1:6)+2*K2(1:6)+2*K3(1:6)+K4(1:6));  % Perform euler intergration step
    
    % Save reaction forces and current derivative in state
    q(ii,end-4:end)   = q_now_p;
   
    % Save full state back in q array
    q         = [q;[q_next 0 0 0 0 0]];
    
    % Correct for intergration drift
    q_now_tmp = q(ii+1,:);
    [q_new,error] = gauss_newton(q_now_tmp,parms);  
    
    % Update the second derivative and the constraint forces
    q_new_tmp        = num2cell(q(ii,1:end-5),1);
    q_update          = feval(EOM,q_new_tmp{:}).';
    
    % Overwrite position coordinates
    q(ii+1,:)       = [q_new(1:6) q_update];
    
    % Create time array
t                   = [t;t(ii)+parms.h];          % Perform Gauss-Newton drift correction
ii                  = ii + 1;                                              % Append counter
t(ii)
q(ii,1)

% Calculate COM coordinates
% Calculate COM coordinates
x_tmp   = feval(parms.X_handle,q_new_tmp{1:3}).';
xdp_tmp = feval(parms.Xp_handle,q_new_tmp{:}).';

% Save x in state
x       = [x;x_tmp];
xdp     = [xdp;xdp_tmp];

end
end

%% Constraint calculation function
function [C,Cd] = constraint_calc(q,parms)

% Get needed angles out
q_now_tmp       = num2cell(q,1);

% Calculate the two needed constraints
C               = [parms.O4A*cos(q(2))+parms.O2A*cos(q(1))           ...
                   parms.O4B*sin(q(2))+parms.BC*sin(q(3))-parms.Yc];    
 
C_test          = feval(parms.C_handle,q_now_tmp{1:3}).';

% Calculate constraint derivative
Cd              = feval(parms.Cd_handle,q_now_tmp{1:3}).';

end

%% Speed correct function
function [q,error] = gauss_newton(q,parms)

% Get rid of the drift by solving a non-linear least square problem by
% means of the Gaus-Newton method
% Calculate the two needed constraints
[C,Cd] = constraint_calc(q,parms);

%% Guass-Newton position correction
n_iter          = 0;                                                                        % Set iteration counter                                                               % Get position data out

% Solve non-linear constraint least-square problem
while (max(abs(C)) > parms.tol)&& (n_iter < parms.nmax)
    q_tmp           = q(1:3);    
    n_iter = n_iter + 1;
    q_del  = Cd*inv(Cd.'*Cd)*-C.';
    q(1:3) = q_tmp+ q_del.';
    
    % Recalculate constraint
    [C,Cd]      = constraint_calc(q,parms);
end

% Calculate the corresponding speeds
q_tmp_vel          = q(4:6);
Dqd_n1             = -Cd*inv(Cd.'*Cd)*Cd.'*q_tmp_vel.';
q(4:6)             = q_tmp_vel + Dqd_n1.';

error = C;
end

%% Calculate (symbolic) Equations of Motion four our setup
function [qdp_handle,C_handle,Cd_handle,X_handle,Xd_handle] = EOM_calc(parms)

% Unpack symbolic variables from parms
x1              = parms.syms.x1;
y1              = parms.syms.y1;               
phi1            = parms.syms.phi1;
phi2            = parms.syms.phi2;

% Generalised state derivative
x1d             = parms.syms.x1d;
y1d             = parms.syms.y1d;
phi1d           = parms.syms.phi1d;
phi2d           = parms.syms.phi2d;

% Create generalized coordinate vectors
q               = [x1;y1;phi1;phi2];
qp              = [x1d;y1d;phi1d;phi2d];

% COM of the bodies expressed in generalised coordinates
x1              = x1;
y1              = y1;
phi1            = phi1;
x2              = parms.b*cos(phi1)+parms.d*cos(phi2);
y2              = parms.b*sin(phi1)+parms.d*sin(phi2);
phi2            = phi2;

% Create mass matrix (Segment 1 and 2)
parms.M         = diag([parms.m1,parms.m1,parms.J1,parms.m2,parms.m2,parms.J2]);

% Put in one state vector
x               = [x1;y1;phi1;x2;y2;phi2];

% Calculate the two needed constraints
C               = [parms.O4A*cos(phi4)+parms.O2A*cos(phi2)           ...
                   parms.O4B*sin(phi4)+parms.BC*sin(phi5)-parms.Yc];

% Compute the jacobian of state and constraints
Jx_q            = simplify(jacobian(x,q.'));
JC_q            = simplify(jacobian(C,q.'));

%% Calculate convective component
Jx_dq           = jacobian(Jx_q*qp,q);
JC_dq           = jacobian(JC_q*qp,q);

% Solve with virtual power
M_bar           = Jx_q.'*parms.M*Jx_q;

% Add forces F=[M2,F3_x,F3_y,M3,F4_x,F4_y,M4,F5_x,F5_y,M5,F6_x];
F               = [parms.T2, 0, -parms.m3*parms.g, 0, 0, -parms.m4*parms.g, 0, 0, -parms.m5*parms.g, 0, parms.F6_x];

% Create system of DAE
A = [M_bar JC_q.'; JC_q zeros(size(JC_q,1))];
B = [Jx_q.'*(F.'-parms.M*Jx_dq*qp); ...
    -JC_dq*qp];

% Calculate result expressed in generalized coordinates
qdp             = A\B;

% % Get result back in COM coordinates
% xdp     = simplify(jacobian(xp,qp.'))*qdp + simplify(jacobian(xp,q.'))*qp;

%% Convert to function handles
% xdp_handle       = matlabFunction(xdp);                                 % Create function handle of EOM in terms of COM positions
qdp_handle         = matlabFunction(simplify(qdp));                          % Create function handle of EOM in terms of generalised coordinates
% matlabFunction(qdp,'file',qdp_cal')

% Constraint function handle
C_handle        = matlabFunction(simplify(C));

% Constraint derivative function handle
Cd              = JC_q;
Cd_handle       = matlabFunction(simplify(Cd));

% Get back to COM coordinates
X_handle        = matlabFunction(simplify(x));
xp              = Jx_q*qp;
xdp             = simplify(jacobian(xp,qp))*qdp(1:3)+simplify(jacobian(xp,q))*qp(1:3);
Xd_handle       = matlabFunction(simplify(xdp));

end