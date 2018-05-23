%% MBD_B: Assignment 7 - Quick return mechanism
%  Rick Staa (4511328)
%  Last edit: 09/05/2018
clear all; close all; clc;
fprintf('--- A7 ---\n');

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
parms.h                     = 0.01;                                     % Intergration step size

%% Model Parameters
% Lengths and distances
parms.O2A                   = 0.2;                                      % Length segment 2 [m]
parms.O4B                   = 0.7;                                      % Length segment 4 [m]
parms.BC                    = 0.6;                                      % Length segment 5 [m]
parms.O4O2                  = 0.3;                                      % Distance between joint 4 and joint 2 [m]
parms.O4G4                  = 0.4;                                      % Distance bewteen COM4 and joint 4 [m]
parms.BG5                   = 0.3;                                      % Distance joint 5 and COM 5 [m]
parms.Yc                    = 0.9;                                      % Height joint C (COM body 6) [m]
parms.O4A                   = sqrt(parms.O2A^2+parms.O4O2^2);           % Distance between joint 4 and joint 3 [m]

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
parms.g                     = 9.81;                                     % [parms.m/s^2]

% Forces
parms.F6_x                  = 1000;                                     % x force on body 6 [N]
parms.T2                    = 0;                                        % Torque around joint 6 [Nm]

%% Calculate Initial states
phi2_init                   = 0;
phi4_init                   = atan2(parms.O4O2,parms.O2A);
phi5_init                   = pi-asin((parms.Yc-parms.O4B*sin(phi4_init))/parms.BC);
phi2d_init                  = (150*pi)/60;
phi4d_init                  = 0;
phi5d_init                  = 0;
q0                          = [phi2_init phi4_init phi5_init phi2d_init phi4d_init phi5d_init];

%% Derive equation of motion
[EOM_qdp,C_handle]          = EOM_calc(parms);                         % Calculate symbolic equations of motion and put in parms struct
parms.C_handle              = C_handle;

%% Calculate motion RK4

%% Runge-Kutta 4th order (RK4)
tic
[t_RK4,q_RK4]                         = RK4_custom(EOM_qdp,time,q0,parms);
toc

% %% Calculate motion with ODE 113
% tic
% opt = odeset('AbsTol',1e-6,'RelTol',1e-6,'Stats','on');
% [t113,q113] = ode113(@(t,q) ODE_func(t,q,EOM_qdp), [0 time], q0',opt);
% t113_mean    = mean(diff(t45));                                            % Caculate mean step size
% disp(t113_mean);
% toc

%% FUNCTIONS

%% Euler numerical intergration function
% This function calculates the motion of the system by means of a euler
% numerical intergration. This function takes as inputs the parameters of
% the system (parms), the EOM of the system (parms.EOM) and the initial
% state.
function [t,q] = RK4_custom(EOM,time,q0,parms)

% Initialise variables
t                   = (0:parms.h:time).';                                  % Create time array
q                   = zeros(length(t),11);                                  % Create empty state array
q(1,1:size(q0,2))   = q0;                                                  % Put initial state in array

% Caculate the motion for the full simulation time by means of a
% Runge-Kutta4 method

% Perform intergration for 0 till time
for ii = 1:(size(t,1)-1)
    
    % Calculate the next state by means of a RK4 method
    q_now_tmp         = num2cell(q(ii,1:end-2),1);                                                % Create cell for feval function
    K1                = [cell2mat(q_now_tmp(1,end-1:end)),feval(EOM,q_now_tmp{:}).'];             % Calculate the second derivative at the start of the step
    q1_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K1);                         % Create cell for feval function
    K2                = [cell2mat(q1_tmp(1,end-1:end)),feval(EOM,q1_tmp{:}).'];                   % Calculate the second derivative halfway the step
    q2_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K2);                         % Refine value calculation with new found derivative
    K3                = [cell2mat(q2_tmp(1,end-1:end)),feval(EOM,q2_tmp{:}).'];                   % Calculate new derivative at the new refined location
    q3_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h)*K3);                             % Calculate state at end step with refined derivative
    K4                = [cell2mat(q3_tmp(1,end-1:end)),feval(EOM,q3_tmp{:}).'];                   % Calculate last second derivative
    q(ii,end-1:end)   = (1/6)*(K1(3:4)+2*K2(3:4)+2*K3(3:4)+K4(3:4));                              % Take weighted sum of K1, K2, K3
    q(ii+1,1:end-2)   = cell2mat(q_now_tmp) + (parms.h/6)*(K1+2*K2+2*K3+K4);                      % Perform euler intergration step
    
    % Correct for intergration drift
    [q,error] = gauss_newton(q,parms)
    
    % Calculate last acceleration
    if ii == (size(t,1)-1)
        q_now_tmp         = num2cell(q(ii+1,1:end-2),1);                                              % Create cell for feval function
        K1                = [cell2mat(q_now_tmp(1,end-1:end)),feval(EOM,q_now_tmp{:}).'];             % Calculate the second derivative at the start of the step
        q1_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K1);                         % Create cell for feval function
        K2                = [cell2mat(q1_tmp(1,end-1:end)),feval(EOM,q1_tmp{:}).'];                   % Calculate the second derivative halfway the step
        q2_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K2);                         % Refine value calculation with new found derivative
        K3                = [cell2mat(q2_tmp(1,end-1:end)),feval(EOM,q2_tmp{:}).'];                   % Calculate new derivative at the new refined location
        q3_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h)*K3);                             % Calculate state at end step with refined derivative
        K4                = [cell2mat(q3_tmp(1,end-1:end)),feval(EOM,q3_tmp{:}).'];                   % Calculate last second derivative
        q(ii,end-1:end)   = (1/6)*(K1(3:4)+2*K2(3:4)+2*K3(3:4)+K4(3:4));                              % Take weighted sum of K1, K2, K3
    end
end
end

%% Speed correct function
function [q,error] = gauss_newton(q,parms)

% Get rid of the drift by solving a non-linear least square problem by
% means of the Gaus-Newton method
% Calculate the two needed constraints
C               = [parms.O4A*cos(phi4)+parms.O2A*cos(phi2)           ...
    parms.O4B*sin(phi4)+parms.BC*sin(phi5)-parms.Yc];
end

%% ODE Function handle
function [qdp] = ODE_func(t,q,EOM_qdp)
q_now = num2cell(q',1);
qdp   = feval(EOM_qdp,q_now{:});
qdp   = [q(3);q(4);qdp];
end

%% Calculate (symbolic) Equations of Motion four our setup
function [qdp_handle,C_handle] = EOM_calc(parms)

% Unpack symbolic variables from varargin
phi2            = parms.syms.phi2;
phi4            = parms.syms.phi4;
phi5            = parms.syms.phi5;
phi2d           = parms.syms.phi2d;
phi4d           = parms.syms.phi4d;
phi5d           = parms.syms.phi5d;

% Create generalized coordinate vectors
q               = [phi2;phi4;phi5];
qp              = [phi2d;phi4d;phi5d];

% COM of the bodies expressed in generalised coordinates
% x2              = 0;
% y2              = parms.O4O2;
x3              = parms.O2A*cos(phi2);
y3              = parms.O4O2+parms.O2A*sin(phi2);
x4              = parms.O4G4*cos(phi4);
y4              = parms.O4G4*sin(phi4);
x5              = parms.O4B*cos(phi4)+parms.BG5*cos(phi5);
y5              = parms.O4B*sin(phi4)+parms.BG5*sin(phi5);
x6              = parms.O4B*cos(phi4)+parms.BC*cos(phi5);
% y6             = parms.O4B*sin(phi4)+parms.BC*sin(phi5);

% Create mass matrix
% x2 = 0, y2 = 0 and y5 =0 also no moments around slider 3 and 6
parms.M         = diag([parms.J2,parms.m3,parms.m3,parms.J3,parms.m4,parms.m4,parms.J4,parms.m5,parms.m5,parms.J5,parms.m6]);

% Put in one state vector
x               = [phi2;x3;y3;phi4;x4;y4;phi4;x5;y5;phi5;x6];

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
qdp_handle      = matlabFunction(simplify(qdp));                                    % Create function handle of EOM in terms of generalised coordinates
% matlabFunction(qdp,'file',qdp_cal')

% Constraint function handle
C_handle        = matlabFunction(simplify(C));

end