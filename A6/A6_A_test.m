%% MBD_B: Assignment 6 - Double pendulum numerical intergration
%  Rick Staa (4511328)
%  Last edit: 02/05/2018
%  In this script I append the acceleration to the generalised state q so
%  q = [phi1 phi2 phi1p phi2p phi1dp phi2dp]. This was done to save space.
clear all; close all; clc;
fprintf('--- A6 ---\n');

%% Set up needed symbolic parameters
% Create needed symbolic variables
syms phi1 phi2 phi1p phi2p

% Put in parms struct for easy function handling
parms.syms.phi1               = phi1;
parms.syms.phi2               = phi2;
parms.syms.phi1p              = phi1p;
parms.syms.phi2p              = phi2p;

%% Intergration parameters
time                          = 3;                                         % Intergration time
parms.h                       = 0.001;                                     % Intergration step size

%% Model Parameters
% Segment 1
parms.L                       = 0.55;                                      % [parms.m]
parms.w                       = 0.05;                                      % [parms.m]
parms.t                       = 0.004;                                     % [parms.m]
parms.p                       = 1180;                                      % [kg/parms.m^3]
parms.m                       = parms.p * parms.w * parms.t * parms.L;     % [kg]
parms.I                       = (1/12) * parms.m * parms.L^2;              % [kg*parms.m^2]

% World parameters
parms.g                       = 9.81;                                      % [parms.m/s^2]

%% Initial state
q0                            = [0.5*pi 0.5*pi 0 0];

%% Derive equation of motion
[EOM_qdp] = EOM_calc(parms);                                       % Calculate symbolic equations of motion and put in parms struct

%% Calculate GLOBAL ERROR of the numerical intergration methods specified in the assignment
% 1). Euler (Euler)
% 2). Heun  (Heun)
% 3). Runge-Kutta 3th order (RK3)
% 4). Runge-Kutta 4th order (RK4)

tic
%% Euler intergration
% Calculate the error per step size for euler

% Loop h and calculate global error
n_range   = 6:1:25;
h_range   = time./(2.^n_range);
q_end_h_euler   = zeros(length(h_range),6);
for kk = 1:length(h_range)
    parms.h                       = h_range(kk);
    [t,q]                         = ODE_custom(EOM_qdp,time,q0,'euler',parms);
    % bar_animate(t,q,parms);                                                               % Animate Bar
    q_end_h_euler(kk,:)           = q(end,:);
end
glob_error                        = abs(q_end_h_euler(2:end,:)-q_end_h_euler(1:end-1,:));   % Calculate global error

% Create error plot
figure;
loglog(h_range(1:end-1),glob_error(:,1),'Color','red','LineWidth',1);hold on;
loglog(h_range(1:end-1),glob_error(:,2),'Color','blue','LineWidth',1);hold on;
line(xlim,[10e-6 10e-6],'Color',[1 0.6471 0],'LineStyle','--','LineWidth',1);
legend('\phi_1','\phi_2','Max error','Location', 'Best');
title('Euler numerical error');
xlabel('Global error [rad]')
ylabel('Intergration step size [s]')

%% Heun intergration
% Calculate the error per step size for heun
% Calculate the error per step size for euler

% Loop h and calculate global error
n_range   = 6:1:20;
h_range   = time./(2.^n_range);
q_end_h_heun   = zeros(length(h_range),6);
for kk = 1:length(h_range)
    parms.h                       = h_range(kk);
    [t,q]                         = ODE_custom(EOM_qdp,time,q0,'heun',parms);
    % bar_animate(t,q,parms);                                                             % Animate Bar
    q_end_h_heun(kk,:)            = q(end,:);
end
glob_error_heun                   = abs(q_end_h_heun(2:end,:)-q_end_h_heun(1:end-1,:));   % Calculate global error

% Create error plot
figure;
loglog(fliplr(h_range(1:end-1)),fliplr(glob_error_heun(:,1)'),'Color','red','LineWidth',1);hold on;
loglog(fliplr(h_range(1:end-1)),fliplr(glob_error_heun(:,2)'),'Color','blue','LineWidth',1);hold on;
line(xlim,[10e-6 10e-6],'Color',[1 0.6471 0],'LineStyle','--','LineWidth',1);
legend('\phi_1','\phi_2','Max error','Location', 'Best');
title('Euler numerical error');
xlabel('Global error [rad]')
ylabel('Intergration step size [s]')

%% Runge-Kutta 3th order (RK3)
% Calculate the error per step size for RK3

% Loop h and calculate global error
n_range   = 6:1:20;
h_range   = time./(2.^n_range);
q_end_h_RK3   = zeros(length(h_range),6);
for kk = 1:length(h_range)
    parms.h                       = h_range(kk);
    [t,q]                         = ODE_custom(EOM_qdp,time,q0,'RK3',parms);
    % bar_animate(t,q,parms);                                                             % Animate Bar
    q_end_h_RK3(kk,:)            = q(end,:);
end
glob_error_RK3                   = abs(q_end_h_RK3(2:end,:)-q_end_h_RK3(1:end-1,:));      % Calculate global error

% Create error plot
figure;
loglog(fliplr(h_range(1:end-1)),fliplr(glob_error_RK3(:,1)'),'Color','red','LineWidth',1);hold on;
loglog(fliplr(h_range(1:end-1)),fliplr(glob_error_RK3(:,2)'),'Color','blue','LineWidth',1);hold on;
line(xlim,[10e-6 10e-6],'Color',[1 0.6471 0],'LineStyle','--','LineWidth',1);
legend('\phi_1','\phi_2','Max error','Location', 'Best');
title('Euler numerical error');
xlabel('Global error [rad]')
ylabel('Intergration step size [s]')

%% Runge-Kutta 4th order (RK4)
% Calclate the error per step size for RK4

% Loop h and calculate global error
n_range   = 6:1:20;
h_range   = time./(2.^n_range);
q_end_h_RK4   = zeros(length(h_range),6);
for kk = 1:length(h_range)
    parms.h                       = h_range(kk);
    [t,q]                         = ODE_custom(EOM_qdp,time,q0,'RK4',parms);
    % bar_animate(t,q,parms);                                                             % Animate Bar
    q_end_h_RK4(kk,:)            = q(end,:);
end
glob_error_RK4                   = abs(q_end_h_RK4(2:end,:)-q_end_h_RK4(1:end-1,:));      % Calculate global error

% Create error plot
figure;
loglog(fliplr(h_range(1:end-1)),fliplr(glob_error_RK4(:,1)'),'Color','red','LineWidth',1);hold on;
loglog(fliplr(h_range(1:end-1)),fliplr(glob_error_RK4(:,2)'),'Color','blue','LineWidth',1);hold on;
line(xlim,[10e-6 10e-6],'Color',[1 0.6471 0],'LineStyle','--','LineWidth',1);
legend('\phi_1','\phi_2','Max error','Location', 'Best');
title('Euler numerical error');
xlabel('Global error [rad]')
ylabel('Intergration step size [s]')
toc;

%% Perform methods at maxstepsize
%% Euler method at step-size 1e-5
tic
parms.h = 1e-5;
[t_euler,q_euler]                         = ODE_custom(EOM_qdp,time,q0,'euler',parms);
toc

%% Heun method at step-size 1e-5
tic
parms.h = 1e-5;
[t_heun,q_heun]                         = ODE_custom(EOM_qdp,time,q0,'heun',parms);
toc

%% Runge-Kutta 3th order (RK4)
tic
parms.h = 1e-4;
[t_RK3,q_RK3]                         = ODE_custom(EOM_qdp,time,q0,'RK3',parms);
toc

%% Runge-Kutta 4th order (RK4)
tic
parms.h = 1e-4;
[t_RK4,q_RK4]                         = ODE_custom(EOM_qdp,time,q0,'RK4',parms);
toc

%% Calculate motion withODE functions
% ODE 23
tic
opt = odeset('AbsTol',1e-6,'RelTol',1e-6,'Stats','on');
[t23,q23] = ode23(@(t,q) ODE_func(t,q,EOM_qdp), [0 time], q0',opt);
t23_mean    = mean(diff(t23));                                             % Caculate mean step size
disp(t23_mean);
toc

% ODE 45
tic
opt = odeset('AbsTol',1e-6,'RelTol',1e-6,'Stats','on');
[t45,q45] = ode45(@(t,q) ODE_func(t,q,EOM_qdp), [0 time], q0',opt);
t45_mean    = mean(diff(t45));                                             % Caculate mean step size
disp(t45_mean);
toc

% ODE 113
tic
opt = odeset('AbsTol',1e-6,'RelTol',1e-6,'Stats','on');
[t113,q113] = ode113(@(t,q) ODE_func(t,q,EOM_qdp), [0 time], q0',opt);
t113_mean    = mean(diff(t45));                                            % Caculate mean step size
disp(t113_mean);
toc

%% FUNCTIONS

%% Bar animate
% This function creates a movie of the double pendulum. I did not find a
% way to do this with the real speed but it gives a nice impression of what
% is happening.
function bar_animate(t,q,parms)
figure;
h=plot(0,0,'MarkerSize',20,'Marker','.','LineWidth',2);
range=1.1*(parms.L+parms.L); axis([-range range -range range]); axis square;
set(gca,'nextplot','replacechildren');
a = tic;
for jj=1:length(q)-1
    if (ishandle(h)==1)
        tic
        phi1 = q(jj,1);
        phi2 = q(jj,2);
        Xcoord=[0,parms.L*cos(phi1),parms.L*cos(phi1)+parms.L*cos(phi2)];
        Ycoord=[0,parms.L*sin(phi1),parms.L*sin(phi1)+parms.L*sin(phi2)];
        set(h,'XData',Xcoord,'YData',Ycoord);
        b = toc(a); % check timer
        if b > (1/30)
            drawnow % update screen every 1/30 seconds
            a = tic; % reset timer after updating
        end
        %             pause(t(jj+1)-t(jj));                                          % Realtime
        toc;
    end
end
drawnow
end

%% ODE Function handle
function [qdp] = ODE_func(t,q,EOM_qdp)
q_now = num2cell(q',1);
qdp   = feval(EOM_qdp,q_now{:});
qdp   = [q(3);q(4);qdp];
end

%% Euler numerical intergration function
% This function calculates the motion of the system by means of a euler
% numerical intergration. This function takes as inputs the parameters of
% the system (parms), the EOM of the system (parms.EOM) and the initial
% state.
function [t,q] = ODE_custom(EOM,time,q0,method,parms)

% Initialise variables
t                   = (0:parms.h:time).';                                  % Create time array
q                   = zeros(length(t),6);                                  % Create empty state array
q(1,1:size(q0,2))   = q0;                                                  % Put initial state in array

% Caculate the motion for the full time by means of the 4 different
% numerical intergration methods
% 1). Euler
% 2). Heun
% 3). Runge-Kutta 3th order
% 4). Runge-Kutta 4th order
% See report for the Workings of each method.

% Euler method
switch method
    
    %% Euler method
    case 'euler'
        
        % Perform the full intergration with eulers method
        for ii = 1:(size(t,1)-1)
            q_now_tmp       = num2cell(q(ii,1:end-2),1);                   % Create cell for feval function
            qdp             = feval(EOM,q_now_tmp{:}).';                   % Calculate the second derivative of the generalised coordinates
            q(ii,end-1:end) = qdp;
            q(ii+1,1:end-2) = q(ii,1:end-2) + parms.h*q(ii,3:end);         % Perform euler intergration step
            
            % Calculate last acceleration
            if ii == (size(t,1)-1)
                q_next            = num2cell(q(ii+1,1:end-2),1);           % Create cell for feval function
                q(ii+1,end-1:end) = feval(EOM,q_next{:}).';                % Calculate the second derivative of the last step
            end
        end
        
        %% Heun method
    case 'heun'
        
        % Perform the full intergration with eulers method
        for ii = 1:(size(t,1)-1)
            % Step 1: Approximate the next state
            q_now           = [q(ii,1:end-2) 0 0];                              % Read out current states
            q_now_tmp       = num2cell(q_now,1);                                % Create cell for feval function
            qdp_now_tmp     = feval(EOM,q_now_tmp{1:end-2}).';                  % Calculate the second derivative of the generalised coordinates
            qdp_now         = [cell2mat(q_now_tmp(end-3:end-2)),qdp_now_tmp];   % Add first derivative
            q_now(end-3:end)= qdp_now;
            q_star          = q(ii,1:end-2) + parms.h*q_now(3:end);             % Make a approximation of the next state by means of a euler step
            
            % Step 2: Calculate the state derivative at next state
            q_star_tmp      = num2cell(q_star,1);                               % Create cell for feval function
            qdp_star_tmp    = feval(EOM,q_star_tmp{:}).';                       % Calculate the second derivative of the generalised coordinates
            qdp_star        = [cell2mat(q_star_tmp(end-1:end)),qdp_star_tmp];   % Add first derivative
            
            % Step3: Calculate the state at the next step using the mean
            % derivative.
            q(ii+1,1:end-2) = q(ii,1:end-2) + (parms.h*0.5)*(qdp_now+qdp_star);              % Calculate state of next step (I use both the approximated new velocity and acceleration)
            
            % Calculate last acceleration
            if ii == (size(t,1)-1)
                q_next_tmp        = num2cell(q(ii+1,1:end-2),1);                % Create cell for feval function
                q(ii+1,end-1:end) = feval(EOM,q_next_tmp{:}).';                 % Calculate the second derivative of the last step
            end
        end
        
        %% Runge-Kutta 3th order
    case 'RK3'
        for ii = 1:(size(t,1)-1)
            q_now_tmp         = num2cell(q(ii,1:end-2),1);                                                % Create cell for feval function
            K1                = [cell2mat(q_now_tmp(1,end-1:end)),feval(EOM,q_now_tmp{:}).'];             % Calculate the second derivative at the start of the step
            q1_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K1);                         % Create cell for feval function
            K2                = [cell2mat(q1_tmp(1,end-1:end)),feval(EOM,q1_tmp{:}).'];                   % Calculate the second derivative halfway the step
            q2_tmp            = num2cell(cell2mat(q_now_tmp) + ((parms.h*0.75))*K2);                      % Refine value calculation with new found derivative
            K3                = [cell2mat(q2_tmp(1,end-1:end)),feval(EOM,q2_tmp{:}).'];                   % Calculate new derivative at the new refined location
            q(ii,end-1:end)   = (1/9)*(2*K1(3:4)+3*K2(1:2)+4*K3(3:4));                                                   % Take weighted sum of K1, K2, K3
            q(ii+1,1:end-2)   = cell2mat(q_now_tmp) + (parms.h/9)*(2*K1+3*K2+4*K3);                       % Perform euler intergration step
            
            % Calculate last acceleration
            if ii == (size(t,1)-1)
                q_now_tmp         = num2cell(q(ii+1,1:end-2),1);                                              % Create cell for feval function
                K1                = [cell2mat(q_now_tmp(1,end-1:end)),feval(EOM,q_now_tmp{:}).'];             % Calculate the second derivative at the start of the step
                q1_tmp            = num2cell(cell2mat(q_now_tmp) + (parms.h*0.5)*K1);                         % Create cell for feval function
                K2                = [cell2mat(q1_tmp(1,end-1:end)),feval(EOM,q1_tmp{:}).'];                   % Calculate the second derivative halfway the step
                q2_tmp            = num2cell(cell2mat(q_now_tmp) + ((parms.h*0.75))*K2);                      % Refine value calculation with new found derivative
                K3                = [cell2mat(q2_tmp(1,end-1:end)),feval(EOM,q2_tmp{:}).'];                   % Calculate new derivative at the new refined location
                q(ii+1,end-1:end) = (1/9)*(2*K1(3:4)+3*K2(3:4)+4*K3(3:4));                                    % Take weighted sum of K1, K2, K3                                               % Take weighted sum of K1, K2, K3
            end
        end
        
        %% Runge-Kutta 4th order
    case 'RK4'
        for ii = 1:(size(t,1)-1)
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
end

%% Calculate (symbolic) Equations of Motion four our setup
function [qdp_handle] = EOM_calc(parms)

syms L
syms m
syms g
syms I
% Unpack symbolic variables from varargin
phi1            = parms.syms.phi1;
phi2            = parms.syms.phi2;
phi1p           = parms.syms.phi1p;
phi2p           = parms.syms.phi2p;

% Create generalized coordinate vectors
q               = [phi1; phi2];
qp              = [phi1p; phi2p];

% COM of the bodies expressed in generalised coordinates
x1              = (L/2)*cos(phi1);
y1              = (L/2)*sin(phi1);
x2              = x1 + (L/2)*cos(phi1) + (L/2) * cos(phi2);
y2              = y1 + (L/2)*sin(phi1) + (L/2) * sin(phi2);

% Calculate derivative of COM expressed in generalised coordinates (We need this for the energy equation)
x               = [x1;y1;phi1;x2;y2;phi2];
Jx_q            = simplify(jacobian(x,q));
xp              = Jx_q*qp;

%% Compute energies
T               = 0.5*xp.'*diag([m;m;I;m;m;I])*xp;           % Kinetic energy
V               = -([m*g 0 0 m*g 0 0]*x);                    % Potential energy

%% Calculate the terms of the jacobian
Q                = 0;                           % Non-conservative forces

% Partial derivatives of Kinetic energy
T_q             = simplify(jacobian(T,q));
T_qp            = simplify(jacobian(T,qp));
T_qpqp          = simplify(jacobian(T_qp,qp));
T_qpq           = simplify(jacobian(T_qp,q));

% Partial derivatives of Potential energy
V_q             = simplify(jacobian(V,q));
V_qp            = simplify(jacobian(V,qp));
V_qpqp          = simplify(jacobian(V_qp,qp));

% Make matrix vector product
M                = T_qpqp;
F                = Q + T_q.' - V_q.' - T_qpq*qp;

% Solve Mqdp=F to get the accelerations
qdp              = M\F;

%% Get back to COM coordinates
% xdp              = simplify(jacobian(xp,qp))*qdp+simplify(jacobian(xp,q))*qp;

%% Convert to function handles
% xdp_handle       = matlabFunction(xdp);                                  % Create function handle of EOM in terms of COM positions
qdp_handle       = matlabFunction(qdp);                                    % Create function handle of EOM in terms of generalised coordinates
% matlabFunction(qdp,'file',qdp_cal')
end