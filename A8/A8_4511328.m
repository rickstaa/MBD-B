%% MBD_B: Assignment 8 - EzyRoller
%  Rick Staa (4511328)
%  Last edit: 29/05/2018

%% - - Pre processing operations --
clear all; close all; clc;
fprintf('--- A8 ---\n');
animate_bool = 0;                           % Set on 1 if you want to see an animation

% Set up needed symbolic parameters
syms x1 y1 phi1 x2 y2 phi2 x1d y1d phi1d x2d y2d phi2d t

% State
parms.syms.x1               = x1;
parms.syms.y1               = y1;
parms.syms.phi1             = phi1;
parms.syms.x2               = x2;
parms.syms.y2               = y2;
parms.syms.phi2             = phi2;
parms.syms.t                = t;

% State derivative
parms.syms.x1d               = x1d;
parms.syms.y1d               = y1d;
parms.syms.phi1d             = phi1d;
parms.syms.x2d               = x2d;
parms.syms.y2d               = y2d;
parms.syms.phi2d             = phi2d;

%% -- Set model/simulation parameters and initial states --
%% Intergration parameters
sim_time                    = 100;                                      % Intergration time
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

% Create mass matrix (Segment 1 and 2)
parms.M                     = diag([parms.m1,parms.m1,parms.J1,parms.m2,parms.m2,parms.J2]);

% Torque and force variables (See assignment)
parms.M0                    = 0.1;
parms.omega                 = pi;

%% World parameters
% Gravity
parms.g                     = 9.81;                                     % [parms.m/s^2]

% %% states for Question 1
% x1_0                        = parms.a;
% y1_0                        = 0;
% phi1_0                      = 0;
% x2_0                        = parms.a+parms.b;
% y2_0                        = parms.d;
% phi2_0                      = pi/2;
%
% % Phi1d
% x1d_0                       = 1;
% y1d_0                       = 0;
% phi1d_0                     = 0;
% x2d_0                       = 0;
% y2d_0                       = 1;
% phi2d_0                     = 0;
%
% % Set forces
% F                           = [0 0 0 0 0 0].';                        % No torque applied
% parms.F                     = F;
% x0                          = [x1_0 y1_0 phi1_0 x2_0 y2_0 phi2_0 x1d_0 y1d_0 phi1d_0 x2d_0 y2d_0 phi2d_0];

%% States Question 2
% In this the generalised coordinates x1_init and y1_init are assumed to be
% defined so that wheel 1 is in the origin.
phi1_0                      = 0;                                        % Angle of first body with horizontal
phi2_0                      = pi;                                       % Angle of second body with horizontal

% Calculate other dependent initial positions and angles
x1_0                        = parms.a*cos(phi1_0);
y1_0                        = parms.b*sin(phi1_0);
x2_0                        = (parms.a+parms.b)*cos(phi1_0)+parms.d*cos(phi2_0);
y2_0                        = (parms.a+parms.b)*sin(phi1_0)+parms.d*sin(phi2_0);

% Velocity initital states (Make sure that the are admissable)
% Phi1d
x1d_0                       = 0;
y1d_0                       = 0;
phi1d_0                     = 0;
x2d_0                       = 0;
y2d_0                       = 0;
phi2d_0                     = 0;

% Create full state for optimization
x0                          = [x1_0 y1_0 phi1_0 x2_0 y2_0 phi2_0 x1d_0 y1d_0 phi1d_0 x2d_0 y2d_0 phi2d_0];

%% Set Forces and torques
% F=[F1_x,F1_y,M1,F2_x,F2_y,M2];
F                           = [0 0 -parms.M0*cos(parms.omega*t) 0 0 parms.M0*cos(parms.omega*t)].';                            % Torque applied

% Store F in function
parms.F                     = F;

%% -- Derive equation of motion --
%% Calculate EOM by means of Newton-Euler equations
EOM_calc(parms);                                                        % Calculate symbolic equations of motion and put in parms struct

%% -- Perform simulation --
%% Calculate movement by mean sof a Runge-Kuta 4th order intergration method
tic
[t,x]                       = RK4_custom(x0,sim_time,parms);
toc

%% -- Post Processing --
%% Calculate com velocities
% xd                        = diff(x)/parms.h;
xdd                         = state_deriv(x,parms);

%% Calculate position of point A B and C
[A,B,C]                     = point_calc(x,parms);

%% Calculate kinetic energy and torque wo[ekin] = ekin_calc(x,parms);
[ekin]                      = ekin_calc(x,parms);
[tw]                        = tw_calc(x,parms);

%% -- ANIMATE --
if animate_bool == 1
    % Adapted from A. Schwab's animation code
    
    % Rename data
    X1 = x(:,1); Y1 = x(:,2); P1 = x(:,3);
    DX1 = x(:,7); DY1 = x(:,8); DP1 = x(:,9);
    X2 = x(:,4); Y2 = x(:,5); P2 = x(:,6);
    DX2 = x(:,10); DY2 = x(:,11); DP2 = x(:,12);
    
    % Rename Points
    XA = A(:,1); YA = A(:,2);
    XB = B(:,1); YB = B(:,2);
    XC = C(:,1); YC = C(:,2);
    
    % Create figure
    figure
    plot(X1,Y1)
    hold on
    plot(XA,YA)
    hold on
    plot(X2,Y2)
    hold on
    plot(XC,YC)
    grid on
    set(gca,'fontsize',16)
    title('Animation EzyRoller')
    axis([min(X1)-parms.a max(X1)+parms.a min(Y1)-parms.a max(Y1)+parms.a]);
    axis equal
    l = plot([X1(1) XA(1)],[Y1(1) YA(1)]);
    k = plot([X2(1) XC(1)],[Y2(1) YC(1)]);
    j = plot([X1(1) XB(1)],[Y1(1) YB(1)]);
    m = plot([XB(1) X2(1)],[YB(1) Y2(1)]);
    set(l,'LineWidth',5);
    set(l,'Color','K')
    set(k,'LineWidth',5);
    set(k,'Color','C')
    set(j,'LineWidth',5);
    set(j,'Color','K')
    set(m,'LineWidth',5);
    set(m,'Color','C')
    nstep = length(t);
    nskip = 10;
    for istep = 2:nskip:nstep
        set(l,'XData',[X1(istep) XA(istep)])
        set(l,'YData',[Y1(istep) YA(istep)])
        set(k,'XData',[X2(istep) XC(istep)])
        set(k,'YData',[Y2(istep) YC(istep)])
        set(j,'XData',[X1(istep) XB(istep)])
        set(j,'YData',[Y1(istep) YB(istep)])
        set(m,'XData',[XB(istep) X2(istep)])
        set(m,'YData',[YB(istep) Y2(istep)])
        drawnow
        pause(1e-10)
    end
    
end

%% - - Create plots - -
%% Plot path of points on the robot
figure;
plot(A(:,1),A(:,2),x(:,1),x(:,2),B(:,1),B(:,2),x(:,4),x(:,5),C(:,1),C(:,2),'linewidth',1.5);
set(gca,'fontsize',18);
title('Travelled path of points on the EzyRoller');
xlabel('x position [m]');
ylabel('y position [m]');
legend('Point A (Wheel 1)','COM 1','Point B (hindge)','COM 2','Point C (wheel 2)','Location', 'Best');

%% Plot linear velocities COM's
figure;
subplot(2,1,1);
plot(t,x(:,7),'b',t,x(:,10),'r','Linewidth',1.5);
title("x velocities of COM's");
xlabel('t [s]');
ylabel('velocity [m/s]');
legend('x velocity (COM 1)','x velocity (COM 2)','Location', 'Best');
subplot(2,1,2);
plot(t,x(:,8),'b',t,x(:,11),'r','Linewidth',1.5);
title("y velocities of COM's");
xlabel('t [s]');
ylabel('velocity [m/s]');
legend('y velocity (COM 1)','y velocity (COM 2)','Location', 'Best');

%% Plot linear magnitude velocities COM's
% Calculate velocity magnitudes
v_com1 = sqrt(x(:,7).^2+x(:,8).^2);
v_com2 = sqrt(x(:,10).^2+x(:,11).^2);

% Plot figure
figure;
plot(t,v_com1,'b',t,v_com2,'r','Linewidth',1.5);
title("velocitie mag of COM's");
xlabel('t [s]');
ylabel('velocity mag [m/s]');
legend('velocity (COM 1)','velocity (COM 2)','Location', 'Best');

%% Plot angular velocities
figure;
plot(t,x(:,9),'b',t,x(:,12),'r','Linewidth',1.5);
title("angular velocities of COM's");
xlabel('t [s]');
ylabel('angular velocity [rad/s]');
legend('angular velocity (COM 1)','angular velocity (COM 2)','Location', 'Best');

%% Plot linear and angular accelerations COM's
figure;
subplot(2,1,1);
plot(t,xdd(:,7),'b',t,xdd(:,10),'r','Linewidth',1.5);
title("x accelerations of COM's");
xlabel('t [s]');
ylabel('accelleration [m/s^2]');
legend('x accelleration (COM 1)','x celleration (COM 2)','Location', 'Best');
subplot(2,1,2);
plot(t,xdd(:,8),'b',t,xdd(:,11),'r','Linewidth',1.5);
title("y accellerations of COM's");
xlabel('t [s]');
ylabel('Accelleration [m/s^2]');
legend('y accelleration (COM 1)','y accelleration (COM 2)','Location', 'Best');

%% Plot angular accelerations-
figure;
plot(t,xdd(:,9),'b',t,xdd(:,12),'r','Linewidth',1.5);
title("Angular velocities of COM's");
xlabel('t [s]');
ylabel('Angular acceleration [rad/s^2]');
legend('Angular acceleration (COM 1)','Angular acceleration (COM 2)','Location', 'Best');

%% Plot reaction forces
figure;
plot(t,x(:,13:end),'Linewidth',1.5);
title("Reaction forces in the constraints");
xlabel('t [s]');
ylabel('Reaction Force [N]');
legend('X reaction force in joint B (FB_x)','Y reaction force in joint B (FB_y)','Wheel A friction force (no slip)','Wheel C friction force (no slip)','Location', 'Best');

%% Plot kinetic energy
figure;
plot(t,ekin(:,1),'-b',t,ekin(:,2),'-r',t,ekin(:,3),'-g','Linewidth',1.5)
set(gca,'fontsize',18);
title("Kinetic energy of the COM's");
xlabel('t [s]');
ylabel('Kinetic energy[Joule]');
legend('Kinetic energy (COM 1)','kinetic energy (COM 2)','Kinetic energy system','Location', 'Best');

%% Plot Kinetic energy plus torque energy
figure;
plot(t,ekin(:,1),'-m',t,ekin(:,2),'-c',t,ekin(:,3),'-b',t,tw,'--r','Linewidth',1.5)
set(gca,'fontsize',18);
title("Kinetic energy of the COM's and Torque Work");
xlabel('t [s]');
ylabel('Energy [Joule]');
legend('Kinetic energy (COM 1)','kinetic energy (COM 2)','Kinetic energy system','Input torque work','Location', 'Best');

%% FUNCTIONS

%% Post processing functions
% These functions are used to calculate quantaties that are not calculated
% during the simulation. This regards quantaties which are not state
% variables

% Calculate second derivative
function [xdd] = state_deriv(x,parms)

% preallocate memory for xdd vector
xdd         = zeros(size(x,1),12);

% Create time vector
time        = 0:parms.h:((parms.h*size(x,1))-parms.h);

% Loop through states
for ii = 1:size(x,1)
    % Set time
    t = time(ii);
    
    % Calculate xdd
    x_now_tmp    = x(ii,1:end-4);
    x_now_input  = num2cell([x(ii,[3 6 7:12]),t],1);
    xdd_tmp      = subs_xdd(x_now_input{:}).';
    xdd(ii,:)    = [x_now_tmp(7:12),xdd_tmp(1:6)];
end
end

% Calculation points on EzyRoller
function [A,B,C] = point_calc(x,parms)

%% Calculate Point A, B, C out of the state
A_x             = x(:,1)-parms.a*cos(x(:,3));
A_y             = x(:,2)-parms.a*sin(x(:,3));
B_x             = x(:,1)+parms.b*cos(x(:,3));
B_y             = x(:,2)+parms.b*sin(x(:,3));
C_x             = x(:,4)+parms.c*cos(x(:,6));
C_y             = x(:,5)+parms.c*sin(x(:,6));

% Put them in their corresponding vector
A               = [A_x A_y];
B               = [B_x B_y];
C               = [C_x C_y];

end

% Calculate kinetic energy of COM's
function [ekin] = ekin_calc(x,parms)

% preallocate memory for ekin vector
ekin            = zeros(size(x,1),1);

% Loop through states
% State is x = [x1 y1 phi1 x2 y2 phi2 x1p y1p phi1p x2p y2p phi2p
for ii = 1:size(x,1)
    ekin(ii,1) = 0.5*x(ii,7:9)*parms.M(1:3,1:3)*x(ii,7:9).';
    ekin(ii,2) = 0.5*x(ii,10:12)*parms.M(4:6,4:6)*x(ii,10:12).';
    ekin(ii,3) = 0.5*x(ii,7:12)*parms.M*x(ii,7:12).';
end
end

% Calculate kinetic energy of COM's
function [tw] = tw_calc(x,parms)

% Calculate the applied torque for the whole movement
% preallocate memory for xdd vector
tw              = zeros(size(x,1),1);

% Create time vector
time            = 0:parms.h:((parms.h*size(x,1))-parms.h);

% Create W vector
for ii = (2:size(x,1))
    tw(ii)     = tw(ii-1) + sum((subs_F(time(ii))).'.*(x(ii,1:6)-x((ii-1),1:6)));
end
end

%% Runge-Kuta numerical intergration function
% This function calculates the motion of the system by means of a
% Runge-Kuta numerical intergration. This function takes as inputs the
% parameters of the system (parms), the EOM of the system (parms.EOM)
% and the initial state.
function [time,x] = RK4_custom(x0,sim_time,parms)

% Initialise variables
time                = (0:parms.h:sim_time).';                                  % Create time array
x                   = zeros(length(time),16);                                 % Create empty state array
x(1,1:length(x0))   = x0;                                                  % Put initial state in array

% Caculate the motion for the full simulation time by means of a
% Runge-Kutta4 method

% Perform intergration till end of set time
for ii = 1:(size(time,1)-1)
    
    % Add time constant
    t = time(ii);
    
    % Perform RK 4
    x_now_tmp           = x(ii,1:end-4);                                                                    % Create cell for subs function function
    x_input             = num2cell([x(ii,[3 6 7:12]),t],1);                                                 % Add time to state
    K1                  = [x_now_tmp(1,end-5:end),subs_xdd(x_input{:}).'];                                  % Calculate the second derivative at the start of the step
    x1_tmp              = x_now_tmp + (parms.h*0.5)*K1(1:end-4);                                            % Create cell for subs function function
    x1_input            = num2cell([x1_tmp([3 6 7:12]),t],1);                                               % Add time to state
    K2                  = [x1_tmp(1,end-5:end),subs_xdd(x1_input{:}).'];                                    % Calculate the second derivative halfway the step
    x2_tmp              = x_now_tmp + (parms.h*0.5)*K2(1:end-4);                                            % Refine value calculation with new found derivative
    x2_input            = num2cell([x2_tmp([3 6 7:12]),t],1);                                               % Add time to state
    K3                  = [x2_tmp(1,end-5:end),subs_xdd(x2_input{:}).'];                                    % Calculate new derivative at the new refined location
    x3_tmp              = x_now_tmp + (parms.h)*K3(1:end-4);                                                % Calculate state at end step with refined derivative
    x3_input            = num2cell([x3_tmp([3 6 7:12]),t],1);                                               % Add time to state
    K4                  = [x3_tmp(1,end-5:end),subs_xdd(x3_input{:}).'];                                    % Calculate last second derivative
    x(ii,end-3:end)     = (1/6)*(K1(end-3:end)+2*K2(end-3:end)+2*K3(end-3:end)+K4(end-3:end));              % Take weighted sum of K1, K2, K3
    x(ii+1,1:end-4)     = x_now_tmp + (parms.h/6)*(K1(1:end-4)+2*K2(1:end-4)+2*K3(1:end-4)+K4(1:end-4));    % Perform euler intergration step
    
    % Calculate last acceleration
    if ii == (size(time,1)-1)
        x_now_tmp           = x(ii+1,1:end-4);                                                                  % Create cell for subs function function
        x_input             = num2cell([x(ii+1,[3 6 7:12]),t],1);                                               % Add time to state
        K1                  = [x_now_tmp(1,end-5:end),subs_xdd(x_input{:}).'];                                  % Calculate the second derivative at the start of the step
        x1_tmp              = x_now_tmp + (parms.h*0.5)*K1(1:end-4);                                            % Create cell for subs function function
        x1_input            = num2cell([x1_tmp([3 6 7:12]),t],1);                                               % Add time to state
        K2                  = [x1_tmp(1,end-5:end),subs_xdd(x1_input{:}).'];                                    % Calculate the second derivative halfway the step
        x2_tmp              = x_now_tmp + (parms.h*0.5)*K2(1:end-4);                                            % Refine value calculation with new found derivative
        x2_input            = num2cell([x2_tmp([3 6 7:12]),t],1);                                               % Add time to state
        K3                  = [x2_tmp(1,end-5:end),subs_xdd(x2_input{:}).'];                                    % Calculate new derivative at the new refined location
        x3_tmp              = x_now_tmp + (parms.h)*K3(1:end-4);                                                % Calculate state at end step with refined derivative
        x3_input            = num2cell([x3_tmp([3 6 7:12]),t],1);                                               % Add time to state
        K4                  = [x3_tmp(1,end-5:end),subs_xdd(x3_input{:}).'];                                    % Calculate last second derivative
        x(ii+1,end-3:end)     = (1/6)*(K1(end-3:end)+2*K2(end-3:end)+2*K3(end-3:end)+K4(end-3:end));            % Take weighted sum of K1, K2, K3
    end
    
    % Correct for intergration drift
    x_now_tmp = x(ii+1,:);
    [x_new,~] = gauss_newton(x_now_tmp,parms);
    
    % Update the constraint forces
    x_new_input       = num2cell([x(ii,[3 6 7:12]),t],1);
    x_update          = subs_xdd(x_new_input{:}).';
    
    % Overwrite position coordinates
    x(ii+1,:)       = [x_new(1:end-4) x_update(end-3:end)];
    
end
end

%% Constraint calculation function
function [C,Cd,D,Dd] = constraint_calc(x,parms)

% Get needed angles out
x_now_tmp       = num2cell(x,1);

%% Calculate position constraint
C               = subs_C(x_now_tmp{1:6}).';

% Calculate constraint derivative
Cd              = subs_Cd(x_now_tmp{[3 6]}).';

%% Calculate velocity constraint
D               = subs_D(x_now_tmp{[3 6:12]}).';

% Calculate velocity constraint derivative
Dd              = subs_Dd(x_now_tmp{[3 6]}).';
end

%% Speed correct function
function [x,error] = gauss_newton(x,parms)

% Get rid of the drift by solving a non-linear least square problem by
% means of the Gaus-Newton method
% Calculate the two needed constraints
[C,Cd,~,~] = constraint_calc(x,parms);

%% Guass-Newton position constraint correction
n_iter          = 0;                                                                        % Set iteration counter                                                               % Get position data out

% Solve non-linear constraint least-square problem
while (max(abs(C)) > parms.tol)&& (n_iter < parms.nmax)
    x_tmp           = x(1:6);
    n_iter = n_iter + 1;
    x_del  = Cd*inv(Cd.'*Cd)*-C.';
    x(1:6) = x_tmp+ x_del.';
    
    % Recalculate constraint
    [C,Cd,~,~]      = constraint_calc(x,parms);
end

% % Calculate the corresponding speeds
% x_tmp_vel          = x(7:12);
% Dxd_n1             = -Cd*inv(Cd.'*Cd)*Cd.'*x_tmp_vel.';
% x(7:12)            = x_tmp_vel + Dxd_n1.';
%

%% Gaus-newton velocity constraint correction
n_iter          = 0;                                                                        % Set iteration counter                                                               % Get position data out

% % Calculate the two needed constraints
% [~,~,D,Dd] = constraint_calc(x,parms);

% % Solve non-linear constraint least-square problem
% while (max(abs(D)) > parms.tol)&& (n_iter < parms.nmax)
%     x_tmp           = x(7:12);
%     n_iter = n_iter + 1;
%     x_del  = Dd*inv(Dd.'*Dd)*-D.';
%     x(7:12) = x_tmp+ x_del.';
%
%     % Recalculate constraint
%     [~,~,D,Dd]      = constraint_calc(x,parms);
% end


% Calculate constraints
[~,Cd,D,Dd]        = constraint_calc(x,parms);
Sd                 = [Cd Dd];

% Calculate new velocities
x_tmp_vel          = x(7:12);
Dxd_n1             = -Sd*inv(Sd.'*Sd)*Sd.'*x_tmp_vel.';
x(7:12)            = x_tmp_vel + Dxd_n1.';

%% Recalculate error
[C,~,D,~]        = constraint_calc(x,parms);
C_error = C;
D_error = D;

% Store full error
error = [C_error D_error];
end

%% Calculate (symbolic) Equations of Motion four our setup
function EOM_calc(parms)

%% -- The code between this lines is done to obtain the latex formulas --
% % Create model parameters in symbolic form
% syms a b c d m1 m2 J1 J2 g;

% Overwrite with real values if you don't want the full symbolic expresion
a               = parms.a;
b               = parms.b;
c               = parms.c;
d               = parms.d;
m1              = parms.m1;
m2              = parms.m2;
J1              = parms.J1;
J2              = parms.J2;
g               = parms.g;

%% -- The code between this lines is done to create the latex formulas --

% Unpack symbolic variables from parms
x1              = parms.syms.x1;
y1              = parms.syms.y1;
phi1            = parms.syms.phi1;
x2              = parms.syms.x2;
y2              = parms.syms.y2;
phi2            = parms.syms.phi2;
t               = parms.syms.t;

% Generalised state derivative
x1d             = parms.syms.x1d;
y1d             = parms.syms.y1d;
phi1d           = parms.syms.phi1d;
x2d             = parms.syms.x2d;
y2d             = parms.syms.y2d;
phi2d           = parms.syms.phi2d;

% Create generalized coordinate vectors
x               = [x1;y1;phi1;x2;y2;phi2];
xd              = [x1d;y1d;phi1d;x2d;y2d;phi2d];

% Calculate Position constraints
C               = [x1+b*cos(phi1)-x2+d*cos(phi2); ...
    y1+b*sin(phi1)-y2+d*sin(phi2)];

% Calculate Velocity constraints
v1              = [x1d y1d 0;x2d y2d 0].';
omega           = [0 0 phi1d;0 0 phi2d].';
R_A_COM         = [-a*cos(phi1) -a*sin(phi1) 0; c*cos(phi2) c*sin(phi2) 0].';
Va              = v1 + cross(omega,R_A_COM);
eA              = [-sin(phi1) cos(phi1) 0; -sin(phi2) cos(phi2) 0].';
D_x             = simplify([Va(:,1).'*eA(:,1);Va(:,2).'*eA(:,2)]);

% Split constraint in matrix vector product
D               = equationsToMatrix(D_x,[x1d y1d phi1d x2d y2d phi2d]);

% Compute the jacobian of the (non-)holonomic constraints
JC_x            = simplify(jacobian(C,x.'));
JD_x            = simplify(jacobian(D_x,xd.'));

% Calculate convective component
JC_xd           = jacobian(JC_x*xd,x);
JD_xd           = jacobian(D*xd,x);

% Create system of DAE
A = [parms.M JC_x.' D.'                                                     ; ...
    JC_x zeros(size(JC_x,1),size(JC_x.',2)) zeros(size(D,1),size(D.',2)); ...
    D zeros(size(D,1),size(JC_x.',2)) zeros(size(D,1),size(D.',2))];
B = [parms.F ;-JC_xd*xd;-JD_xd*xd];

% Calculate result expressed in generalized coordinates
xdd             = A\B;

%% Convert to function handles
matlabFunction(simplify(xdd),'vars',[x1 y1 phi1 x2 y2 phi2],'vars',[phi1 phi2 x1d y1d phi1d x2d y2d phi2d t],'File','subs_xdd');

% Position constraint function handle
matlabFunction(simplify(C),'vars',[x1 y1 phi1 x2 y2 phi2],'File','subs_C');

% Position constraint derivative function handle
matlabFunction(simplify(JC_x),'File','subs_Cd');

% Velocity constraint  function handle
matlabFunction(simplify(D_x),'vars',[phi1 phi2 x1d y1d phi1d x2d y2d phi2d],'File','subs_D');

% Velocity constraint derivative function handle
matlabFunction(simplify(JD_x),'File','subs_Dd');

% Force torque volocity handle
matlabFunction(parms.F,'File','subs_F');

end