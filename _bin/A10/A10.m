%% MBD_B: Assignment 9 - Euler angles and the human arm
%  Rick Staa (4511328)
% clear all; close all; clc;
fprintf('--- A10 ---\n');

%% Simulation settings
EOM_calc_bool               = 0;                                         % Set on 1 if you want to recalculate the EOM

%% Intergration parameters
parms.sim.sim_time          = 60;                                        % Intergration time
parms.sim.dt                = 1e-3;                                      % Intergration step size
parms.sim.nmax              = 10;                                        % Max number of coordinate projections
parms.sim.tol               = 1e-12;                                     % Maximum allowable drift

%% Model Parameters
% Vehicle parameters
m               = 420;                          % mass vehicle [kg]
J               = diag([170 120 140]);          % moment of inertia of vehicle [kg]
r               = 1;                            % Radius of vehicle [m]
c               = [-0.01;0.01;-0.1];            % Relative COM position to vehicle center [m]
p_s             = [0;r;0];                      % Connection of bungie cords relative 2 sphere center
A               = pi*r^2;                       % Frontal area of the vehicle [m]

% Paramters of support collumns
h               = 25;                           % Hight of supporting collumns [m]
w               = 18;                           % Width between the supporting collums [m]
zeta            = 0.1;                          % Damping ratio

%% World parameters
g               = 9.81;                         % N/kg
rho             = 1.25;                         % kg/m^3
cd              = 0.5;                          % Drag coefficient

%% put parameters in struct
parms.m         = m;
parms.J         = J;
parms.r         = r;
parms.c         = c;
parms.A         = A;
parms.h         = h;
parms.w         = w;
parms.zeta      = zeta;
parms.g         = g;
parms.rho       = rho;
parms.cd        = cd;
parms.p_s       = p_s;

% calculate the missing spring parameters
[parms.k,parms.l_0,parms.b] = spring_param_calc(parms);

% Create xtra symbolic variables
syms x y z q0 q1 q2 q3 x_d y_d z_d omega_x omega_y omega_z            % In this q0 = lambda 0 this was done for code Readability

% Put symbolic variables in struct
parms.syms.x        = x;
parms.syms.y        = y;
parms.syms.z        = z;
parms.syms.q0       = q0;
parms.syms.q1       = q1;
parms.syms.q2       = q2;
parms.syms.q3       = q3;
parms.syms.x_d      = x_d;
parms.syms.y_d      = y_d;
parms.syms.z_d      = z_d;
parms.syms.omega_x  = omega_x;
parms.syms.omega_y  = omega_y;
parms.syms.omega_z  = omega_z;

%% Set Initial states
% Set euler parameters (In the initial state the axis of ration can said to
% be alighed with the axis through the spring attachment sites.
n                   = [0;1;0];                                           % Axis through spring attachment sites
phi                 = 0;                                                 % No rotation

% Calculate initial states
x0                  = 0;
y0                  = 0;
z0                  = 0;
q0                  = cos(0.5*phi);
q1                  = sin(0.5*phi)*n(1);
q2                  = sin(0.5*phi)*n(2);
q3                  = sin(0.5*phi)*n(3);
x_d                 = 0;
y_d                 = 0;
z_d                 = 0;
omega_x             = 0;
omega_y             = 0;
omega_z             = 0;
x0                  = [x0;y0;z0;q0;q1;q2;q3;x_d;y_d;z_d;omega_x;omega_y;omega_z];

%% Calculate equations of motion
if (EOM_calc_bool == 1)
    EOM_calc(parms);
end

%% Calculate movement by mean sof a Runge-Kuta 4th order intergration method
[t,x,error,r_axis]         = RK4_custom(x0,parms);

%% Play sound
load gong
sound(y,Fs)

%% Optimize k b and c values
% Try to find the optimal values that get the highest height

% % Create arrays
% flag = 1;
% h_max_val = 0;
% l_0_array = 13:-0.5:8;
% k_array   = 500:50:1000;
% parms.sim.sim_time = 20; % Set to 20 seconds that should be more than enough
% for ii = 1:length(l_0_array)
%     % set l_0
%     parms.l_0 = l_0_array(ii);
%     
%     for jj = 1:length(k_array)
%         % Set k
%         parms.k = k_array(jj);
%         
%         % Calculate the b that corresponds to a zeta of 0.1 and the given k
%         parms.b           = 2*parms.zeta*sqrt(parms.k*m);
%         
%         % Plot for new values
%         [~,x_opt,~,~] = RK4_custom(x0,parms);
%         
%         % Get height and save in array
%         if max(x_opt(:,3)) > h_max_val
%             h_max_val   = max(x_opt(:,3));
%             h_max.k   = parms.k;
%             h_max.b   = parms.b;
%             h_max.l_0 = parms.l_0;
%             
%             % Exit loop if
%             if h_max_val > (parms.h+0.5*parms.w)
%                 flag = 0;
%             end
%             
%             % break out of loop if flag = 0;
%             if flag == 0
%                 break
%             end
%         end
%         % break out of loop if flag = 0;
%         if flag == 0
%             break
%         end
%     end
%     % break out of loop if flag = 0;
%     if flag == 0
%         break
%     end
% end


%% Calculate energies
%% Energies

% Potential Energy
E_p     = m*g*x(:,3);   % Global axis is at the ground and z is relative to N frame

% Spring Energy
l_s1    = subs_l_s1(parms.l_0,x(:,1),x(:,2),x(:,3),x(:,4),x(:,5), ...
    x(:,6),x(:,7));
l_s2    = subs_l_s2(parms.l_0,x(:,1),x(:,2),x(:,3),x(:,4),x(:,5), ...
    x(:,6),x(:,7));

E_e     = 0.5*parms.k*(l_s1 - parms.l_0).^2 + 0.5*parms.k*(l_s2 - parms.l_0).^2;

% Translational Kinetic Energy
E_t     = 0.5*m*(x(:,8).^2 + x(:,9).^2 + x(:,10));

% Rotational Kinetic Energy
E_r     = 0.5*parms.J(1,1)*x(:,11).^2 + 0.5*parms.J(2,2)*x(:,12).^2 + 0.5*parms.J(3,3)*x(:,13).^2;

%% Create plots

%% Plots for the vehicle movement
% Plot x y z positions and velocities of the vehicle
figure;
subplot(3,1,1)
plot(t,x(:,1),t,x(:,8),'linewidth',2);
set(gca,'fontsize',14);
title('Linear properties of the COM of vehicle')
xlabel('time [s]');
ylabel('x');
legend('Position [m]','Velocity [m/s]');
subplot(3,1,2)
plot(t,x(:,2),t,x(:,9),'linewidth',2);
set(gca,'fontsize',14);
xlabel('time [s]');
ylabel('y');
legend('Position [m]','Velocity [m/s]');
subplot(313)
plot(t,x(:,3),t,x(:,10),'linewidth',2);
set(gca,'fontsize',14);
xlabel('time [s]');
ylabel('z');
legend('Position [m]','Velocity [m/s]');

% Plot angular velocities of the vehicle
figure(2)
subplot(3,1,1)
plot(t,x(:,11),'linewidth',2);
set(gca,'fontsize',14);
title('Angular Velocities expressed in the Body Frame')
xlabel('time [s]');
ylabel('\omega_x [rad/s]');
subplot(3,1,2)
plot(t,x(:,12),'linewidth',2);
set(gca,'fontsize',14);
xlabel('time');
ylabel('\omega_y [rad/s]');
subplot(3,1,3)
plot(t,x(:,13),'linewidth',2);
set(gca,'fontsize',14);
xlabel('time [s]');
ylabel('\omega_z [rad/s]');

% Plot trajectory of vehicle in 3D
figure;
plot3(x(:,1),x(:,2),x(:,3));
hold on;
plot3(x(1,1),x(1,2),x(1,3),'g*','linewidth',8); % Plot beginning
plot3(x(end,1),x(end,2),x(end,3),'r*','linewidth',8);
set(gca,'fontsize',14);
title('Trajectory of vehicle in 3D')
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
legend('vehicle trajectory','start','end');
axis equal;

% Plot spring force magnitude
figure;
subplot(4,1,1)
plot(t,subs_F_spring(parms.k,parms.l_0,x(:,1),x(:,2),x(:,3),x(:,4), ...
    x(:,5),x(:,6),x(:,7)),'linewidth',2);
set(gca,'fontsize',14);
title('Forces (in N frame) and Moments (in B frame)')
ylabel('F_{spring} [N]');
subplot(4,1,2)
plot(t,subs_F_damp(parms.b,x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6), ...
    x(:,7),x(:,8),x(:,9),x(:,10)),'linewidth',2);
set(gca,'fontsize',14);
ylabel('F_{damping} [N]');
subplot(4,1,3)
plot(t,sum(subs_F_drag(x(:,8)',x(:,9)',x(:,10)').^2),'linewidth',2);
set(gca,'fontsize',14);
ylabel('F_{drag} [N]');
subplot(4,1,4)
plot(t,sum(subs_M(parms.k,parms.l_0,parms.b,x(:,1)',x(:,2)',x(:,3)', ...
    x(:,4)',x(:,5)',x(:,6)',x(:,7)',x(:,8)',x(:,9)',x(:,10)').^2) ...
    ,'linewidth',2);
set(gca,'fontsize',14);
ylabel('M_{total} [Nm]');
xlabel('time [s]');

% Plot the accuracy of the solution
figure;
plot(t,abs(error));
hold on;
plot([0 parms.sim.sim_time],[1e-12 1e-12],'linewidth',2,'linestyle','--');
set(gca,'fontsize',14);
ylim([0 1.15e-12]);
xlabel('time [s]');
ylabel('Intergration drift (Error on normality constraint)');
title('Intergration accuracy of the solution against the time')
legend('Intergration error','max allowed error line','location','best');

%% Plot for the energies

% First plot the individual energies
figure;
subplot(4,1,1)
plot(t,E_p,'linewidth',2);
set(gca,'fontsize',14);
title('Individual energies')
ylabel('E_p [J]');
subplot(4,1,2)
plot(t,E_e,'linewidth',2);
set(gca,'fontsize',14);
ylabel('E_e [J]');
xlabel('time [s]');
subplot(4,1,3)
plot(t,E_t,'linewidth',2);
set(gca,'fontsize',14);
ylabel('E_t [J]');
xlabel('time [s]');
subplot(4,1,4)
plot(t,E_r,'linewidth',2);
set(gca,'fontsize',14);
ylabel('E_r [J]');
xlabel('time [s]');

% Now plot interesting combinations of energy
figure;
plot(t,E_p,'linewidth',2,'linestyle','--');hold on % Gravitational
plot(t,E_e,'linewidth',2,'linestyle','--');hold on % Gravitational
plot(t,E_t,'linewidth',2,'linestyle','--');hold on;
plot(t,E_r,'linewidth',2,'linestyle','--');hold on;
plot(t,E_p+E_e,'linewidth',2,'linestyle','--'); hold on; % Full potential energy
plot(t,E_t+E_r,'linewidth',2,'linestyle','--'); hold on;
plot(t,E_p+E_e+E_t+E_r,'linewidth',2,'linestyle','-');
set(gca,'fontsize',14);
title('All energies in one plot');
legend('E_p','E_e','E_t','E_r','E_p+E_e','E_t+E_r','Combined energy');
ylabel('Energie [J]');
xlabel('time [s]');

%% Plot xomponents of the axis through the spring attachment sites expressed in the global fram
figure;
subplot(3,1,1);
set(gca,'fontsize',14)
title('Check for looping')
plot(t,r_axis(:,1),'r');
ylabel('r_{cx} [m]');
subplot(3,1,2);
plot(t,r_axis(:,2),'b');
ylabel('r_{cy} [m]');
subplot(3,1,3);
plot(t,r_axis(:,3),'g');
ylabel('r_{cz} [m]');

%% FUNCTIONS

%% Runge-Kuta numerical intergration function
% This function calculates the motion of the system by means of a
% Runge-Kuta numerical intergration. This function takes as inputs the
% parameters of the system (parms), the EOM of the system (parms.EOM)
% and the initial state.
function [time,x,error,r_axis] = RK4_custom(x0,parms)

% Initialise variables
time                = (0:parms.sim.dt:parms.sim.sim_time).';               % Create time array
x                   = zeros(length(time),length(x0));                      % Create empty state array
x(1,1:length(x0))   = x0;                                                  % Put initial state in array
error               = zeros(length(time),1);

% preallocate memory for rotation axis
R_tmp = subs_R_B_N(x0(4),x0(5),x0(6),x0(7));
r_axis              = zeros(length(time),3);                                
r_axis(1,:)         = (R_tmp*parms.c)';
   
% Caculate the motion for the full simulation time by means of a
% Runge-Kutta4 method

% Perform intergration till end of set time
for ii = 1:(size(time,1)-1)
    
    % Perform RK 4
    x_now_tmp           = x(ii,:);                                                          % Create cell for subs function function
    x_input             = num2cell([parms.k,parms.l_0,parms.b,x(ii,:)],1);                  % Add time to state
    K1                  = subs_Xdd(x_input{:}).';                                           % Calculate the second derivative at the start of the step
    x1_tmp              = x_now_tmp + (parms.sim.dt*0.5)*K1;                                % Create cell for subs function function
    x1_input            = num2cell([parms.k,parms.l_0,parms.b,x1_tmp],1);                   % Add time to state
    K2                  = subs_Xdd(x1_input{:}).';                                          % Calculate the second derivative halfway the step
    x2_tmp              = x_now_tmp + (parms.sim.dt*0.5)*K2;                                % Refine value calculation with new found derivative
    x2_input            = num2cell([parms.k,parms.l_0,parms.b,x2_tmp],1);                   % Add time to state
    K3                  = subs_Xdd(x2_input{:}).';                                          % Calculate new derivative at the new refined location
    x3_tmp              = x_now_tmp + (parms.sim.dt)*K3;                                    % Calculate state at end step with refined derivative
    x3_input            = num2cell([parms.k,parms.l_0,parms.b,x3_tmp],1);                   % Add time to state
    K4                  = subs_Xdd(x3_input{:}).';                                          % Calculate last second derivative
    x(ii+1,:)           = x_now_tmp + (parms.sim.dt/6)*(K1+2*K2+2*K3+K4);                  % Perform euler intergration step
    
    % Correct for intergration drift (Renormalise the axis of rotation)
    
    %% Coordinate projection method (Gaus-newton method)
    % Correct for intergration drift
    x_now_tmp = x(ii+1,:);
    [x_new,error_tmp] = gauss_newton(x_now_tmp,parms);
    
    % Calculate the roation of the vector pointing from the center to one
    % of the springs
    
    % For Question f
    R_tmp = subs_R_B_N(x_new(4),x_new(5),x_new(6),x_new(7));
    r_axis(ii,:) = (R_tmp*parms.c)';
    
    % Overwrite position coordinates
    x(ii+1,:)       = x_new;
    error(ii+1,:)   = error_tmp;
end
end

%% Speed correct function
function [x,error] = gauss_newton(x,parms)

% Get rid of the drift by solving a non-linear least square problem by
% means of the Gaus-Newton method

%% Gaus-newton velocity constraint correction
n_iter          = 0;                             % Set iteration counter                                                               % Get position data out

% % Calculate the two needed constraints
D                   = subs_D(x(4),x(5),x(6),x(7));
Dd                  = subs_Dd(x(4),x(5),x(6),x(7));

% Solve drift using gaus-newton iteration
while (max(abs(D)) > parms.sim.tol)&& (n_iter < parms.sim.nmax)
    x_tmp           = x;
    n_iter          = n_iter + 1;
    x_del           = Dd.'*inv(Dd*Dd.')*-D;
    x               = x_tmp + x_del.';
    
    % Recalculate constraint
    D                   = subs_D(x(4),x(5),x(6),x(7));
    Dd                  = subs_Dd(x(4),x(5),x(6),x(7));
end


% Store full error
error = D;
end

%% Calculate (symbolic) Equations of Motion four our setup
function EOM_calc(parms)

%% Get parameters and variables

% create symbolic variables
syms k b l_0;

% Unpack variables for clarity
m               = parms.m;
J               = parms.J;
r               = parms.r;
c               = parms.c;
A               = parms.A;
h               = parms.h;
w               = parms.w;
g               = parms.g;
rho             = parms.rho;
c_d             = parms.cd;
p_s             = parms.p_s;

% Unpack symbolic variables from parms
x               = parms.syms.x;
y               = parms.syms.y;
z               = parms.syms.z;
q0              = parms.syms.q0;
q1              = parms.syms.q1;
q2              = parms.syms.q2;
q3              = parms.syms.q3;
x_d             = parms.syms.x_d;
y_d             = parms.syms.y_d;
z_d             = parms.syms.z_d;
omega_x         = parms.syms.omega_x;
omega_y         = parms.syms.omega_y;
omega_z         = parms.syms.omega_z;

% Create small generalised spring state
x_state         = [x;y;z];
xd_state        = [x_d;y_d;z_d];
omega_state     = [omega_x;omega_y;omega_z];

%% Calculate Rotation Matrix

R_B_N = [q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3),      2*(q1*q3-q0*q2);     ...
    2*(q2*q1-q0*q3),     q0^2-q1^2+q2^2-q3^2 , 2*(q2*q3-q0*q1);     ...
    2*(q3*q1-q0*q2),     2*(q3*q2-q0*q1),      q0^2-q1^2-q2^2+q3^2];

matlabFunction(R_B_N,'File','subs_R_B_N');

%% Calculate Spring contribution to virtual work
r_c             = [x;y;z] + R_B_N*-parms.c;
r_s1            = r_c + R_B_N*-parms.p_s;
r_s2            = r_c + R_B_N*parms.p_s;

% Create vetors along the spring
l_s1            = r_s1 - [0; -w/2; h];
l_s2            = r_s2 - [0; w/2; h];

% Calculate delta spring lenght
l_s2            = sqrt(sum(l_s2.^2)) - l_0;
l_s1            = sqrt(sum(l_s1.^2)) - l_0;

% Calculate spring forces
sigma_s1        = k*l_s1;
sigma_s2        = k*l_s2;

%% Damping Components
sigma_c1        = b*jacobian(l_s1,x_state)*xd_state;
sigma_c2        = b*jacobian(l_s2,x_state)*xd_state;

%% Air Drag
sigma_a         = 0.5*rho*A*c_d*sqrt(sum(xd_state.^2))*xd_state;

%% Add the contributions of all forces together

% External intertial forces
F_ext           = [0;0;-m*g];

% Spring, damper and drag forces
F               = [F_ext - sigma_s1*jacobian(l_s1,x_state)' - sigma_s2*jacobian(l_s2,x_state)' - ...
    sigma_c1*jacobian(l_s1,x_state)' - sigma_c2*jacobian(l_s2,x_state)' - sigma_a];

%% Now finaly calculate the linear accelerations by means of gaussian elimination
M               = diag([m m m]);
xdd_lin         = inv(M)*F;

%% Now lets find the angular accelerations
% Calculate forces
F1              = sigma_s1*jacobian(l_s1,x_state).' + sigma_c1*jacobian(l_s1,x_state).';
F2              = sigma_s2*jacobian(l_s2,x_state).' + sigma_c2*jacobian(l_s2,x_state).';

% Calculate moment arms
r1              = -parms.c-parms.p_s;
r2              = -parms.c+parms.p_s;

% Calculate moments
M               = cross(r1,R_B_N'*F1) + cross(r2,R_B_N'*F2);

%% Calculate angular accelerations

omega_d_state   = inv(J)*(M - cross(omega_state,J*omega_state));

%% Calculate the derivative of the omegas and put them in one big state variable
X_state         = [x,y,z,q0,q1,q2,q3,x_d,y_d,z_d,omega_x,omega_y,omega_z].';

dlambda         = 0.5*[q0 -q1 -q2 -q3;...
    q1 q0 -q3 q2;...
    q2 q3 q0 -q1;...
    q3 -q2 q1 q0]*[0;omega_state];
Xdd             = [xd_state;dlambda;xdd_lin;omega_d_state];

%% Calculate the constraint
D               = q0^2 + q1^2 + q2^2 + q3^2 - 1;
Dd              = jacobian(D,X_state);

% Save all the symbolic expressions in function files
matlabFunction(D,'File','subs_D','vars',[q0 q1 q2 q3]);
matlabFunction(Dd,'File','subs_Dd','vars',[q0 q1 q2 q3]);
matlabFunction(l_s1,'File','subs_l_s1','vars',[l_0,x,y,z,q0,q1,q2,q3]);
matlabFunction(l_s2,'File','subs_l_s2','vars',[l_0,x,y,z,q0,q1,q2,q3]);
matlabFunction((sigma_s1+sigma_s2),'File','subs_F_spring','vars',[k,l_0,x,y,z,q0,q1,q2,q3]);
matlabFunction((sigma_c1+sigma_c2),'File','subs_F_damp','vars',[b,x,y,z,q0,q1,q2,q3,x_d,y_d,z_d]);
matlabFunction((sigma_a),'File','subs_F_drag');
matlabFunction(M,'File','subs_M','vars',[k,l_0,b,x,y,z,q0,q1,q2,q3,x_d,y_d,z_d]);
matlabFunction(X_state,'File','subs_X_state','Vars',[x,y,z,q0,q1,q2,q3,x_d,y_d,z_d,omega_x,omega_y,omega_z]);
matlabFunction(Xdd,'File','subs_Xdd','Vars',[k,l_0,b,x,y,z,q0,q1,q2,q3,x_d,y_d,z_d,omega_x,omega_y,omega_z]);

end

%% Calculate spring-damper paramters

function [k,l0,b] = spring_param_calc(parms)
% Unpack variables for clarity
m               = parms.m;
J               = parms.J;
r               = parms.r;
c               = parms.c;
A               = parms.A;
zeta            = parms.zeta;
h               = parms.h;
w               = parms.w;
g               = parms.g;
rho             = parms.rho;
cd              = parms.cd;
p_s             = parms.p_s;

%% A) Calculate k and L0 paramters
% This can be done by using the principle of conservation of energy and
% newtons second law of motion.

% Create symbolic variables
syms k l0

% Calculate needed spring lengths
l_t       = sqrt((0.5*w)^2 + (0.5*w-r)^2);
l_b       = sqrt(h^2 + (0.5*w-r)^2);

% Conservation of energy equation
eq_con      = k*(l_t - l0)^2+m*g*(0.5*w+h)-k*(l_b-l0)^2;

% Force equation (Newtons second law of motion)
eq_F        = 2*k*(l_b-l0)*cos(atan2((0.5*w-r),h))-5.8*m*g;

% Put equations in one vector
eq          = [eq_con; eq_F];

% Solve these two equations for two unknowns with the symbolic toolbox
sol         = solve(eq,[k;l0]);
k           = double(sol.k);
l0          = double(sol.l0);

% Calculate damping ratio
b           = 2*zeta*sqrt(k*m);
end