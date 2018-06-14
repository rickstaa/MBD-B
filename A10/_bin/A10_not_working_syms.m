%% MBD_B: Assignment 9 - Euler angles and the human arm
%  Rick Staa (4511328)
clear all; close all; clc;
fprintf('--- A10 ---\n');

%% Simulation settings
EOM_calc_bool               = 1;                                        % Set on 1 if you want to recalculate the EOM

%% Intergration parameters
parms.sim.sim_time          = 5;                                        % Intergration time
parms.sim.dt                = 1e-3;                                     % Intergration step size

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
[parms.k,parms.l0,parms.b] = spring_param_calc(parms);

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

%% Play sound 
load chirp
sound(y,Fs)

%% Calculate movement by mean sof a Runge-Kuta 4th order intergration method
[t,x]                       = RK4_custom(x0,parms);

%% FUNCTIONS

%% Runge-Kuta numerical intergration function
% This function calculates the motion of the system by means of a
% Runge-Kuta numerical intergration. This function takes as inputs the
% parameters of the system (parms), the EOM of the system (parms.EOM)
% and the initial state.
function [time,x] = RK4_custom(x0,parms)

% Initialise variables
time                = (0:parms.sim.dt:parms.sim.sim_time).';               % Create time array
x                   = zeros(length(time),length(x0));                      % Create empty state array
x(1,1:length(x0))   = x0;                                                  % Put initial state in array

% Caculate the motion for the full simulation time by means of a
% Runge-Kutta4 method

% Perform intergration till end of set time
for ii = 1:(size(time,1)-1)
    
    % Add time constant
    t = time(ii);
    
    % Perform RK 4
    x_now_tmp           = x(ii,:);                                                                  % Create cell for subs function function
    x_input             = num2cell(x(ii,:),1);                                                      % Add time to state
    K1                  = [x_now_tmp(1,8:10),subs_xdd(x_input{:}).'];                                % Calculate the second derivative at the start of the step
    x1_tmp              = x_now_tmp + (parms.h*0.5)*K1;                                             % Create cell for subs function function
    x1_input            = num2cell(x1_tmp,1);                                                       % Add time to state
    K2                  = [x1_tmp(1,8:10),subs_xdd(x1_input{:}).'];                                  % Calculate the second derivative halfway the step
    x2_tmp              = x_now_tmp + (parms.h*0.5)*K2;                                             % Refine value calculation with new found derivative
    x2_input            = num2cell(x2_tmp,1);                                                       % Add time to state
    K3                  = [x2_tmp(1,8:10),subs_xdd(x2_input{:}).'];                                  % Calculate new derivative at the new refined location
    x3_tmp              = x_now_tmp + (parms.h)*K3;                                                 % Calculate state at end step with refined derivative
    x3_input            = num2cell(x3_tmp,1);                                                       % Add time to state
    K4                  = [x3_tmp(1,8:10),subs_xdd(x3_input{:}).'];                                  % Calculate last second derivative
    x(ii+1,:)           = x_now_tmp + (parms.h/6)*(K1+2*K2+2*K3+K4);                                % Perform euler intergration step
    
    % Correct for intergration drift (Renormalise the axis of rotation)
    q_next              = x(ii+1,4:7);
    x(ii+1,4:7)         = q_next/norm(q_next);
    
end
end

%% Calculate (symbolic) Equations of Motion four our setup
function EOM_calc(parms)

% Unpack symbolic variables from varargin
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

% Calculate the rotation matrix
% See page 160 of the reader
R_B2N           = [q0^2 + q1^2 - q2^2 - q3^3,       2*(q1*q2-q0*q3),        2*(q1*q3-q0*q2)     ; ...
                   2*(q2*q1-q0*q3),                 q0^2-q1^2+q2^2-q3^2,    2*(q2*q3-q0*q1)     ; ...
                   2*(q3*q1-q0*q2),                 2*(q3*q2-q0*q1),        q0^2-q1^2-q2^2+q3^2]; ...

% Create mass matrix
parms.M         = diag([parms.m,parms.m,parms.m]);

% Express vehicle points in the global (inertial frame N)
r_c             = [x;y;z] + R_B2N*-parms.c;                               % Vehicle center in N frame. -parms.c is the relative position of center w.r.t. COM in N frame
r_s1            = r_c + R_B2N*-parms.p_s;                                 % Left spring attachment point in N frame
r_s2            = r_c + R_B2N*parms.p_s;                                  % Right spring attachment point in N frame

%% Get the virtual power of the spring element
% Create small generalised spring state
x_state         = [x;y;z];
xd_state        = [x_d;y_d;z_d];
w_state         = [omega_x;omega_y;omega_z];

% Calculate vectors pointing along the springs downwards
l_s1            = r_s1 - [0;-0.5*parms.w;parms.h];
l_s2            = r_s2 - [0;0.5*parms.w;parms.h];

% Calculat spring lenghts
s1_length       = sqrt(sum(l_s1.^2));
s2_length       = sqrt(sum(l_s2.^2));

% Calculate delta of lengths
s1_del          = s1_length-parms.l0;
s2_del          = s2_length-parms.l0;

% Calculate spring force sigmas
sigma_s1         = parms.k*s1_del;
sigma_s2         = parms.k*s2_del;

% Calculate jacobian of spring vector to the state
Js1_q            = jacobian(s1_del,x_state);
Js2_q            = jacobian(s2_del,x_state);

% Calculate spring Forces in N frame
Fs1             = (sigma_s1*Js1_q).';
Fs2             = (sigma_s2*Js2_q).';
% Fs            = [Fs1;Fs2];

% % Get spring forces in N frame
% matlabFunction(Fs,'vars',{q0,q1,q2,q3,x,y,z},'File','subs_Fs'); % Create function handle for spring force vector

%% Get the virtual power of the damper element

% Calculate the derivative of the spring change
s1_length_t     = Js1_q*xd_state;
s2_length_t     = Js2_q*xd_state;

% Calculate the damping sigmas
sigma_b1        = parms.b*s1_length_t;
sigma_b2        = parms.b*s2_length_t;

% Calculate damping forces in N frame
Fb1             = (sigma_b1*Js1_q).';
Fb2             = (sigma_b2*Js2_q).';
% Fb              = [Fb1;Fb2];

% % Get damping forces in N frame
% matlabFunction(Fb,'vars',{x,y,z,q0,q1,q2,q3,x_d,y_d,z_d},'File','subs_Fb'); % Create function handle for damping force vector

%% Get the virtual power of the air drag
v               = [x_d;y_d;z_d];
D               = 0.5*parms.rho*parms.A*parms.cd*sqrt(sum(v.^2))*v;

% % Get damping forces in N frame
% matlabFunction(D,'vars',{x,y,z,q0,q1,q2,q3,x_d,y_d,z_d},'File','subs_D'); % Create function handle for damping force vector

%% Other External forces
F_ext           = [0; 0; -parms.m*parms.g];

%% Now the linear acceleartion part of the virtual power expression (d'alambert forces)
xdd         = parms.M\(F_ext-Fs1-Fs2-Fb1-Fb2-D);

%% Now the angular acceleration part of the virtual power expression (newton-Euler equation)
% Create omega state
omega_B         = [omega_x;omega_y;omega_z];

% Calculate moments (Drag and gravity are in COM so no moment)
Fs1_B           = R_B2N.'*(Fs1-Fb1);                                      % Forces expressed in body fixed fram COM
Fs2_B           = R_B2N.'*(Fs2-Fb2);                                      % Forces expressed in body fixed fram COM
r_Fs1_B         = -parms.c-parms.p_s;                                     % Arm from COM to spring 1 attachment
r_Fs2_B         = -parms.c+parms.p_s;                                     % Arm form COM to spring 2 attachment
M_1             = cross(r_Fs1_B,Fs1_B);                                   % Moment spring 1 exerted on body in frame B
M_2             = cross(r_Fs2_B,Fs2_B);                                   % Moment spring 2 exerted on body in frame B

% Calculate angular acceleration of the body
omega_dd        = parms.J\(M_1-M_2 -cross(omega_B,parms.J*omega_B));

% Calculate the derivatives of the euler parameters
Q               = [q0 -q1 -q2  -q3; ...
                   q1  q0 -q3   q2; ...
                   q2  q3  q0  -q1; ...
                   q3 -q2  q1   q0];
qd              = 0.5*Q*[0;w_state]; 
               
% Put state derivative in one vector
Xdd             = [xd_state;qd;xdd;omega_dd];

% Get damping forces in N frame
tic;
matlabFunction(Xdd,'vars',{x,y,z,q0,q1,q2,q3,x_d,y_d,z_d,omega_x,omega_y,omega_z},'File','subs_xdd');                                    % Create function handle for damping force vector
toc;

end

%% Calculate spring-damper paramters

function [k,l0,b] = spring_param_calc(parms)

%% A) Calculate k and L0 paramters
% This can be done by using the principle of conservation of energy and
% newtons second law of motion.

% Create symbolic variables
syms k l0

% Calculate needed spring lengths
l_t       = sqrt((0.5*parms.w)^2 + (0.5*parms.w-parms.r)^2);
l_b       = sqrt(parms.h^2 + (0.5*parms.w-parms.r)^2);

% Conservation of energy equation
eq_con      = k*(l_t - l0)^2+parms.m*parms.g*(0.5*parms.w+parms.h)-k*(l_b-l0)^2;

% Force equation (Newtons second law of motion)
eq_F        = 2*k*(l_b-l0)*cos(atan2((0.5*parms.w-parms.r),parms.h))-5.8*parms.m*parms.g;

% Put equations in one vector
eq          = [eq_con; eq_F];

% Solve these two equations for two unknowns with the symbolic toolbox
sol         = solve(eq,[k;l0]);
k           = double(sol.k);
l0          = double(sol.l0);

% Calculate damping ratio
b           = 2*parms.zeta*sqrt(k*parms.m);
end