%% MBD_B: Assignment 9 - Euler angles and the human arm
%  Rick Staa (4511328)
% clear all; close all; clc;
fprintf('--- A10 ---\n');

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
EOM_calc(parms);
    
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