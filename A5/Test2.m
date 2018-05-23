%% ME41055 - Multibody Dynamics B - HW Set 5
% Submitted by Prajish Sekoor Lakshmana Sankar #4743873
% Due 27th March, 2018
 
% Virtual Power Method
% clc; clear all;
 
%% Initialize variables
 
syms l m I g M
 
I = (1/12)*m*(l^2); % mass moment of inertia of each leg about its CM
 
%% Compute the derivatives of constraint matrix
 
% Use symbolic toolbox
syms x1 y1 phi1 x2 y2 phi2 x3 y3 phi3
syms dx1 dy1 dphi1 dx2 dy2 dphi2 dx3 dy3 dphi3
syms alpha beta gamma dalpha dbeta 
 
 
% Defining the genaralized coordinates
q = [alpha beta].';
dq = [dalpha dbeta].';
 
% Expressing coordinates of CM using generalized coordinates
% Here, x = [x1 y1 phi1 x2 y2 phi2 x3 y3 phi3]';
x = [0.5*l*sin(alpha);...
     0.5*l*cos(alpha);...
     pi/2 - alpha;...
     l*sin(alpha);...
     l*cos(alpha);...
     0;...
     l*sin(alpha) + 0.5*l*sin(pi/2 + beta);...
     l*cos(alpha) + 0.5*l*cos(pi/2 + beta);...
     pi/2 + beta];
 
T = simplify(jacobian(x,q)); % Computing the jacobian
 
dx = T*dq; % This is dx/dt=dx/dq*(dq)
 
%% Solving for accelerations 
Q     = 0;                                          % Extra non-potential force terms
M_bar = T.'*diag([m,m,I,M,M,0,m,m,I])*T; 
 
g_conv = jacobian(T*dq,q)*dq;
 
F = [m*g*sin(gamma) -m*g*cos(gamma) 0 M*g*sin(gamma) -M*g*cos(gamma) 0 ...
    m*g*sin(gamma) -m*g*cos(gamma) 0];
 
Z = simplify(T.'*(F.' - diag([m,m,I,M,M,0,m,m,I])*g_conv));
 
ddq = inv(M_bar)*Z;
 
%% Converting back to CM coordinates
 
ddx = simplify(jacobian(dx,dq.'))*ddq + simplify(jacobian(dx,q.'))*dq;
 
%% Substitute values
 
ddx = (subs(ddx, [l,m,g,M,gamma],...
    [0.8,12,9.81,36,pi/12]));
 
ddx = double(subs(ddx, [alpha,dalpha,beta,dbeta],...
    [deg2rad(30),-pi,deg2rad(45),eps]))  % eps for 0, to avoid singularity 
%% END %%
