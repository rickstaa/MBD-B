%% Runge-Kuta numerical intergration function
% This function calculates the motion of the system by means of a
% Runge-Kuta numerical intergration. This function takes as inputs the
% parameters of the system (parms), the EOM of the system (parms.EOM)
% and the initial state.
function [time,x,error] = RK4_custom(x0,parms)

% Initialise variables
time                = (0:parms.sim.dt:parms.sim.sim_time).';               % Create time array
x                   = zeros(length(time),length(x0));                      % Create empty state array
x(1,1:length(x0))   = x0;                                                  % Put initial state in array
error               = zeros(length(time),1);

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
    
    % Overwrite position coordinates
    x(ii+1,:)       = x_new;
    error(ii+1,:)   = error_tmp;
end
end