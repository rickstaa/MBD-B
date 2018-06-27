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