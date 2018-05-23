function [t,q] = gaus_newton(t,q,varargin)

% Get rid of the drift by solving a non-linear least square problem by
% means of the Gaus-Newton method
C       = feval(parms.C_handle,q(1),q(2),q(3));

end