function out1 = subs_Cd(phi1,phi2)
%SUBS_CD
%    OUT1 = SUBS_CD(PHI1,PHI2)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    01-Jun-2018 15:35:06

out1 = reshape([1.0,0.0,0.0,1.0,sin(phi1).*(-1.0./2.0),cos(phi1).*(1.0./2.0),-1.0,0.0,0.0,-1.0,sin(phi2).*(-1.0./8.0),cos(phi2).*(1.0./8.0)],[2,6]);
