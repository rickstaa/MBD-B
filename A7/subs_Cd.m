function out1 = subs_Cd(phi2,phi4,phi5)
%SUBS_CD
%    OUT1 = SUBS_CD(PHI2,PHI4,PHI5)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    27-Jun-2018 12:36:49

t2 = sin(phi2);
t3 = t2.*1.2e1;
t4 = t3+1.3e1;
t5 = cos(phi4);
out1 = reshape([t2.*(-1.0./5.0)-1.0./sqrt(t4).*t5.*cos(phi2).*(3.0./5.0),0.0,sqrt(t4).*sin(phi4).*(1.0./1.0e1),t5.*(7.0./1.0e1),0.0,cos(phi5).*(3.0./5.0)],[2,3]);
