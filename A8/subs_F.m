function out1 = subs_F(t)
%SUBS_F
%    OUT1 = SUBS_F(T)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    31-May-2018 13:13:01

t2 = t.*pi;
t3 = cos(t2);
out1 = [0.0;0.0;t3.*(-1.0./1.0e1);0.0;0.0;t3.*(1.0./1.0e1)];
