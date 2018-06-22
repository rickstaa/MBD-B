function l_s2 = subs_l_s2(l_0,x,y,z,q0,q1,q2,q3)
%SUBS_L_S2
%    L_S2 = SUBS_L_S2(L_0,X,Y,Z,Q0,Q1,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    22-Jun-2018 10:48:02

t3 = q0.^2;
t4 = q1.^2;
t5 = q2.^2;
t6 = q3.^2;
t2 = t3.*(1.0./1.0e2)+t4.*(1.0./1.0e2)-t5.*(1.0./1.0e2)-t6.*(1.0./1.0e2)+x-q0.*q2.*(1.0./5.0)-q0.*q3.*(9.9e1./5.0e1)+q1.*q2.*(9.9e1./5.0e1)+q1.*q3.*(1.0./5.0);
t7 = t3.*(9.9e1./1.0e2)-t4.*(9.9e1./1.0e2)+t5.*(9.9e1./1.0e2)-t6.*(9.9e1./1.0e2)+y-q0.*q1.*(1.0./5.0)-q0.*q3.*(1.0./5.0e1)+q1.*q2.*(1.0./5.0e1)+q2.*q3.*(1.0./5.0)-9.0;
t8 = t3.*(1.0./1.0e1)-t4.*(1.0./1.0e1)-t5.*(1.0./1.0e1)+t6.*(1.0./1.0e1)+z-q0.*q1.*(9.9e1./5.0e1)-q0.*q2.*(1.0./5.0e1)+q1.*q3.*(1.0./5.0e1)+q2.*q3.*(9.9e1./5.0e1)-2.5e1;
l_s2 = -l_0+sqrt(t2.^2+t7.^2+t8.^2);
