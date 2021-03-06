function out1 = subs_F_spring(k,l_0,x,y,z,q0,q1,q2,q3)
%SUBS_F_SPRING
%    OUT1 = SUBS_F_SPRING(K,L_0,X,Y,Z,Q0,Q1,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    26-Jun-2018 10:50:30

t3 = q0.^2;
t4 = q1.^2;
t5 = q2.^2;
t6 = q3.^2;
t9 = q0.*q2.*(1.0./5.0);
t10 = q1.*q3.*(1.0./5.0);
t11 = t3.*(1.0./1.0e2);
t12 = t4.*(1.0./1.0e2);
t13 = t5.*(1.0./1.0e2);
t14 = t6.*(1.0./1.0e2);
t2 = -t9+t10+t11+t12-t13-t14+x-q0.*q3.*(9.9e1./5.0e1)+q1.*q2.*(9.9e1./5.0e1);
t16 = q0.*q1.*(1.0./5.0);
t17 = q0.*q3.*(1.0./5.0e1);
t18 = q1.*q2.*(1.0./5.0e1);
t19 = q2.*q3.*(1.0./5.0);
t7 = t3.*(9.9e1./1.0e2)-t4.*(9.9e1./1.0e2)+t5.*(9.9e1./1.0e2)-t6.*(9.9e1./1.0e2)-t16-t17+t18+t19+y-9.0;
t21 = q0.*q2.*(1.0./5.0e1);
t22 = q1.*q3.*(1.0./5.0e1);
t23 = t3.*(1.0./1.0e1);
t24 = t4.*(1.0./1.0e1);
t25 = t5.*(1.0./1.0e1);
t26 = t6.*(1.0./1.0e1);
t8 = -t21+t22+t23-t24-t25+t26+z-q0.*q1.*(9.9e1./5.0e1)+q2.*q3.*(9.9e1./5.0e1)-2.5e1;
t15 = -t9+t10+t11+t12-t13-t14+x+q0.*q3.*(1.01e2./5.0e1)-q1.*q2.*(1.01e2./5.0e1);
t20 = t3.*(-1.01e2./1.0e2)+t4.*(1.01e2./1.0e2)-t5.*(1.01e2./1.0e2)+t6.*(1.01e2./1.0e2)-t16-t17+t18+t19+y+9.0;
t27 = t21-t22-t23+t24+t25-t26-z-q0.*q1.*(1.01e2./5.0e1)+q2.*q3.*(1.01e2./5.0e1)+2.5e1;
out1 = -k.*(l_0-sqrt(t2.^2+t7.^2+t8.^2))-k.*(l_0-sqrt(t15.^2+t20.^2+t27.^2));
