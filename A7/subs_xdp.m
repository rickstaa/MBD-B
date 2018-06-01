function out1 = subs_xdp(phi2,phi4,phi5,phi2d,phi4d,phi5d)
%SUBS_XDP
%    OUT1 = SUBS_XDP(PHI2,PHI4,PHI5,PHI2D,PHI4D,PHI5D)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    01-Jun-2018 15:51:25

t2 = cos(phi4);
t3 = cos(phi5);
t4 = cos(phi2);
t5 = t2.^2;
t6 = t4.^2;
t7 = t3.^2;
t8 = sin(phi4);
t9 = sin(phi5);
t10 = sqrt(1.3e1);
t11 = sin(phi2);
t12 = phi2d.^2;
t13 = phi4d.^2;
t14 = phi5d.^2;
t15 = t5.*5.782e3;
t16 = t7.*2.03379e5;
t17 = t5.*t6.*t7.*2.352e3;
t18 = t2.*t3.*t6.*t8.*t9.*2.352e3;
t21 = t5.*t6.*5.782e3;
t22 = t6.*t7.*8.34e3;
t23 = t5.*t7.*1.97391e5;
t24 = t2.*t3.*t8.*t9.*2.352e3;
t19 = t15+t16+t17+t18-t21-t22-t23-t24;
t20 = 1.0./t19;
t25 = phi4.*2.0;
t26 = sin(t25);
t27 = phi5.*2.0;
t28 = sin(t27);
t29 = t13.*t28.*5.88e3;
t30 = t7.*t8.*4.2e6;
t31 = t3.*t8.*t14.*1.008e4;
t32 = t2.*t6.*t8.*t13.*5.782e4;
t33 = t2.*t6.*t9.*t14.*4.956e4;
t34 = t2.*t7.*t8.*t13.*1.97391e6;
t35 = t2.*t3.*t6.*t9.*4.2e6;
t36 = t3.*t5.*t6.*t9.*t13.*2.352e4;
t37 = t4.*t7.*t8.*t10.*t12.*3.0006e5;
t38 = t29+t30+t31+t32+t33+t34+t35+t36+t37+t2.*t7.*2.23668e5-t13.*t26.*2.891e4-t2.*t3.*t9.*4.2e6-t2.*t6.*t7.*2.23668e5-t6.*t7.*t8.*4.2e6-t2.*t9.*t14.*4.956e4-t3.*t5.*t9.*t13.*2.352e4-t3.*t6.*t8.*t14.*1.008e4-t3.*t6.*t9.*t13.*1.176e4-t2.*t6.*t7.*t8.*t13.*2.352e4-t4.*t7.*t8.*t10.*t11.*2.943e3;
t39 = t5.*2.94e7;
t40 = t3.*t5.*t6.*t9.*1.565676e6;
t41 = t2.*t3.*t4.*t8.*t9.*t10.*t11.*2.0601e4;
out1 = [t20.*(t4.*t7.*3.8259e4-t4.*t5.*t7.*3.8259e4-t7.*t10.*t11.*4.2e6+t4.*t5.*t11.*t12.*1.1564e5+t5.*t7.*t10.*t11.*4.2e6+t4.*t7.*t11.*t12.*1.668e5+t2.*t10.*t11.*t13.*5.782e4-t3.*t10.*t11.*t14.*1.008e4-t2.*t7.*t8.*t10.*t11.*2.23668e5-t4.*t5.*t7.*t11.*t12.*4.704e4+t2.*t7.*t10.*t11.*t13.*5.988e4+t3.*t5.*t10.*t11.*t14.*1.008e4+t2.*t3.*t8.*t9.*t10.*t11.*4.2e6+t2.*t8.*t9.*t10.*t11.*t14.*4.956e4-t3.*t8.*t9.*t10.*t11.*t13.*1.176e4-t2.*t3.*t4.*t8.*t9.*t11.*t12.*4.704e4).*(-1.0./2.0e1);t20.*(t7.*t10.*4.2e6-t4.*t7.*t11.*3.8259e4-t5.*t7.*t10.*4.2e6+t4.*t7.*t12.*3.90078e6-t6.*t7.*t10.*4.2e6-t2.*t10.*t13.*5.782e4+t3.*t10.*t14.*1.008e4+t2.*t7.*t8.*t10.*2.23668e5+t4.*t5.*t7.*t11.*3.8259e4-t4.*t5.*t7.*t12.*3.90078e6+t5.*t6.*t7.*t10.*4.2e6+t2.*t6.*t10.*t13.*5.782e4-t2.*t7.*t10.*t13.*5.988e4-t3.*t5.*t10.*t14.*1.008e4-t3.*t6.*t10.*t14.*1.008e4-t2.*t3.*t8.*t9.*t10.*4.2e6-t2.*t6.*t7.*t8.*t10.*2.23668e5+t2.*t6.*t7.*t10.*t13.*5.988e4+t3.*t5.*t6.*t10.*t14.*1.008e4-t2.*t8.*t9.*t10.*t14.*4.956e4+t3.*t8.*t9.*t10.*t13.*1.176e4+t2.*t3.*t6.*t8.*t9.*t10.*4.2e6+t2.*t6.*t8.*t9.*t10.*t14.*4.956e4-t3.*t6.*t8.*t9.*t10.*t13.*1.176e4).*(-1.0./1.0e2);t20.*(t6.*t7.*3.8259e4-t5.*t6.*t7.*3.8259e4+t5.*t11.*t12.*1.1564e5+t7.*t11.*t12.*4.06758e6-t4.*t7.*t10.*t11.*4.2e6-t5.*t7.*t11.*t12.*3.94782e6+t4.*t5.*t7.*t10.*t11.*4.2e6+t2.*t4.*t10.*t11.*t13.*5.782e4-t3.*t4.*t10.*t11.*t14.*1.008e4-t2.*t4.*t7.*t8.*t10.*t11.*2.23668e5-t2.*t3.*t8.*t9.*t11.*t12.*4.704e4+t2.*t4.*t7.*t10.*t11.*t13.*5.988e4+t3.*t4.*t5.*t10.*t11.*t14.*1.008e4+t2.*t3.*t4.*t8.*t9.*t10.*t11.*4.2e6+t2.*t4.*t8.*t9.*t10.*t11.*t14.*4.956e4-t3.*t4.*t8.*t9.*t10.*t11.*t13.*1.176e4).*(-1.0./1.0e2);t20.*t38.*(-1.0./1.0e1);t20.*(t7.*4.2e6-t5.*t7.*4.2e6-t6.*t7.*4.2e6-t2.*t13.*5.782e4+t3.*t14.*1.008e4+t2.*t7.*t8.*2.23668e5+t5.*t6.*t7.*4.2e6+t2.*t6.*t13.*5.782e4-t2.*t7.*t13.*5.988e4-t3.*t5.*t14.*1.008e4-t3.*t6.*t14.*1.008e4-t2.*t3.*t8.*t9.*4.2e6-t2.*t6.*t7.*t8.*2.23668e5+t2.*t6.*t7.*t13.*5.988e4+t3.*t5.*t6.*t14.*1.008e4-t4.*t7.*t10.*t11.*2.943e3-t2.*t8.*t9.*t14.*4.956e4+t3.*t8.*t9.*t13.*1.176e4+t4.*t7.*t10.*t12.*3.0006e5+t2.*t3.*t6.*t8.*t9.*4.2e6+t4.*t5.*t7.*t10.*t11.*2.943e3-t4.*t5.*t7.*t10.*t12.*3.0006e5+t2.*t6.*t8.*t9.*t14.*4.956e4-t3.*t6.*t8.*t9.*t13.*1.176e4).*(1.0./2.5e1);t20.*(t5.*t7.*2.23668e5+t2.*t7.*t8.*4.2e6-t3.*t5.*t9.*4.2e6-t5.*t6.*t7.*2.23668e5-t5.*t9.*t14.*4.956e4+t7.*t8.*t13.*2.03379e6-t2.*t6.*t7.*t8.*4.2e6+t3.*t5.*t6.*t9.*4.2e6+t2.*t3.*t8.*t14.*1.008e4-t2.*t3.*t9.*t13.*1.176e4+t5.*t6.*t9.*t14.*4.956e4-t6.*t7.*t8.*t13.*8.34e4-t2.*t3.*t6.*t8.*t14.*1.008e4+t2.*t3.*t6.*t9.*t13.*1.176e4-t2.*t4.*t7.*t8.*t10.*t11.*2.943e3+t2.*t4.*t7.*t8.*t10.*t12.*3.0006e5).*(-1.0./2.5e1);t20.*t38.*(-1.0./1.0e1);t20.*(t7.*5.88e7+t39+t40+t41-t5.*t6.*2.94e7-t5.*t7.*8.82e7-t6.*t7.*5.88e7-t2.*t13.*7.2716e5-t3.*t14.*1.206162e7+t2.*t7.*t8.*3.131352e6-t3.*t5.*t9.*1.565676e6+t5.*t6.*t7.*8.82e7+t2.*t6.*t13.*7.2716e5-t2.*t7.*t13.*9.2064e5+t3.*t5.*t14.*1.135542e7+t3.*t6.*t14.*3.5928e5-t2.*t3.*t8.*t9.*8.82e7-t2.*t6.*t7.*t8.*3.131352e6+t2.*t6.*t7.*t13.*9.2064e5+t3.*t5.*t6.*t14.*3.4692e5-t4.*t7.*t10.*t11.*4.1202e4-t2.*t8.*t9.*t14.*6.2328e5-t3.*t8.*t9.*t13.*1.407189e7+t4.*t7.*t10.*t12.*4.20084e6+t2.*t3.*t6.*t8.*t9.*8.82e7+t4.*t5.*t7.*t10.*t11.*4.1202e4-t4.*t5.*t7.*t10.*t12.*4.20084e6+t2.*t6.*t8.*t9.*t14.*6.2328e5+t3.*t6.*t8.*t9.*t13.*4.1916e5-t2.*t3.*t4.*t8.*t9.*t10.*t12.*2.10042e6).*(1.0./2.0e2);t20.*(t5.*t7.*1.565676e6+t2.*t7.*t8.*2.94e7-t3.*t5.*t9.*2.94e7-t5.*t6.*t7.*1.565676e6-t5.*t9.*t14.*3.4692e5+t7.*t8.*t13.*1.423653e7-t2.*t6.*t7.*t8.*2.94e7+t3.*t5.*t6.*t9.*2.94e7+t2.*t3.*t8.*t14.*7.056e4-t2.*t3.*t9.*t13.*8.232e4+t5.*t6.*t9.*t14.*3.4692e5-t6.*t7.*t8.*t13.*5.838e5-t2.*t3.*t6.*t8.*t14.*7.056e4+t2.*t3.*t6.*t9.*t13.*8.232e4-t2.*t4.*t7.*t8.*t10.*t11.*2.0601e4+t2.*t4.*t7.*t8.*t10.*t12.*2.10042e6).*(-1.0./2.0e2);t20.*(t3.*t5.*5.21892e5-t5.*t9.*9.8e6-t14.*t26.*1.176e4+t14.*t28.*2.03379e6+t2.*t3.*t8.*9.8e6-t3.*t5.*t6.*5.21892e5+t5.*t6.*t9.*9.8e6-t2.*t9.*t13.*2.744e4+t3.*t8.*t13.*4.74551e6-t2.*t3.*t6.*t8.*9.8e6+t2.*t6.*t8.*t14.*2.352e4+t2.*t6.*t9.*t13.*2.744e4-t3.*t6.*t8.*t13.*1.946e5+t2.*t7.*t8.*t14.*4.704e4-t3.*t5.*t9.*t14.*3.94782e6-t3.*t6.*t9.*t14.*1.668e5-t2.*t6.*t7.*t8.*t14.*4.704e4+t3.*t5.*t6.*t9.*t14.*4.704e4-t2.*t3.*t4.*t8.*t10.*t11.*6.867e3+t2.*t3.*t4.*t8.*t10.*t12.*7.0014e5).*(1.0./2.0e1);t20.*(t7.*2.94e7+t39+t40+t41-t5.*t6.*2.94e7-t5.*t7.*5.88e7-t6.*t7.*2.94e7-t2.*t13.*3.2242e5-t3.*t14.*1.213218e7+t2.*t7.*t8.*1.565676e6-t3.*t5.*t9.*1.565676e6+t5.*t6.*t7.*5.88e7+t2.*t6.*t13.*3.2242e5-t2.*t7.*t13.*5.0148e5+t3.*t5.*t14.*1.142598e7+t3.*t6.*t14.*4.2984e5-t2.*t3.*t8.*t9.*5.88e7-t2.*t6.*t7.*t8.*1.565676e6+t2.*t6.*t7.*t13.*5.0148e5+t3.*t5.*t6.*t14.*2.7636e5-t4.*t7.*t10.*t11.*2.0601e4-t2.*t8.*t9.*t14.*2.7636e5-t3.*t8.*t9.*t13.*1.415421e7+t4.*t7.*t10.*t12.*2.10042e6+t2.*t3.*t6.*t8.*t9.*5.88e7+t4.*t5.*t7.*t10.*t11.*2.0601e4-t4.*t5.*t7.*t10.*t12.*2.10042e6+t2.*t6.*t8.*t9.*t14.*2.7636e5+t3.*t6.*t8.*t9.*t13.*5.0148e5-t2.*t3.*t4.*t8.*t9.*t10.*t12.*2.10042e6).*(1.0./1.0e2)];
