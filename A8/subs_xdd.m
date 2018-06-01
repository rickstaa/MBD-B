function out1 = subs_xdd(phi1,phi2,x1d,y1d,phi1d,x2d,y2d,phi2d,t)
%SUBS_XDD
%    OUT1 = SUBS_XDD(PHI1,PHI2,X1D,Y1D,PHI1D,X2D,Y2D,PHI2D,T)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    01-Jun-2018 15:35:05

t2 = t.*pi;
t3 = cos(t2);
t4 = cos(phi1);
t5 = t4.^2;
t6 = sin(phi2);
t7 = sin(phi1);
t8 = cos(phi2);
t9 = t3.*t7.*(1.0./7.0);
t10 = t3.*t6.*(4.0./7.0);
t11 = phi1d.^2;
t12 = t8.^2;
t13 = phi2d.^2;
t14 = phi1.*2.0;
t15 = sin(t14);
t16 = phi2.*2.0;
t17 = sin(t16);
t18 = cos(t14);
out1 = [t9+t10-phi1d.*y1d.*(2.0./7.0)-t3.*t5.*t6.*(6.0./3.5e1)+phi1d.*t5.*y1d.*(2.0./7.0)+t3.*t4.*t7.*t8.*(6.0./3.5e1)-phi1d.*t4.*t7.*x1d.*(2.0./7.0);t3.*t4.*(-1.0./7.0)-t3.*t8.*(2.0./5.0)-t3.*t5.*t8.*(6.0./3.5e1)+phi1d.*t5.*x1d.*(2.0./7.0)-t3.*t4.*t6.*t7.*(6.0./3.5e1)+phi1d.*t4.*t7.*y1d.*(2.0./7.0);t3.*(-2.0./7.0)-t3.*t4.*t8.*(8.0./7.0)-t3.*t6.*t7.*(8.0./7.0)-phi1d.*t4.*x1d.*(1.0e1./7.0)-phi1d.*t7.*y1d.*(1.0e1./7.0);t9+t10-t4.*t11.*(1.0./4.0)-t8.*t13.*(1.0./8.0)+phi1d.*y1d.*(3.0./1.4e1)-phi2d.*y2d.*(1.0./2.0)-t3.*t5.*t6.*(1.3e1./3.5e1)+t3.*t6.*t12.*(1.3e1./3.5e1)+t3.*t7.*t12.*(1.0./7.0)-t4.*t11.*t12.*(1.0./4.0)+phi1d.*t15.*x1d.*(3.0./2.8e1)-phi2d.*t17.*x2d.*(1.0./4.0)-phi1d.*t5.*y1d.*(3.0./1.4e1)+phi1d.*t12.*y1d.*(3.0./1.4e1)+phi2d.*t12.*y2d.*(1.0./2.0)-t3.*t4.*t6.*t8.*(1.0./7.0)-t3.*t5.*t6.*t12.*(2.6e1./3.5e1)-t6.*t7.*t8.*t11.*(1.0./4.0)-phi1d.*t5.*t12.*y1d.*(3.0./1.4e1)+t3.*t4.*t7.*t8.*t12.*(2.6e1./3.5e1)-phi1d.*t5.*t6.*t8.*x1d.*(3.0./1.4e1)+phi1d.*t4.*t7.*t12.*x1d.*(3.0./1.4e1)-phi1d.*t4.*t6.*t7.*t8.*y1d.*(3.0./1.4e1);t3.*t4.*(-2.0./7.0)+t3.*t8.*(6.0./3.5e1)-t7.*t11.*(1.0./2.0)-t6.*t13.*(1.0./8.0)-t3.*t5.*t8.*(3.9e1./3.5e1)+t3.*t4.*t12.*(1.0./7.0)-t3.*t8.*t12.*(1.3e1./3.5e1)+t7.*t11.*t12.*(1.0./4.0)-phi1d.*t5.*x1d.*(3.0./7.0)+phi2d.*t12.*x2d.*(1.0./2.0)-phi1d.*t15.*y1d.*(3.0./1.4e1)+phi1d.*t17.*y1d.*(3.0./2.8e1)+phi2d.*t17.*y2d.*(1.0./4.0)-t3.*t4.*t6.*t7.*(2.6e1./3.5e1)+t3.*t6.*t7.*t8.*(1.0./7.0)+t3.*t5.*t8.*t12.*(2.6e1./3.5e1)-t4.*t6.*t8.*t11.*(1.0./4.0)+phi1d.*t5.*t12.*x1d.*(3.0./1.4e1)+t3.*t4.*t6.*t7.*t12.*(2.6e1./3.5e1)-phi1d.*t5.*t6.*t8.*y1d.*(3.0./1.4e1)+phi1d.*t4.*t7.*t12.*y1d.*(3.0./1.4e1)+phi1d.*t4.*t6.*t7.*t8.*x1d.*(3.0./1.4e1);t3.*(1.08e2./3.5e1)+t3.*t18.*cos(t16).*(5.2e1./3.5e1)+t3.*t4.*t8.*(8.0./7.0)+t3.*t6.*t7.*(8.0./7.0)-t4.*t6.*t11.*2.0+t7.*t8.*t11.*2.0+t3.*t15.*t17.*(5.2e1./3.5e1)+phi1d.*t8.*x1d.*(6.0./7.0)+phi2d.*t8.*x2d.*4.0+phi1d.*t6.*y1d.*(6.0./7.0)+phi2d.*t6.*y2d.*4.0+phi1d.*t6.*t15.*x1d.*(6.0./7.0)+phi1d.*t8.*t18.*x1d.*(6.0./7.0)+phi1d.*t8.*t15.*y1d.*(6.0./7.0)-phi1d.*t6.*t18.*y1d.*(6.0./7.0);t3.*t6.*(-2.0./5.0);t3.*t8.*(2.0./5.0);t3.*(1.0./7.0)+t3.*t4.*t8.*(6.0./3.5e1)+t3.*t6.*t7.*(6.0./3.5e1)-phi1d.*t4.*x1d.*(2.0./7.0)-phi1d.*t7.*y1d.*(2.0./7.0);t3.*(2.0./5.0)];