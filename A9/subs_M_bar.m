function out1 = subs_M_bar(beta,gamma)
%SUBS_M_BAR
%    OUT1 = SUBS_M_BAR(BETA,GAMMA)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    26-Jun-2018 20:58:22

t2 = cos(gamma);
t3 = beta.*2.0;
t4 = cos(t3);
t5 = sin(gamma);
t6 = t2.^2;
t7 = sin(beta);
t8 = t5.*2.0;
t9 = t8-3.0;
t10 = cos(beta);
t11 = t5.*3.0;
t12 = t11-2.0;
out1 = reshape([t4.*(2.1e1./1.0e2)-t5.*(9.0./5.0e1)+t6.*(3.0./5.0e1)-t4.*t5.*(9.0./5.0e1)-t4.*t6.*(3.0./5.0e1)+2.1e1./1.0e2,t2.*t7.*t9.*(-3.0./5.0e1),t10.*t12.*(-3.0./5.0e1),t2.*t7.*t9.*(-3.0./5.0e1),t5.*(-9.0./2.5e1)+t5.^2.*(3.0./2.5e1)+3.0./1.0e1,0.0,t10.*t12.*(-3.0./5.0e1),0.0,3.0./2.5e1],[3,3]);
