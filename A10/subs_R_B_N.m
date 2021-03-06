function R_B_N = subs_R_B_N(q0,q1,q2,q3)
%SUBS_R_B_N
%    R_B_N = SUBS_R_B_N(Q0,Q1,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    26-Jun-2018 10:50:29

t2 = q1.*q2.*2.0;
t3 = t2-q0.*q3.*2.0;
t4 = q0.^2;
t5 = q1.^2;
t6 = q2.^2;
t7 = q3.^2;
t8 = q1.*q3.*2.0;
t9 = t8-q0.*q2.*2.0;
t10 = q2.*q3.*2.0;
t11 = t10-q0.*q1.*2.0;
R_B_N = reshape([t4+t5-t6-t7,t3,t9,t3,t4-t5+t6-t7,t11,t9,t11,t4-t5-t6+t7],[3,3]);
