function out1 = subs_qdp(phi1,phi2,phi1p,phi2p)
%SUBS_QDP
%    OUT1 = SUBS_QDP(PHI1,PHI2,PHI1P,PHI2P)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    25-Jun-2018 17:05:37

t3 = phi1-phi2;
t2 = cos(t3);
t4 = sin(t3);
t5 = t2.^2;
t6 = t5.*3.12694103621336e40;
t7 = t6-5.55900628660153e40;
t8 = 1.0./t7;
t9 = sin(phi2);
t10 = phi1p.^2;
t11 = sin(phi1);
t12 = phi2p.^2;
out1 = [t8.*(t11.*2.880850071868245e38-t2.*t9.*1.440425035934122e38+t4.*t12.*5.383851646372866e36+t2.*t4.*t10.*8.075777469559298e36).*3.872e3;t8.*(t9.*-3.84113342915766e38+t2.*t11.*4.321275107802367e38+t4.*t10.*2.153540658549146e37+t2.*t4.*t12.*8.075777469559298e36).*-3.872e3];
