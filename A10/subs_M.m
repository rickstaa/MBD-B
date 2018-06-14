function M = subs_M(k,l_0,b,x,y,z,q0,q1,q2,q3,x_d,y_d,z_d)
%SUBS_M
%    M = SUBS_M(K,L_0,B,X,Y,Z,Q0,Q1,Q2,Q3,X_D,Y_D,Z_D)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    14-Jun-2018 12:02:49

t3 = q0.^2;
t4 = q1.^2;
t5 = q2.^2;
t6 = q3.^2;
t9 = q0.*q2.*(1.0./5.0);
t10 = q0.*q3.*(9.9e1./5.0e1);
t11 = q1.*q2.*(9.9e1./5.0e1);
t12 = q1.*q3.*(1.0./5.0);
t13 = t3.*(1.0./1.0e2);
t14 = t4.*(1.0./1.0e2);
t15 = t5.*(1.0./1.0e2);
t16 = t6.*(1.0./1.0e2);
t2 = -t9-t10+t11+t12+t13+t14-t15-t16+x;
t18 = q0.*q1.*(1.0./5.0);
t19 = q0.*q3.*(1.0./5.0e1);
t20 = q1.*q2.*(1.0./5.0e1);
t21 = q2.*q3.*(1.0./5.0);
t22 = t3.*(9.9e1./1.0e2);
t23 = t4.*(9.9e1./1.0e2);
t24 = t5.*(9.9e1./1.0e2);
t25 = t6.*(9.9e1./1.0e2);
t7 = -t18-t19+t20+t21+t22-t23+t24-t25+y-9.0;
t27 = q0.*q1.*(9.9e1./5.0e1);
t28 = q0.*q2.*(1.0./5.0e1);
t29 = q1.*q3.*(1.0./5.0e1);
t30 = q2.*q3.*(9.9e1./5.0e1);
t31 = t3.*(1.0./1.0e1);
t32 = t4.*(1.0./1.0e1);
t33 = t5.*(1.0./1.0e1);
t34 = t6.*(1.0./1.0e1);
t8 = -t27-t28+t29+t30+t31-t32-t33+t34+z-2.5e1;
t17 = t2.^2;
t26 = t7.^2;
t35 = t8.^2;
t36 = t17+t26+t35;
t37 = 1.0./sqrt(t36);
t38 = y.*2.0;
t39 = q1.*q2.*(1.0./2.5e1);
t40 = q2.*q3.*(2.0./5.0);
t41 = t3.*(9.9e1./5.0e1);
t42 = t4.*(9.9e1./5.0e1);
t43 = t5.*(9.9e1./5.0e1);
t44 = t6.*(9.9e1./5.0e1);
t46 = q0.*q1.*(2.0./5.0);
t47 = q0.*q3.*(1.0./2.5e1);
t45 = t38+t39+t40+t41-t42+t43-t44-t46-t47-1.8e1;
t48 = conj(q0);
t49 = conj(q1);
t50 = conj(q2);
t51 = conj(q3);
t59 = q0.*q3.*(1.01e2./5.0e1);
t60 = q1.*q2.*(1.01e2./5.0e1);
t52 = -t9+t12+t13+t14-t15-t16+t59-t60+x;
t62 = t3.*(1.01e2./1.0e2);
t63 = t4.*(1.01e2./1.0e2);
t64 = t5.*(1.01e2./1.0e2);
t65 = t6.*(1.01e2./1.0e2);
t53 = -t18-t19+t20+t21-t62+t63-t64+t65+y+9.0;
t67 = q0.*q1.*(1.01e2./5.0e1);
t68 = q2.*q3.*(1.01e2./5.0e1);
t54 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t55 = x.*2.0;
t56 = q1.*q3.*(2.0./5.0);
t57 = t3.*(1.0./5.0e1);
t58 = t4.*(1.0./5.0e1);
t61 = t52.^2;
t66 = t53.^2;
t69 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t70 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t71 = z.*2.0;
t72 = q1.*q3.*(1.0./2.5e1);
t73 = t3.*(1.0./5.0);
t74 = t4.*(1.0./5.0);
t75 = t5.*(1.0./5.0);
t76 = t6.*(1.0./5.0);
t77 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t78 = t4.*(1.01e2./5.0e1);
t79 = t6.*(1.01e2./5.0e1);
t83 = t3.*(1.01e2./5.0e1);
t84 = t5.*(1.01e2./5.0e1);
t80 = t38+t39+t40-t46-t47+t78+t79-t83-t84+1.8e1;
t81 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t82 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t85 = t48.^2;
t86 = t49.^2;
t87 = t50.^2;
t88 = t51.^2;
t89 = t85-t86+t87-t88;
t90 = q1.*q2.*(9.9e1./2.5e1);
t103 = q0.*q2.*(2.0./5.0);
t105 = t5.*(1.0./5.0e1);
t106 = t6.*(1.0./5.0e1);
t116 = q0.*q3.*(9.9e1./2.5e1);
t91 = t55+t56+t57+t58+t90-t103-t105-t106-t116;
t92 = b.*t37.*t91.*x_d.*(1.0./2.0);
t93 = b.*t37.*t45.*y_d.*(1.0./2.0);
t94 = q2.*q3.*(9.9e1./2.5e1);
t98 = q0.*q1.*(9.9e1./2.5e1);
t99 = q0.*q2.*(1.0./2.5e1);
t95 = t71+t72+t73-t74-t75+t76+t94-t98-t99-5.0e1;
t96 = b.*t37.*t95.*z_d.*(1.0./2.0);
t97 = t92+t93+t96;
t100 = sqrt(t36);
t101 = l_0-t100;
t102 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t104 = q0.*q3.*(1.01e2./2.5e1);
t107 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t108 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t109 = q0.*q1.*(1.01e2./2.5e1);
t110 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t114 = q2.*q3.*(1.01e2./2.5e1);
t111 = t71+t72+t73-t74-t75+t76-t99+t109-t114-5.0e1;
t112 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t113 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t115 = t85-t86-t87+t88;
t117 = t37.*t91.*t97.*(1.0./2.0);
t156 = k.*t37.*t91.*t101.*(1.0./2.0);
t118 = t117-t156;
t119 = t48.*t50.*2.0;
t171 = t49.*t51.*2.0;
t120 = t119-t171;
t121 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t126 = q1.*q2.*(1.01e2./2.5e1);
t122 = t55+t56+t57+t58-t103+t104-t105-t106-t126;
t123 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t124 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t125 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t127 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t128 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t129 = t48.*t51.*2.0;
t179 = t49.*t50.*2.0;
t130 = t129-t179;
t131 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t132 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t133 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t134 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t135 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t136 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t137 = t37.*t45.*t97.*(1.0./2.0);
t178 = k.*t37.*t45.*t101.*(1.0./2.0);
t138 = t137-t178;
t139 = t48.*t49.*2.0;
t147 = t50.*t51.*2.0;
t140 = t139-t147;
t141 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t142 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t143 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t144 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t145 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t146 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t148 = t37.*t95.*t97.*(1.0./2.0);
t164 = k.*t37.*t95.*t101.*(1.0./2.0);
t149 = t148-t164;
t150 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t151 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t152 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t153 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t154 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t155 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t157 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t158 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t159 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t160 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t161 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t162 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t163 = t85+t86-t87-t88;
t165 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t166 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t167 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t168 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t169 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t170 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t172 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t173 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t174 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t175 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t176 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t177 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t180 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t181 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t182 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t183 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t184 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t185 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t186 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t187 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t188 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t189 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t190 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t191 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t192 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t193 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t194 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t195 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t196 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t197 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t198 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t199 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t200 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t201 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t202 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t203 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t204 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t205 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t206 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t207 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t208 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t209 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t210 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t211 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t212 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t213 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t214 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t215 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t216 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t217 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t218 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t219 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t220 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t221 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t222 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t223 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t224 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t225 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t226 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t227 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t228 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t229 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t230 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t231 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t232 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
t233 = t28-t29-t31+t32+t33-t34-t67+t68-z+2.5e1;
M = [t89.*t138.*(-1.0./1.0e1)-t118.*t120.*(9.9e1./1.0e2)+t118.*t130.*(1.0./1.0e1)+t115.*t149.*(9.9e1./1.0e2)-t138.*t140.*(9.9e1./1.0e2)+t140.*t149.*(1.0./1.0e1)-t89.*(t80.*1.0./sqrt(t61+t66+t77.^2).*(b.*t80.*y_d.*1.0./sqrt(t61+t66+t69.^2).*(1.0./2.0)+b.*x_d.*1.0./sqrt(t61+t66+t54.^2).*(t5.*(-1.0./5.0e1)-t6.*(1.0./5.0e1)+t55+t56+t57+t58+t104-q0.*q2.*(2.0./5.0)-q1.*q2.*(1.01e2./2.5e1)).*(1.0./2.0)+b.*z_d.*1.0./sqrt(t61+t66+t70.^2).*(t71+t72+t73-t74-t75+t76+t109-q0.*q2.*(1.0./2.5e1)-q2.*q3.*(1.01e2./2.5e1)-5.0e1).*(1.0./2.0)).*(1.0./2.0)-k.*t80.*(l_0-sqrt(t61+t66+t81.^2)).*1.0./sqrt(t61+t66+t82.^2).*(1.0./2.0)).*(1.0./1.0e1)-t115.*(t111.*1.0./sqrt(t61+t66+t110.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t102.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t107.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t108.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t111.*(l_0-sqrt(t61+t66+t112.^2)).*1.0./sqrt(t61+t66+t113.^2).*(1.0./2.0)).*(1.01e2./1.0e2)+t120.*(t122.*1.0./sqrt(t61+t66+t125.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t121.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t123.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t124.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t122.*(l_0-sqrt(t61+t66+t127.^2)).*1.0./sqrt(t61+t66+t128.^2).*(1.0./2.0)).*(1.01e2./1.0e2)+t140.*(t80.*1.0./sqrt(t61+t66+t144.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t141.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t142.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t143.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t80.*(l_0-sqrt(t61+t66+t145.^2)).*1.0./sqrt(t61+t66+t146.^2).*(1.0./2.0)).*(1.01e2./1.0e2)+t130.*(t122.*1.0./sqrt(t61+t66+t134.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t131.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t132.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t133.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t122.*(l_0-sqrt(t61+t66+t135.^2)).*1.0./sqrt(t61+t66+t136.^2).*(1.0./2.0)).*(1.0./1.0e1)+t140.*(t111.*1.0./sqrt(t61+t66+t153.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t150.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t151.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t152.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t111.*(l_0-sqrt(t61+t66+t154.^2)).*1.0./sqrt(t61+t66+t155.^2).*(1.0./2.0)).*(1.0./1.0e1);t118.*t120.*(1.0./1.0e2)-t115.*t149.*(1.0./1.0e2)-t130.*t138.*(1.0./1.0e1)-t120.*t149.*(1.0./1.0e1)+t138.*t140.*(1.0./1.0e2)+t118.*t163.*(1.0./1.0e1)-t115.*(t111.*1.0./sqrt(t61+t66+t168.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t165.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t166.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t167.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t111.*(l_0-sqrt(t61+t66+t169.^2)).*1.0./sqrt(t61+t66+t170.^2).*(1.0./2.0)).*(1.0./1.0e2)+t163.*(t122.*1.0./sqrt(t61+t66+t160.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t157.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t158.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t159.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t122.*(l_0-sqrt(t61+t66+t161.^2)).*1.0./sqrt(t61+t66+t162.^2).*(1.0./2.0)).*(1.0./1.0e1)+t140.*(t80.*1.0./sqrt(t61+t66+t183.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t180.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t181.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t182.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t80.*(l_0-sqrt(t61+t66+t184.^2)).*1.0./sqrt(t61+t66+t185.^2).*(1.0./2.0)).*(1.0./1.0e2)+t120.*(t122.*1.0./sqrt(t61+t66+t175.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t172.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t173.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t174.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t122.*(l_0-sqrt(t61+t66+t176.^2)).*1.0./sqrt(t61+t66+t177.^2).*(1.0./2.0)).*(1.0./1.0e2)-t130.*(t80.*1.0./sqrt(t61+t66+t189.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t186.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t187.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t188.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t80.*(l_0-sqrt(t61+t66+t190.^2)).*1.0./sqrt(t61+t66+t191.^2).*(1.0./2.0)).*(1.0./1.0e1)-t120.*(t111.*1.0./sqrt(t61+t66+t195.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t192.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t193.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t194.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t111.*(l_0-sqrt(t61+t66+t196.^2)).*1.0./sqrt(t61+t66+t197.^2).*(1.0./2.0)).*(1.0./1.0e1);t89.*t138.*(1.0./1.0e2)-t118.*t130.*(1.0./1.0e2)+t130.*t138.*(9.9e1./1.0e2)+t120.*t149.*(9.9e1./1.0e2)-t118.*t163.*(9.9e1./1.0e2)-t140.*t149.*(1.0./1.0e2)+t89.*(t80.*1.0./sqrt(t61+t66+t207.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t204.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t205.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t206.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t80.*(l_0-sqrt(t61+t66+t208.^2)).*1.0./sqrt(t61+t66+t209.^2).*(1.0./2.0)).*(1.0./1.0e2)-t130.*(t80.*1.0./sqrt(t61+t66+t219.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t216.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t217.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t218.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t80.*(l_0-sqrt(t61+t66+t220.^2)).*1.0./sqrt(t61+t66+t221.^2).*(1.0./2.0)).*(1.01e2./1.0e2)+t163.*(t122.*1.0./sqrt(t61+t66+t201.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t198.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t199.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t200.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t122.*(l_0-sqrt(t61+t66+t202.^2)).*1.0./sqrt(t61+t66+t203.^2).*(1.0./2.0)).*(1.01e2./1.0e2)-t130.*(t122.*1.0./sqrt(t61+t66+t213.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t210.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t211.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t212.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t122.*(l_0-sqrt(t61+t66+t214.^2)).*1.0./sqrt(t61+t66+t215.^2).*(1.0./2.0)).*(1.0./1.0e2)-t140.*(t111.*1.0./sqrt(t61+t66+t225.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t222.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t223.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t224.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t111.*(l_0-sqrt(t61+t66+t226.^2)).*1.0./sqrt(t61+t66+t227.^2).*(1.0./2.0)).*(1.0./1.0e2)-t120.*(t111.*1.0./sqrt(t61+t66+t231.^2).*(b.*t122.*x_d.*1.0./sqrt(t61+t66+t228.^2).*(1.0./2.0)+b.*t80.*y_d.*1.0./sqrt(t61+t66+t229.^2).*(1.0./2.0)+b.*t111.*z_d.*1.0./sqrt(t61+t66+t230.^2).*(1.0./2.0)).*(1.0./2.0)-k.*t111.*(l_0-sqrt(t61+t66+t232.^2)).*1.0./sqrt(t61+t66+t233.^2).*(1.0./2.0)).*(1.01e2./1.0e2)];