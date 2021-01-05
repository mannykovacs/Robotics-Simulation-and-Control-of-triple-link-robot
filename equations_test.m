%% ESE 446-Task1
syms q1_dd q2_dd q3_dd g tau1 tau2 tau3
a0 = 0;
alpha0 = 0;
d1 = 0;
theta1 = q1;
 
T_01 = [cos(theta1) -sin(theta1) 0 a0;
        sin(theta1)*cos(alpha0) cos(theta1)*cos(alpha0) -sin(alpha0) -sin(alpha0)*d1;
        sin(theta1)*sin(alpha0) cos(theta1)*sin(alpha0) cos(alpha0) cos(alpha0)*d1;
        0 0 0 1];
 
alpha1 = 0;
a1 = L1;
d2 = 0;
theta2 = q2;
T_12 = [cos(theta2) -sin(theta2) 0 a1;
        sin(theta2)*cos(alpha1) cos(theta2)*cos(alpha1) -sin(alpha1) -sin(alpha0)*d2;
        sin(theta2)*sin(alpha1) cos(theta2)*sin(alpha1) cos(alpha1) cos(alpha0)*d2;
        0 0 0 1];   
  
alpha2 = 0;
a2 = L2;
d3 = 0;
theta3 = q3;
T_23 = [cos(theta3) -sin(theta3) 0 a2;
        sin(theta3)*cos(alpha2) cos(theta3)*cos(alpha2) -sin(alpha2) -sin(alpha2)*d3;
        sin(theta3)*sin(alpha2) cos(theta3)*sin(alpha2) cos(alpha2) cos(alpha2)*d3;
        0 0 0 1];  
 
alpha3 = 0;
a3 = L3;
d4 = 0;
theta4 = 0;
T_3E = [cos(theta4) -sin(theta4) 0 a3;
        sin(theta4)*cos(alpha3) cos(theta4)*cos(alpha3) -sin(alpha3) -sin(alpha3)*d4;
        sin(theta4)*sin(alpha3) cos(theta4)*sin(alpha3) cos(alpha3) cos(alpha3)*d4;
        0 0 0 1]  
T_02 = T_01*T_12;
T_03 = T_01*T_12*T_23;
T_0E = T_01*T_12*T_23*T_3E  % T_0E represents the representation of the end effector in the base frame.
 
 
P_01 = T_01(1:3,4);
P_12 = T_12(1:3,4); % P are the position or translation vectors of the robots.
P_23 = T_23(1:3,4);
 
R_01 = T_01(1:3,1:3); % Rotation matrices derived fropm the Transformation matrices for each frame
R_12 = T_12(1:3,1:3);
R_23 = T_23(1:3,1:3);
 
R_02 = R_01*R_12;
R_03 = R_01*R_12*R_23;
 
P_0C1 = T_01*[L1/2 0 0 1]'; % Position vector of the center of mass of each link
P_0C2 = T_02*[L2/2 0 0 1]';
P_0C3 = T_03*[L3/2 0 0 1]';
 
P_0C1 = P_0C1(1:3);
P_0C2 = P_0C2(1:3);
P_0C3 = P_0C3(1:3);
 
% the velocity Jacobian matrix found by differentiating position vector
% with respect to position of the link
 
Jv1 = [diff(P_0C1,q1) diff(P_0C1,q2) diff(P_0C1,q3)]; 
Jv2 = [diff(P_0C2,q1) diff(P_0C2,q2) diff(P_0C2,q3)];
Jv3 = [diff(P_0C3,q1) diff(P_0C3,q2) diff(P_0C3,q3)];
 
Jv1 = simplify(Jv1);
Jv2 = simplify(Jv2);
Jv3 = simplify(Jv3);
 
%%Angular Jacobian found using Z elements of the rotation matrix
Jw1 = [R_01(1,3) 0 0;
       R_01(2,3) 0 0;
       R_01(3,3) 0 0];
Jw2 = [R_01(1,3) R_02(1,3) 0;
       R_01(2,3) R_02(2,3) 0;
       R_01(3,3) R_02(3,3) 0];
Jw3 = [R_01(1,3) R_02(1,3) R_03(1,3);
       R_01(2,3) R_02(2,3) R_03(2,3);
       R_01(3,3) R_02(3,3) R_03(3,3)]; 
%Considering rotation around z to find the inertial properties.   
 
I_C1 = [0 0 0 ;
    0 0 0;
    0 0 I1];
I_C2 = [0 0 0 ;
    0 0 0;
    0 0 I2];
I_C3 = [0 0 0 ;
    0 0 0;
    0 0 I3];
 
mv1 = Jv1'*Jv1*m1;
mv2 = Jv2'*Jv2*m2;
mv3 = Jv3'*Jv3*m3;
 
mw1 = Jw1'*I_C1*Jw1;
mw2 = Jw2'*I_C2*Jw2;
mw3 = Jw3'*I_C3*Jw3;
 
M_q = mv1 + mv2 + mv3 + mw1+ mw2 + mw3;
%M_q = simplify(M_q) % Overall Mass Matrix
 
 
m_111 = diff(M_q(1,1), q1);
m_112 = diff(M_q(1,1), q2);
m_113 = diff (M_q(1,1), q3);
m_121 = diff (M_q(1,2), q1);
m_122 = diff (M_q(1,2), q2);
m_123 = diff (M_q(1,2), q3);
m_131 = diff (M_q(1,3), q1);
m_132 = diff (M_q(1,3), q2);
m_133 = diff (M_q(1,3), q3);
 
m_211 = diff(M_q(2,1), q1);
m_212 = diff(M_q(2,1), q2);
m_213 = diff (M_q(2,1), q3);
m_221 = diff (M_q(2,2), q1);
m_222 = diff (M_q(2,2), q2);
m_223 = diff (M_q(2,2), q3);
m_231 = diff (M_q(2,3), q1);
m_232 = diff (M_q(2,3), q2);
m_233 = diff (M_q(2,3), q3);
 
m_311 = diff(M_q(3,1), q1);
m_312 = diff(M_q(3,1), q2);
m_313 = diff (M_q(3,1), q3);
m_321 = diff (M_q(3,2), q1);
m_322 = diff (M_q(3,2), q2);
m_323 = diff (M_q(3,2), q3);
m_331 = diff (M_q(3,3), q1);
m_332 = diff (M_q(3,3), q2);
m_333 = diff (M_q(3,3), q3);
 
% centrifugal and corliolis matrices 
 
C_q = 1/2*[(m_111+m_111-m_111) (m_122+m_122-m_221) (m_133+m_133-m_331);
       (m_211+m_211-m_112) (m_222) (m_233+m_233-m_332);
       (m_311+m_311-m_113) (m_322+m_322-m_233) (m_333)];
   
B_q = [(m_112+m_121-m_121) (m_113+m_131-m_131) (m_123+m_132-m_231);
       (m_212+m_221-m_122) (m_213+m_231-m_132) (m_223+m_232-m_232);
       (m_312+m_321-m_123) (m_313+m_331-m_133) (m_323+m_332-m_233)];
   
C_q = simplify(C_q)
B_q = simplify(B_q)
 
V_q = C_q + B_q;
 
%% with values
L1 = 4;
L2 = 3;
L3 = 2;
m1 = 20;
m2 = 15;
m3 = 10;
I1 = 0.5;
I2 = 0.2;
I3 = 0.1;
q1 = 10;
q2 = 20;
q3 = 30;
% constructing gravity vector
 
G_q = -(Jv1'*[0; -m1*g; 0] + Jv2'*[0; -m2*g; 0] + Jv3'*[0; -m3*g; 0]);
G_q = simplify(G_q)
 
% Finding tau
 
%tau = [tau1; tau2; tau3] == M_q*[q1_dd; q2_dd; q3_dd] + V_q + G_q ;
 
% solving for q_dd (acceleration) 
 
%force_diff = solve(tau, q1_dd, q2_dd, q3_dd);
%q1_dd = simplify(force_diff.q1_dd);
%q2_pp = simplify(force_diff.q2_dd);
%q3_pp = simplify(force_diff.q3_dd);
 
%solution_tau = solve(tau, tau1, tau2, tau3);
%tau1 = simplify(solution_tau.tau1)
%tau2 = simplify(solution_tau.tau2)
%tau3 = simplify(solution_tau.tau3)
 
   


