% six bar linkage
%Static equilibrium
clc;
clear;
% define the joints
A = [7 4 0];
B = [5 16 0];
C = [25 25 0];
D = [23 10 0];
E = [18 35 0];
F = [43 32 0];
G = [45 17 0];

%joint values

new_B_x(1) = B(1);
new_B_y(1) = B(2);
new_C_x(1) = C(1);
new_C_y(1) = C(2);
new_E_x(1) = E(1);
new_E_y(1) = E(2);
new_F_x(1) = F(1);
new_F_y(1) = F(2);

% define the lengths of the links
lAB = norm(B - A);
lBC = norm(C - B);
lCD = norm(D - C);
lDE = norm(E - D);
lEF = norm(F - E);
lFG = norm(G - F);
lCE = norm(E-C);
%weight of links
WAB = [0 -1 0];
WBEC = [0 -1 0];
WCD = [0 -1 0];
WED = [0 -1 0];
WEF = [0 -1 0];
WFG = [0 -1 0];
%center of mass of each link
S1 = (A+B)/2;
S2 = (B+C+E)/3;
S3 = (C+D)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;
syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin
ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
InputTorque = [0 0 Tin];
%Applied Force
AppliedForce = [50 0 0];
%Static Equilibrium Conditions for link AB
eqn1 = ForceA + ForceB + WAB == [0 0 0];
%Sum of moments = 0 with respect to CoM of Link AB
%S1A x FA + S1B x FB + InputTorque = 0
eqn2 = cross(A-S1, ForceA) + cross(B-S1, ForceB) + InputTorque == [0 0 0];
%Equations for Link BEC
%Sum of Force = 0
% -Fb + fc + fe + wbec = 0
eqn3 = -ForceB + ForceC + ForceE + WBEC == [0 0 0];
%SUm of Moments equals 0
eqn4 = cross(B-S2, ForceB) + cross(C-S2, ForceC) + cross(E-S2, ForceE) == [0 0 0];
%Equations for Link CD
%Static Equilibrium Conditions for link CD
%Sum of Force = 0
eqn5 = -ForceC + ForceD + WCD == [0 0 0];
%Sum of Moments = 0
%Sum of Moments = 0 for Link CD
eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == [0 0 0];
%Equations for Link EF
%Sum of Force = 0 for Link EF
eqn7 = -ForceE + ForceF + WEF == [0 0 0];
%Sum of Moments = 0 for Link EF
eqn8 = cross(E-S4, -ForceE) + cross(F-S4, ForceF) == [0 0 0];
%Equations for Link FG
%Sum of Force = 0 for Link FG
eqn9 = -ForceF + ForceG + WFG + AppliedForce == [0 0 0];
%Sum of Moments = 0 for Link FG
eqn10 = cross(F-S5, -ForceF) + cross(G-S5, ForceG) == [0 0 0];
%Solving the 10 equations
eqnMatrix = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];
staticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);
Force_Ax = double(staticSolution.FAx);
Force_Ay = double(staticSolution.FAy);
Force_By = double(staticSolution.FBy);
Force_Cx = double(staticSolution.FCx);
Force_Cy = double(staticSolution.FCy);
Force_Dx = double(staticSolution.FDx);
Force_Dy = double(staticSolution.FDy);
Force_Ex = double(staticSolution.FEx);
Force_Ey = double(staticSolution.FEy);
Force_Fx = double(staticSolution.FFx);
Force_Fy = double(staticSolution.FFy);
Force_Gx = double(staticSolution.FGx);
Force_Gy = double(staticSolution.FGy);
StaticTorqueIn = double(staticSolution.Tin);
% Angular velocities calculations
% Loop ABCDA
syms wBEC wCD
omega_AB = [0 0 1];
omega_BEC = [0 0 wBEC];
omega_CD = [0 0 wCD];
eqn11 = cross(omega_AB, B-A) + cross(omega_BEC, C-B) + cross(omega_CD, D-C) == [0 0 0];
loop1Solution = solve(eqn11, [wBEC wCD]);
angularVelocity_BEC = double(loop1Solution.wBEC);
angularVelocity_CD = double(loop1Solution.wCD);
% Second loop DCEFGD
syms wEF wFG
omega_AB = [0 0 1];
omega_BEC = [0 0 angularVelocity_BEC];
omega_CD = [0 0 angularVelocity_CD];
omega_EF = [0 0 wEF];
omega_FG = [0 0 wFG];
eqn12 = cross(omega_CD,C-D) + cross(omega_BEC, E-C) + cross(omega_EF, F-E) + cross(omega_FG, G-F) == [0 0 0];
loop1Solution = solve(eqn12, [wEF wFG]);
angularVelocity_EF = double(loop1Solution.wEF);
angularVelocity_FG = double(loop1Solution.wFG);
% Replace the symbolic omega_EF/omega_FG placeholders with their solved
% numeric values now that angularVelocity_EF/angularVelocity_FG are known
% (mirrors what was already done for omega_BEC and omega_CD above).
omega_EF = [0 0 angularVelocity_EF];
omega_FG = [0 0 angularVelocity_FG];
% Angular Acceleration
% Loop 1 ABCDA
syms aBEC aCD
alpha_AB = [0 0 0];
alpha_BEC = [0 0 aBEC];
alpha_CD = [0 0 aCD];
a_B_A = cross(alpha_AB, B-A) + cross(omega_AB, cross(omega_AB, B-A));
a_C_B = cross(alpha_BEC, C-B) + cross(omega_BEC, cross(omega_BEC, C-B));
a_D_C = cross(alpha_CD, D-C) + cross(omega_CD, cross(omega_CD, D-C));
eqn13 = a_B_A + a_C_B + a_D_C == [0 0 0];
loop1AccSolution = solve(eqn13, [aBEC aCD]);
alphaBEC = double(loop1AccSolution.aBEC);
alphaCD = double(loop1AccSolution.aCD);
% Angular Acceleration for Loop 2 DCEFGD
syms aEF aFG
alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];
angVel_EF = [0 0 angularVelocity_EF];
anVel_FG = [0 0 angularVelocity_FG];
alpha__BEC = [0 0 alphaBEC];
alpha__CD = [0 0 alphaCD];
a_C_D = cross(alpha__CD, C-D) + cross(omega_CD, cross(omega_CD, C-D));
a_E_C = cross(alpha__BEC, E-C) + cross(omega_BEC, cross(omega_BEC, E-C));
a_F_E = cross(alpha_EF, F-E) + cross(angVel_EF, cross(angVel_EF, F-E));
a_G_F = cross(alpha_FG, G-F) + cross(anVel_FG, cross(anVel_FG, G-F));
eqn14 = a_C_D + a_E_C + a_F_E + a_G_F == [0 0 0];
loop2AccSolution = solve(eqn14, [aEF, aFG]);
alphaFG = double(loop2AccSolution.aFG);
alphaEF = double(loop2AccSolution.aEF);
% Replace the symbolic alpha_EF/alpha_FG placeholders with their solved
% numeric values now that alphaEF/alphaFG are known (mirrors what was
% already done for alpha__BEC and alpha__CD above).
alpha_EF = [0 0 alphaEF];
alpha_FG = [0 0 alphaFG];
% Velocity at Joint B
vB_A = cross(omega_AB, B-A);
% VE_A = V_E_B + V_B_A Joint E
v_E_B = cross(omega_BEC, E-B);
vE_A = v_E_B + vB_A;
%V_S4/G = V_S4_F + V_F_G
V_S4_F = cross(omega_EF, S4-F);
V_F_G = cross(omega_FG, F-G);
vS4_G = V_S4_F + V_F_G;
% Velocity of S1
vS1_A = cross(omega_AB, S1 - A);
% Velocity of S2
vS2_D = cross(omega_CD, S2 - D);
%Velocity of S3
vS3_G = cross(omega_FG, S3 - D);
%Velocity of S4
% VS5_A = VB_A + VS5_B
VS5_B = cross(omega_BEC, S5 - B);
VS5_A = vB_A + VS5_B;
%Velocity at Joint C
VC_D = cross(omega_CD, C -D);
%Velocity at Joint F
VF_G = cross(omega_FG, F - G);
% Acceleration at Joint B
aB_A = cross(alpha_AB, B-A) + cross(omega_AB, cross(omega_AB, B-A));
% Acceleration at Joint E
% aE_A = aE_B + aB_A
aE_B = cross(alpha_BEC, E-B) + cross(omega_BEC, cross(omega_BEC, E-B));
aE_A = aE_B + aB_A;
% Acceleration at Joint C
aC_D = cross(alpha_CD, C-D) + cross(omega_CD, cross(omega_CD, C-D));
% Acceleration at Joint F
aF_G = cross(alpha_FG, F-G) + cross(omega_FG, cross(omega_FG, F-G));
% Acceleration at S1
aS1_A = cross(alpha_AB, S1-A) + cross(omega_AB, cross(omega_AB, S1 - A));
% Acceleration at S2
aS2_D = cross(alpha__CD, S2 - D) + cross(omega_CD, cross(omega_CD, S2 - D));
% Acceleration at S3
aS3_G = cross(alpha_FG, S3 - G) + cross(omega_FG, cross(omega_FG, S3 - G));
% Acceleration at S4
% aS4_G = aF_G + aS4_F
aS4_F = cross(alpha_EF, S4 - F) + cross(omega_EF, cross(omega_EF, S4 - F));
aS4_G = aF_G + aS4_F;
% Acceleration at S5
% aS5_A = aB_A + aS5_B
aS5_B = cross(alpha__BEC, S5 - B) + cross(omega_BEC, cross(omega_BEC, S5 - B));
aS5_A = aB_A + aS5_B;
% Newton's Second Law Implementation
MassAB = 1;
MassBEC = 1;
MassCD = 1;
MassEF = 1;
MassFG = 1;
% Mass Moment of Inertia
J_AB = 1;
J_BEC = 1;
J_CD = 1;
J_EF = 1;
J_FG = 1;
syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin
% Define Forces
NForceA = [NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD = [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF = [NFFx NFFy 0];
NForceG = [NFGx NFGy 0];
NInputTorque = [0 0 NTin];
% Equations for Link AB
% Sum of Forces
eqn15 = NForceA + NForceB + WAB == MassAB * aS1_A;
% sum of moments = 0
eqn16 = cross(A - S1, NForceA) + cross(B-S1, NForceB) + NInputTorque == J_AB * alpha_AB;
% Equations for Link BEC
% Sum of Forces
eqn17 = -NForceB + NForceC + NForceE + WBEC == MassBEC * aS2_D;
% Sum of Moments
eqn18 = cross(B - S2, -NForceB) + cross(C - S2, NForceC) + cross(E - S2, NForceE) == J_BEC * alpha__BEC;
% Equations for Link CD
% Sum of Forces
eqn19 = -NForceC + NForceD + WCD == MassFG * aS3_G;
% Sum of Moments
eqn20 = cross(C-S3, -NForceC) + cross(D-S3, NForceD) == J_CD * alpha__CD;
%Equations for link EF
%Sum of forces
eqn21 = -NForceE + NForceF + WEF == MassEF * aS4_G;
% Sum of Moments
eqn22 = cross(E-S4, -NForceE) + cross(F-S4, NForceF) == J_EF * [0 0 alphaEF];
%Equations for link FG
% Sum of Forces
eqn23 = -NForceF + NForceG + WFG + AppliedForce == MassFG * aS5_A;
% Sum of Moments
eqn24 = cross(F-S5, -NForceF) + cross(G - S5, NForceG) == J_FG * [0 0 alphaFG];
%Solving equations
NeqnMatrix = [eqn15, eqn16, eqn17, eqn18, eqn19, eqn20, eqn21, eqn22, eqn23, eqn24];
DynamicSolution = solve(NeqnMatrix, [NFAx, NFAy, NFBx, NFBy, NFCx, NFCy, NFDx, NFDy, NFEx, NFEy, NFFx, NFFy, NFGx, NFGy, NTin]);
%Extract the forces from the dynamic solution
NForce_Ax = double(DynamicSolution.NFAx);
NForce_Ay = double(DynamicSolution.NFAy);
NForce_Bx = double(DynamicSolution.NFBx);
NForce_By = double(DynamicSolution.NFBy);
NForce_Cx = double(DynamicSolution.NFCx);
NForce_Cy = double(DynamicSolution.NFCy);
NForce_Dx = double(DynamicSolution.NFDx);
NForce_Dy = double(DynamicSolution.NFDy);
NForce_Ex = double(DynamicSolution.NFEx);
NForce_Ey = double(DynamicSolution.NFEy);
NForce_Fx = double(DynamicSolution.NFFx);
NForce_Fy = double(DynamicSolution.NFFy);
NForce_Gx = double(DynamicSolution.NFGx);
NForce_Gy = double(DynamicSolution.NFGy);
NInputTorque = double(DynamicSolution.NTin);
%Circle intersections technique
%joint coordinates have been defined
% length of links also defined
%compute initial angle of input link AB
initial_theta = atan2(B(2) - A(2), B(1) - A(1));
if (initial_theta < 0)
inputAngle = 2*pi + initial_theta;
else
inputAngle = initial_theta;
end
for theta = 1:1:360
% New position of joint B
B_new = A + [lAB*cos(inputAngle + deg2rad(theta)) lAB*sin(inputAngle + deg2rad(theta)) 0];
% New position of C
% with B_new as center, BC as radius
% with D as center and DC as radius
[Cx, Cy] = circcirc(B_new(1), B_new(2), lBC, D(1), D(2), lCD);

%checking if there is a Nan

circIntersect_x = any(isnan(vpa(Cx)));
circIntersect_y = any(isnan(vpa(Cy)));

if circIntersect_x && circIntersect_y == 0 
    C_1 = [Cx(1) Cy(1) 0];
    C_2 = [Cx(2) Cy(2) 0];

    %distance to determine C_1 or C_2 is correct

    dist1 = norm(C_1-C);
    dist2 = norm(C_2-C);

    if(dist1<dist2)
        C_new = vpa(C_1);

    else 

        C_new = vpa(C_2);
    end

    %New position of joint E using B_new and C_new

    [Ex,Ey] = circcirc(B_new(1), B_new(2), lBE, C_new(1), C_new(2), lCE);

    %checking if there is a Nan

    circIntersect_x = any(isnan(vpa(Ex)));
    circIntersect_y = any(isnan(vpa(Ey)));

    if circIntersect_x && circIntersect_y == 0 
        E_1 = [Ex(1) Ey(1) 0];
        E_2 = [Ex(2) Ey(2) 0];

        %distance to determine E_1 or E_2 is correct

        dist1 = norm(E_1-E);
        dist2 = norm(E_2-E);

        if(dist1<dist2)
            E_new = vpa(E_1);

        else 

            E_new = vpa(E_2);
        end
    else
        fprintf('New position of E cannot be determined at angle: %d degree',theta);
    end
        %New position of joint F using E_new and G

        [Fx,Fy] = circcirc(E_new(1), E_new(2), lEF, G(1), G(2), lFG);

        %checking if there is a Nan

        circIntersect_x = any(isnan(vpa(Fx)));
        circIntersect_y = any(isnan(vpa(Fy)));

        if circIntersect_x && circIntersect_y == 0 
            F_1 = [Fx(1) Fy(1) 0];
            F_2 = [Fx(2) Fy(2) 0];

            %distance to determine F_1 or F_2 is correct

            dist1 = norm(F_1-F);
            dist2 = norm(F_2-F);

            if(dist1<dist2)
                F_new = vpa(F_1);

            else 

                F_new = vpa(F_2);
            end
        else
            fprintf('New position of F cannot be determined at angle: %d degree',theta);
        end

        %store values for plotting

        new_B_x(theta+1) = B_new(1);
        new_B_y(theta+1) = B_new(2);
        new_C_x(theta+1) = C_new(1);
        new_C_y(theta+1) = C_new(2);
        new_E_x(theta+1) = E_new(1);
        new_E_y(theta+1) = E_new(2);
        new_F_x(theta+1) = F_new(1);
        new_F_y(theta+1) = F_new(2);

        B=B_new;
        C=C_new;
        E=E_new;
        F=F_new;

else
    fprintf('New position of C cannot be determined at angle: %d degree',theta);
end
end