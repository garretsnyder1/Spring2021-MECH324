format compact
clear all; close all; clc
% Turn off warnings
warning('off');

%% define constants and inputs
omega_2_rpm = 100; %Motor RPM from experiment 
alpha_2 = 0; % rad/sec^2 

% link dimensions in meters
a = 3.175 * .01; %driver
b = 1.905 * .01; %coupler
c = 5.08 * .01; %output rocker
d = 2.54 * .01; %ground
w = .4 * .01; %link width

%distances and angles from attachment points to CENTERS OF MASS of all links
CG3 =  0.0222;  %meters
CG2 = 0.0191; %meters
CG4 = 0.015875; %meters
CG3_to_P = 0.01; %Dist from CG to point P on coupler (meters)

delta_3_deg = 9.51; %defined from theta 3 - it'll be close to 0
delta_2_deg = 18.91;  % defined from theta 2, should be along same line
delta_4_deg = 0;  % defined from theta 4, should be along same line
delta_CG_P_deg = 9.51; %defined from CG to point P on the coupler link
%NOTE: If a square, delta_CG_P_deg = delta_3_deg b/c CG of link 3 is in the center

%convert deg to rad
delta_3 = delta_3_deg*pi/180;
delta_2 = delta_2_deg*pi/180;
delta_4 = delta_4_deg*pi/180;
delta_CG_P = delta_CG_P_deg*pi/180;

%mass of links
m2 = .0016; %kg
m3 = .0023; %kg
m4 = .00139; %kg

% External forces on the fourbar
Fpx = 13.34; %Fp is force applied at point P % in newtons
Fpy = 0;
T4 = 0; %external torque acting on link 4 resisting motion (output link) % in nm

%% defining constants
%convert rpm to rad/s
omega_2 = omega_2_rpm*2*pi/60; % rad/s  

% Mass Moments of Inertia
IG2=1/12*m2*(a^2+w^2); %FIXME kg*m^2
IG3=1/12*m3*(d^2+d^2); %FIXME kg*m^2
IG4=1/12*m4*(c^2+w^2); %FIXME kg*m^2


%link length ratios
k1 = d/a;
k2 = d/c;
k3 = (a^2-b^2+c^2+d^2)/(2*a*c);
k4 = d/b;
k5 = (c^2-d^2-a^2-b^2)/(2*a*b);

%% Test equation_function.m
test_function();

%% for loop to calculate Matrix X for each degree of theta_2
for i=1:360 % iterate for all values of theta 2 
    theta_2_deg = i;
    %convert degrees to rads
    theta_2 = theta_2_deg*(pi/180);
    %A,B,C,D,E,F for velcity calculations
    A_v = cos(theta_2)-k1-k2*cos(theta_2)+k3;
    B_v = -2*sin(theta_2);
    C_v = k1-(k2+1)*cos(theta_2)+k3;
    D_v = cos(theta_2)-k1+k4*cos(theta_2)+k5;
    E_v = -2*sin(theta_2);
    F_v = k1+(k4-1)*cos(theta_2)+k5;

%% calculate theta_3, theta_4, omega_3,omega_4
    theta_3 = 2*atan((-E_v - sqrt(E_v^2-4*D_v*F_v))/(2*D_v));
    theta_4 = 2*atan((-B_v - sqrt(B_v^2-4*A_v*C_v))/(2*A_v));

    omega_3 = (a*omega_2/b)*sin(theta_4-theta_2)/sin(theta_3-theta_4);
    omega_4 = (a*omega_2/c)*sin(theta_2-theta_3)/sin(theta_4-theta_3);
    
%% calculate alpha_3 and alpha_4
    %A,B,C,D,E,F for accel calculations
    A_a = c.*sin(theta_4);
    B_a = b.*sin(theta_3);
    C_a = a.*alpha_2.*sin(theta_2)+a.*omega_2.^2.*cos(theta_2)+b.*omega_3.^2.*cos(theta_3)-c.*omega_4.^2.*cos(theta_4);
    D_a = c.*cos(theta_4);
    E_a = b.*cos(theta_3);
    F_a = a.*alpha_2.*cos(theta_2)-a.*omega_2.^2.*sin(theta_2)-b.*omega_3.^2.*sin(theta_3)+c.*omega_4.^2.*sin(theta_4);

    %calculate alpha_3 and alpha_4
    alpha_3 = (C_a.*D_a-A_a.*F_a)./(A_a.*E_a-B_a.*D_a); 
    alpha_4 = (C_a.*E_a-B_a.*F_a)./(A_a.*E_a-B_a.*D_a);
    
%% Equation function
% open the equation_function.m file to input the equations
    [R12x,R12y,R14x,R14y,R32x,R32y,R23x,R23y,R43x,R43y,R34x,R34y,...
        AG2x,AG2y,AG3x,AG3y,Rpx,Rpy,AG4x,AG4y] = ...
    equation_function10(a,b,c,CG2,CG3,CG4,CG3_to_P,theta_2,theta_3,theta_4,...
        delta_2,delta_3,delta_4,delta_CG_P,omega_2,omega_3,omega_4,alpha_2,alpha_3,alpha_4);

%% calculate matrix of reaction forces and torque (See eqn 11.9 for A and B)
    %create matrix of equation coefficients (left hand side values)
     Matrix_A= [1 0 1 0 0 0 0 0 0;
                0 1 0 1 0 0 0 0 0;
                -R12y R12x -R32y R32x 0 0 0 0 1;
                0 0 -1 0 1 0 0 0 0;
                0 0 0 -1 0 1 0 0 0;
                0 0 R23y -R23x -R43y R43x 0 0 0;
                0 0 0 0 -1 0 1 0 0;
                0 0 0 0 0 -1 0 1 0;
                0 0 0 0 R34y -R34x -R14y R14x 0];
    % create matrix of bounding conditions (right hand side of equations)
     Matrix_B= [m2*AG2x; 
                m2*AG2y; 
                IG2*alpha_2; 
                m3*AG3x-Fpx; 
                m3*AG3y-Fpy;
                IG3*alpha_3 - Rpx*Fpy + Rpy*Fpx; 
                m4*AG4x;
                m4*AG4y;
                IG4*alpha_4-T4];

     Matrix_X = linsolve(Matrix_A, Matrix_B);
     % values of Matrix_X= ([F12x F12y F32x F32y F43x F43y F14x F14y T12]);
     Matrix_of_results(:,i) = Matrix_X;
end

%find the maximum toruqe
max_T = max(Matrix_of_results(9,:))
min_T = min(Matrix_of_results(9,:))

%% convert input torque to input power required (Torque*Speed)
P_in = Matrix_of_results(9,:)*omega_2; %[W]
P_in_hp = P_in*0.00134102; %[hp]

%% plot results
thetas = 1:360;

% %%make plots for forces on each link and torque for all angles of theta
figure(5)
subplot(4,1,1)
plot(thetas,Matrix_of_results(9,:))
title('Necessary Torque on Input Link (Nm)')
 xlabel('Angle of Input Link (deg)')
ylabel('Torque')

subplot(4,1,2)
plot(thetas,Matrix_of_results(1,:),thetas,Matrix_of_results(2,:))
title('Force on Input Link (N)')
 xlabel('Angle of Input Link (deg)')
ylabel('Force')
legend('F12x','F12y')
inmax1=max(Matrix_of_results(1,50:300))
inmin1=min(Matrix_of_results(1,50:300))
inmax2=max(Matrix_of_results(2,50:300))
inmin2=min(Matrix_of_results(2,50:300))
subplot(4,1,3)
plot(thetas,Matrix_of_results(3,:),thetas,Matrix_of_results(4,:))
title('Force on Coupler (N)')
 xlabel('Angle of Input Link (deg)')
ylabel('Force')
legend('F32x','F32y')

subplot(4,1,4)
plot(thetas,Matrix_of_results(5,:),thetas,Matrix_of_results(6,:))
title('Force on Output Link (N)')
xlabel('Angle of Output Link (deg)')
ylabel('Force')
legend('F43x','F43y')

figure(6)
plot(thetas,P_in_hp)
title('Input Power Required')
xlabel('Angle of Input Link (deg)')
ylabel('Power [hp]')
