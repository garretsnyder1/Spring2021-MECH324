% This is the main function for generating a 4 bar linkage coupler curve X
% You (the student) need to change the coupler point location
% and how fast you ran the motor during the experiment (omega_2_rpm) X

clear all
close all
clc

%% user inputs
omega_2_rpm = 100;
alpha_2 = 0;
p = 0.0864;
a = 0.1283;
b = 0.5350;
c = 0.2775;
d = 0.5;
% units of meters for distance and degrees for the angle.

conversion_factor = 0.0254;
%Change the above to 0.08333 to convert everything to feet

% distance and angle from input link to coupler point
p = 4.8635 * conversion_factor; %CHANGEME, give value in inches but keep the conv. factor
delta_3_deg = 19.9812; %CHANGEME, degrees

%% defining constants
theta_2_deg = linspace(0,360);
%convert degrees to rads (you will want to work in radians)
% notice the .* you will want to use this for all scalar-matrix
% multipication (.*) division (./) and powers (.^)
theta_2 = theta_2_deg.*(pi/180);
%convert rpm to rad/s
omega_2 = omega_2_rpm*(2*pi/60);
%convert deg to rad
delta_3 = delta_3_deg*(pi/180);

%link lengths (be aware of your units)
a = 1.5543 * conversion_factor; %driver
b = 6.36 * conversion_factor; %coupler
c = 3.2955 * conversion_factor; %output rocker
d = 6 * conversion_factor;    %ground

%% Test function
test_equations();

%% Equation function - open equation_function.m to address all FIXME's until everything has passed
[k1, k2, k3, k4, k5, A_v, B_v, C_v, D_v, E_v, F_v, theta_3, theta_4, omega_3, omega_4,...
    Ax, Ay, Px, Py, V_Px, V_Py, A_a, B_a, C_a, D_a, E_a, F_a, A_Px, A_Py,...
    Prot, V_Prot, A_Prot] = equation_function(...
    alpha_2, p,a,b,c,d, theta_2, omega_2, delta_3);
 

% plotting results, comment/uncomment all the the lines below to plot/not plot everything

%Unrotated results

% figure
% subplot(2,2,1)
% plot(theta_2_deg, V_Px)
% hold on
% plot(theta_2_deg, V_Py)
% title('Velocity [m/sec]')
% xlabel('\theta')
% legend('X','Y,')
% 
% subplot(2,2,2)
% plot(theta_2_deg, A_Px)
% hold on
% plot(theta_2_deg, A_Py)
% title('Acceleration [m/s^2]')
% xlabel('\theta')
% legend('X','Y,')
% 
% subplot(2,2,3)
% plot(theta_2_deg,Px)
% hold on
% plot(theta_2_deg,Py)
% title('Displacement [m]')
% xlabel('\theta')
% legend('X','Y,')
% 
% subplot(2,2,4)
% plot(Px,Py)
% title('Trace of coupler point')
% xlabel('x [m]')
% ylabel('y [m]')
% 
% % create other plots as needed
% 
% 
% % Plot rotated results
% figure
% subplot(3,1,1)
% plot(theta_2_deg, V_Prot(1,:))
% hold on
% plot(theta_2_deg, V_Prot(2,:))
% title('Rotated Velocity')
% xlabel('\theta')
% 
% subplot(3,1,2)
% plot(theta_2_deg, A_Prot(1,:))
% hold on
% plot(theta_2_deg, A_Prot(2,:))
% title('Rotated Acceleration')
% xlabel('\theta')
% legend('X','Y,')
% 
% subplot(3,1,3)
% plot(theta_2_deg, Prot(1,:),theta_2_deg, Prot(2,:))
% title('Rotated Position')
% xlabel('\theta')
