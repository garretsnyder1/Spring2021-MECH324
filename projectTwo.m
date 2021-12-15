clc
clear
%% Linkage dependent inputs

a=3.175; %input
b=1.905; %coupler
c=5.08; %output
d=2.54; %ground
omega_2=30; %input angular velocity (deg/sec)
alpha_2=0; %input angular acceleration (deg/sec^2)
delta_1=257; %angle from L3 to coupler POF 1
delta_2=149; %angle from L3 to coupler POF 2
dist_1=1.016; %distance from L2 to coupler POF 1
dist_2=0.762; %distance from L2 to coupler POF 2


%% Calculating impossible theta values

n=1000;
theta_2=linspace(0,360,n);
diag=sqrt(d^2 + a^2 - 2*d*a*cosd(theta_2));
phi=acosd((b^2+diag.^2-c^2)./(2*b.*diag));

for i=1:n %remove impossible thetas
    check=phi(i);
    if ~isreal(check)
        theta_2(i)=0;
    end
end

theta_2=theta_2(theta_2~=0);

%% Theta 3 and 4

k1 = d/a;
k2 = d/c;
k3 = ((a^2)-(b^2)+(c^2)+(d^2))/(2*a*c);
k4 = d/b;
k5 = ((c^2)-(a^2)-(b^2)-(d^2))/(2*a*b);

Av = cosd(theta_2)-k1-k2*cosd(theta_2)+k3;
Bv = -2*sind(theta_2);
Cv = k1-(k2+1)*cosd(theta_2)+k3;
Dv = cosd(theta_2)-k1+k4*cosd(theta_2)+k5;
Ev = -2*sind(theta_2);
Fv = k1+(k4-1)*cosd(theta_2)+k5;

theta_3 = 2.*atand((-Ev-sqrt((Ev.^2)-4.*Dv.*Fv))./(2.*Dv));
theta_4 = 2.*atand((-Bv-sqrt((Bv.^2)-4.*Av.*Cv))./(2.*Av));

%% Graphing Full Position Analysis and Useful Position Analysis

figure
for i=1:1:4
    part=(floor(length(theta_2)/3));
    if i==1
        step=1;
    else
        step=part*(i-1);
    end
 
    curve1x=a*cosd(theta_2)+dist_1*cosd(theta_3+delta_1); %coupler curves
    curve1y=a*sind(theta_2)+dist_1*sind(theta_3+delta_1);
    curve2x=a*cosd(theta_2)+dist_2*cosd(theta_3+delta_2);
    curve2y=a*sind(theta_2)+dist_2*sind(theta_3+delta_2);
    
    coupler1x=a*cosd(theta_2(step))+dist_1*cosd(theta_3(step)+delta_1); %coupler locations
    coupler1y=a*sind(theta_2(step))+dist_1*sind(theta_3(step)+delta_1);
    coupler2x=a*cosd(theta_2(step))+dist_2*cosd(theta_3(step)+delta_2);
    coupler2y=a*sind(theta_2(step))+dist_2*sind(theta_3(step)+delta_2);
    
    ogX=[0,a*cosd(theta_2(step)),coupler1x,coupler2x,c*cosd(theta_4(step))+d,d]; %original positions
    ogY=[0,a*sind(theta_2(step)),coupler1y,coupler2y,c*sind(theta_4(step)),0];
    
    x2=a*cosd(theta_2); %line following path of point A
    y2=a*sind(theta_2);
    
    x3=c*cosd(theta_4)+d; %line following path of point B
    y3=c*sind(theta_4);
    
    subplot(2,2,i)
    hold on
    plot(ogX,ogY,'k-')
    plot(ogX(1),ogY(1),'ro')
    plot(ogX(2),ogY(2),'ro')
    plot(ogX(3),ogY(3),'go')
    plot(ogX(4),ogY(4),'go')
    plot(ogX(5),ogY(5),'ro')
    plot(ogX(6),ogY(6),'ro')
    plot(curve1x,curve1y,'g:','Linewidth',1.5)
    plot(curve2x,curve2y,'g:','Linewidth',1.5)
    plot(x2,y2,'b:','Linewidth',2)
    plot(x3,y3,'c:','Linewidth',2)
    xlabel('X'); ylabel('Y')
    sgtitle('Positions of Linkage at 4 Input Angles')
    hold off
end

figure
for i=1:1:4
    if i==1
        step=360;
    elseif i==2
        step=240;
    elseif i==3
        step=120;
    elseif i==4
        step=1;
    end
    
    coupler1x=a*cosd(theta_2(step))+dist_1*cosd(theta_3(step)+delta_1); %coupler locations
    coupler1y=a*sind(theta_2(step))+dist_1*sind(theta_3(step)+delta_1);
    coupler2x=a*cosd(theta_2(step))+dist_2*cosd(theta_3(step)+delta_2);
    coupler2y=a*sind(theta_2(step))+dist_2*sind(theta_3(step)+delta_2);
    
    ogX=[0,a*cosd(theta_2(step)),coupler1x,coupler2x,c*cosd(theta_4(step))+d,d]; %original positions
    ogY=[0,a*sind(theta_2(step)),coupler1y,coupler2y,c*sind(theta_4(step)),0];
    
    x2=a*cosd(theta_2); %line following path of point A
    y2=a*sind(theta_2);
    
    x3=c*cosd(theta_4)+d; %line following path of point B
    y3=c*sind(theta_4);
    
    subplot(2,2,i)
    hold on
    plot(ogX,ogY,'k-')
    plot(ogX(1),ogY(1),'ro')
    plot(ogX(2),ogY(2),'ro')
    plot(ogX(3),ogY(3),'go')
    plot(ogX(4),ogY(4),'go')
    plot(ogX(5),ogY(5),'ro')
    plot(ogX(6),ogY(6),'ro')
    plot(curve1x,curve1y,'g:','Linewidth',1.5)
    plot(curve2x,curve2y,'g:','Linewidth',1.5)
    plot(x2,y2,'b:','Linewidth',2)
    plot(x3,y3,'c:','Linewidth',2)
    xlabel('X'); ylabel('Y')
    sgtitle('Summary of Useful Linkage Positions (In Order)')
    hold off
end

%% Graphing Velocity Analysis

figure

omega_3 = ((a*omega_2)/(b)).*((sind(theta_4-theta_2))./(sind(theta_3-theta_4))); %calculating omega 3 and 4
omega_4 = ((a*omega_2)/(c)).*((sind(theta_2-theta_3))./(sind(theta_4-theta_3)));

omega_2plot=ones(1,630)*omega_2; %making plotable arrays
vA=omega_2plot*a;
vB=omega_4*c;

vC1 = dist_1.*omega_3.*(-sind(theta_3+delta_1)+1i.*cosd(theta_3+delta_1));
vC2 = dist_2.*omega_3.*(-sind(theta_3+delta_2)+1i.*cosd(theta_3+delta_2));
vC1x=real(vC1);
vC1y=imag(vC1);
vC2x=real(vC2);
vC2y=imag(vC2);
couplerVel_1=sqrt(vC1x.^2+vC1y.^2);
couplerVel_2=sqrt(vC2x.^2+vC2y.^2);

hold on
plot(theta_2,omega_3,'g-')
plot(theta_2,omega_4,'c-')
plot(theta_2,vA,'r-')
plot(theta_2,vB,'b-')
plot(theta_2,couplerVel_1,'m-')
plot(theta_2,couplerVel_2,'k-')
legend('\omega3','\omega4','Linear Velocity of A','Linear Velocity of B','Linear Velocity of Coupler Point 1','Linear Velocity of Coupler Point 2')
xlabel('\theta (deg)'); ylabel('Velocity (deg/s)')
hold off

sgtitle('Link Angular Velocity as \theta Changes')

%% Graphing Acceleration Analysis

figure

A_a = c.*sind(theta_4);
B_a = b.*sind(theta_3);
C_a = a.*alpha_2.*sind(theta_2)+a.*(omega_2^2).*cosd(theta_2)+b.*(omega_3.^2).*cosd(theta_3)-c.*(omega_4.^2).*cosd(theta_4);
D_a = c.*cosd(theta_4);
E_a = b.*cosd(theta_3);
F_a = a.*alpha_2.*cosd(theta_2)-a.*(omega_2^2)*sind(theta_2)-b.*(omega_3.^2).*sind(theta_3)+c.*(omega_4.^2).*sind(theta_4);

alpha_3 = (C_a.*D_a-A_a.*F_a)./(A_a.*E_a-B_a.*D_a);
alpha_4 = (C_a.*E_a-B_a.*F_a)./(A_a.*E_a-B_a.*D_a);

alpha_2plot=ones(1,630)*alpha_2; %making plotable arrays
aAC=omega_2plot.^2.*a;
aAT=alpha_2plot.*a;
aBC=omega_4.^2.*c;
aBT=alpha_4.*c;
aA=sqrt(aAC.^2+aAT.^2);
aB=sqrt(aBC.^2+aBT.^2);

hold on
plot(theta_2,alpha_3,'g-')
plot(theta_2,alpha_4,'c-')
plot(theta_2,aA,'r-')
plot(theta_2,aB,'b-')
legend('\alpha3','\alpha4','Linear Acceleration of A','Linear Acceleration of B')
xlabel('\theta (deg)'); ylabel('Acceleration (deg/s^2)')
hold off
sgtitle('Link Angular Acceleration as \theta Changes')
