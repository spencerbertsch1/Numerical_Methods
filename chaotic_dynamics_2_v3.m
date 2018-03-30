%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL EXAM: PY331
%~~~~~~~~~~~~~~~~~~~
% Program #2: Double Pendulum and Poincare Map
%
% Author: Spencer Bertsch
% Date: May 9, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all variables
clc
close all
clear
%%%%%%%%%%%%%%%%%%%%

%Input from user: 
disp('Enter "1" to see a movie of the double pendulum solved with the Euler Method.')
disp('Enter "2" to see a movie of the double pendulum with lines showing the traced paths.')
disp('Enter "3" to see plots of the two positions and the two angles over time.')
disp('Enter "4" to see a movie of the paths traced out by the masses.')
disp('Enter "5" to see Potential, Kinetic, and Total energy of the inner mass over time.')
disp('Enter "6" to see Poincare Map_1: Theta1 vs dTheta1dt')
disp('Enter "7" to see Poincare Map_2: dTheta1dt vs dTheta21dt')
disp('Enter "8" to see a proof of concept output from the ODE45 differential equation solver.')
problemnumber = input('Enter a number: ');

% Declare Variables 

l1 = 2; %<----- SET LENGTH OF INNER ROD --1

l2 = 1.5; %<----- SET LENGTH OF OUTER ROD --1.4

m1 = 2.5; %<----- SET MASS 1 --2.5

m2 = 1; %<----- SET MASS 2 --1.6

g = 9.8; %<----- ACCELERATION DUE TO EARTH'S GRAVITY (m/sec^2) 

%Initial Conditions

u10 = 1.6; %<----- SET INITIAL ANGLE 1 --1.6

u20 = 0; %<----- SET INITIAL VELOCITY 1 -- 0

v10 = 2.5; %<----- SET INITIAL ANGLE 2 --2.2

v20 = 0; %<----- SET INITIAL VELOCITY 2 --0 

Pe_0 = abs((g)*(l1)*( ((u10).^2) /2)); %Initial Potential energy

Ke_0 = abs(0.5*(u20*l1)^2); %Initial Kinetic energy

Etotal_0 = Pe_0 + Ke_0; 

%time vector%
dt = 0.025; %Time Step (remember making time step smaller will reduce error)
%(and you can change the frame rate of the plot movies to account for the
%slower movie speed) 
tStop = 15; %Record Length
t = 0:dt:tStop; %Time Vector


%% Numerical Solution Euler Method 

u1 = zeros(1, length(t)); %u1 = theta1(t)
u1(1) = u10; %set initial theta1

u2 = zeros(1, length(t)); %u2 = theta1Prime(t), or angular velocity 1
u2(1) = u20; %set initial velocity of first mass

v1 = zeros(1, length(t)); %v1 = theta2(t)
v1(1) = v10; %set initial theta2

v2 = zeros(1, length(t)); %v1 = theta2(t)
v2(1) = v20; %set initial theta2

PeE = zeros(1,length(t)); % Initializes the delta x vector
PeE(1) = Pe_0; % Set initial delta x to 0

KeE = zeros(1,length(t)); % Initializes the delta x vector
KeE(1) = Ke_0; % Set initial delta x to 0

EtotalE = zeros(1,length(t)); % Initializes the delta x vector
EtotalE(1) = Etotal_0; % Set initial delta x to 0


%Euler Method
for i=2:length(t)
 
    a = (m1+m2)*l1 ;
    
    b = m2*l2*cos(u1(i-1)-v1(i-1)) ;
    
    c = m2*l1*cos(u1(i-1)-v1(i-1)) ;
    
    d = m2*l2 ;
    
    e = -m2*l2*v2(i-1)* v2(i-1)*sin( u1(i-1) - v1(i-1) )-g*(m1+m2)*sin(u1(i-1)) ;
    
    f = m2*l1*u2(i-1)*u2(i-1)*sin(u1(i-1) - v1(i-1))-m2*g*sin(v1(i-1)) ;
    
                           %%%%%%%%%%%%%%
    
    du2dt= (e*d-b*f)/(a*d-c*b) ;
    
    u2(i) = u2(i-1) + dt * (du2dt);
    
    dv2dt= (a*f-c*e)/(a*d-c*b) ;
    
    v2(i) = v2(i-1) + dt * (dv2dt);
    
    du1dt = (u2(i)); 
    
    u1(i) = u1(i-1) + dt * (du1dt);
    
    dv1dt = (v2(i)); 
    
    v1(i) = v1(i-1) + dt * (dv1dt);
    
                           %%%%%%%%%%%%%%
    
    %Calculate energy in the system as a function of time 
    PeE(i) = abs((g)*(l1)*( ((u1(i)).^2) /2)); %Potential energy
    KeE(i) = abs(0.5*(u2(i)*l1)^2); %Kinetic energy
    EtotalE(i) = PeE(i) + KeE(i); %total energy
   
end

% Now that we have the angular positions and angular velocities as functions
% of time, we can find the x and y positions relatively easily using
% trigonometry 
x1 = l1*sin(u1);
y1 = -l1*cos(u1);
x2 = l1*sin(u1)+l2*sin(v1);
y2 = -l1*cos(u1)-l2*cos(v1);

%% Numerical Solution 5th Order Runga Kutta using "ode45" 

t_stop=40;
y0=[1.6 0 2.2 0]; %[theta_1 dtheta_1/dt theta_2 dtheta_2/dt]
[tt,yy]=ode45(@solution, [0 ,t_stop],y0);

%Convert radial position into x and y coordinates
xx1=l1*sin(yy(:,1));
yy1=-l1*cos(yy(:,1));
xx2=l1*sin(yy(:,1))+l2*sin(yy(:,3));
yy2=-l1*cos(yy(:,1))-l2*cos(yy(:,3));

%% Plots
switch problemnumber 
    case 1 %represents the first case in the switch statement and the first plot

                            %%%%CASE 1%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[388   323   561   474]);
for pn = 1:2:length(t)
    figure(1)
    clf
    plot(0, 0,'black.','markersize',20);
    hold on
    plot(x1(pn),y1(pn),'r.','markersize',(20*m1)); %size of ball is proportional to its mass
    plot(x2(pn),y2(pn),'b.','markersize',(20*m2)); %size of ball is proportional to its mass
    axis square
    
    line([0 x1(pn)], [0 y1(pn)],'Linewidth',1.5);
    axis([-(l1+l2) l1+l2 -(l1+l2) l1+l2]);
    line([x1(pn) x2(pn)], [y1(pn) y2(pn)],'linewidth',1.5);
    h=gca; 
    get(h,'fontSize') ;
    set(h,'fontSize',12)
    xlabel('X','fontSize',12);
    ylabel('Y','fontSize',12);
    title('Double Pendulum','fontsize',14)
    pause(0.03)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


case 2 %second plot
                           %%%%CASE 2%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[839   249   423   423]);
for pn = 1:2:length(t)
    figure(1)
    clf
    plot(0, 0,'black.','markersize',20);
    hold on
    plot(x1(pn),y1(pn),'r.','markersize',(20*m1)); %size of ball is proportional to its mass
    plot(x2(pn),y2(pn),'b.','markersize',(20*m2)); %size of ball is proportional to its mass
    plot(x1(1:pn),y1(1:pn),'b',x2(1:pn),y2(1:pn),'g','linewidth',1.5)
    axis square

    line([0 x1(pn)], [0 y1(pn)],'Linewidth',1.5);
    axis([-(l1+l2) l1+l2 -(l1+l2) l1+l2]);
    line([x1(pn) x2(pn)], [y1(pn) y2(pn)],'linewidth',1.5);
    h=gca; 
    get(h,'fontSize');
    set(h,'fontSize',12)
    xlabel('X','fontSize',12);
    ylabel('Y','fontSize',12);
    title('Double Pendulum','fontsize',14)
    pause(0.03)

end

figure('Position',[99   150   698   596]);
   plot(x1,y1,'b','linewidth',1.5)
   hold on
   plot(x2,y2,'g','linewidth',1.5)
   legend('Position 1','Position2')
   xlabel('X Position','fontSize',14);
   ylabel('Y Position','fontSize',14);
   title('Position of Mass1 and Mass2','fontsize',14)
   %axis square
   grid on 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       


   case 3 %third plot
                              %%%%CASE 3%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure('Position',[40   154   649   563]);
   plot(t,u1,'r','linewidth',1.25)
   hold on 
   plot(t,u2,'b','linewidth',1.25)
   grid
   legend('\Theta_1','\Theta_2')
   xlabel('Time (sec)','fontSize',14);
   ylabel('Angle \Theta','fontSize',14);
   title('Angle: \Theta_1 and \Theta_2 vs Time','fontsize',14)

   figure('Position',[715   155   649   563]);
   plot(t,v1,'r','linewidth',1.25)
   hold on 
   plot(t,v2,'b','linewidth',1.25)
   grid
   legend('d\Theta_1/dt','d\Theta_2/dt')
   xlabel('Time (sec)','fontSize',14);
   ylabel('Angular Velocity (d\Theta/dt)','fontSize',14);
   title('Angular Velocity: (d\Theta_1/dt) and (d\Theta_2/dt) vs Time','fontsize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


case 4 %fourth plot
                             %%%%CASE 4%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
figure('Position',[343    33   768   760]);
for pn = 1:1:length(t)
    clf
    plot(x1(1:pn),y1(1:pn),'b','linewidth',1.5)
    hold on
    plot(x2(1:pn),y2(1:pn),'g-o','linewidth',1.5)
    grid
    axis([-(l1+l2) l1+l2 -(l1+l2) l1+l2]);
    h=gca; 
    get(h,'fontSize');
    set(h,'fontSize',12)
    xlabel('X Position','fontSize',12);
    ylabel('Y Position','fontSize',12);
    title('X and Y Position over Time','fontsize',14)
    pause(0.03)
end
    legend('Mass 1','Mass 2')
    
    case 5 %ENERGY IN SYSTEM PLOT 
figure('Position',[97   22   1216   776]); %make initial size large       
plot(t,PeE,'g','linewidth',1.5)
hold on
plot(t,KeE,'b','linewidth',1.5)
hold on 
plot(t,EtotalE,'r','linewidth',1.5)
title('System Energy: Mass 1')
xlabel('Time (sec)')
ylabel('Energy (Joules)')
set(gca,'fontsize',18)
h = legend('PE Euler','KE Euler', 'Total Energy Euler','location','northeast');
set(h,'FontSize',14);
grid
    
    case 6 %Poincare Map #1
        
        % Note: 
        %We need a different time vector so we need to redo the integration
        %for the two Poincre maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%new vars%%% 
l1=1; 
l2=2; 
m1=2; 
m2=1; 
g=9.8; 
%Initial Conditions
u10 = 1.6; 
u20 = 0; 
v10 = 2.2; 
v20 = 0; 
%time vector%
dt = 0.015; %Time Step
tStop = 130; %Record Length
t = 0:dt:tStop; %Time Vector

u1 = zeros(1, length(t)); %u1 = theta1(t)
u1(1) = u10; %set initial theta1
u2 = zeros(1, length(t)); %u2 = theta1Prime(t), or angular velocity 1
u2(1) = u20; %set initial velocity of first mass
v1 = zeros(1, length(t)); %v1 = theta2(t)
v1(1) = v10; %set initial theta2
v2 = zeros(1, length(t)); %v1 = theta2(t)
v2(1) = v20; %set initial theta2

%Euler Method
for i=2:length(t)
 
    a = (m1+m2)*l1 ;
    b = m2*l2*cos(u1(i-1)-v1(i-1)) ;
    c = m2*l1*cos(u1(i-1)-v1(i-1)) ;
    d = m2*l2 ;
    e = -m2*l2*v2(i-1)* v2(i-1)*sin( u1(i-1) - v1(i-1) )-g*(m1+m2)*sin(u1(i-1)) ;
    f = m2*l1*u2(i-1)*u2(i-1)*sin(u1(i-1) - v1(i-1))-m2*g*sin(v1(i-1)) ;
    
                           %%%%%%%%%%%%%%
    
    du2dt= (e*d-b*f)/(a*d-c*b) ;
    u2(i) = u2(i-1) + dt * (du2dt);
    dv2dt= (a*f-c*e)/(a*d-c*b) ;
    v2(i) = v2(i-1) + dt * (dv2dt);
    du1dt = (u2(i)); 
    u1(i) = u1(i-1) + dt * (du1dt);
    dv1dt = (v2(i)); 
    v1(i) = v1(i-1) + dt * (dv1dt);
   
end

% find the x and y positions 
x1 = l1*sin(u1);
y1 = -l1*cos(u1);
x2 = l1*sin(u1)+l2*sin(v1);
y2 = -l1*cos(u1)-l2*cos(v1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %Poincare Plot #1
        figure('Position',[97   22   1216   776]); %make initial size large
        U1 = u1; 
        plot(u2,u1,'black.','markersize',(6))
        xlabel('\Theta_1','fontSize',15)
        ylabel('d\Theta_1dt','fontSize',15)
        title('Poincare Map: d\Theta_1dt vs \Theta_1','fontSize',20)
        h = legend('Point plotted every 0.015 seconds');
        set(h,'FontSize',11);
        grid
        
    case 7 %Poincare Map and Chaotic Fractal Creation #2
        %%%declare some new variables for seconds map%%% 
        
l1=1; 
l2=1.4; 
m1=2.5; 
m2=1.6; 
g=9.8; 
%Initial Conditions
u10 = 1.6; 
u20 = 0; 
v10 = 2.2; 
v20 = 0; 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%time vector%
dt = 0.005; %Time Step
tStop = 250; %Record Length
t = 0:dt:tStop; %Time Vector

u1 = zeros(1, length(t)); %u1 = theta1(t)
u1(1) = u10; %set initial theta1
u2 = zeros(1, length(t)); %u2 = theta1Prime(t), or angular velocity 1
u2(1) = u20; %set initial velocity of first mass
v1 = zeros(1, length(t)); %v1 = theta2(t)
v1(1) = v10; %set initial theta2
v2 = zeros(1, length(t)); %v1 = theta2(t)
v2(1) = v20; %set initial theta2

%Euler Method
for i=2:length(t)
 
    a = (m1+m2)*l1 ;
    b = m2*l2*cos(u1(i-1)-v1(i-1)) ;
    c = m2*l1*cos(u1(i-1)-v1(i-1)) ;
    d = m2*l2 ;
    e = -m2*l2*v2(i-1)* v2(i-1)*sin( u1(i-1) - v1(i-1) )-g*(m1+m2)*sin(u1(i-1)) ;
    f = m2*l1*u2(i-1)*u2(i-1)*sin(u1(i-1) - v1(i-1))-m2*g*sin(v1(i-1)) ;
    
                           %%%%%%%%%%%%%%
    
    du2dt= (e*d-b*f)/(a*d-c*b) ;
    u2(i) = u2(i-1) + dt * (du2dt);
    dv2dt= (a*f-c*e)/(a*d-c*b) ;
    v2(i) = v2(i-1) + dt * (dv2dt);
    du1dt = (u2(i)); 
    u1(i) = u1(i-1) + dt * (du1dt);
    dv1dt = (v2(i)); 
    v1(i) = v1(i-1) + dt * (dv1dt);
   
end

% find the x and y positions 
x1 = l1*sin(u1);
y1 = -l1*cos(u1);
x2 = l1*sin(u1)+l2*sin(v1);
y2 = -l1*cos(u1)-l2*cos(v1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        %Poincare Plot #2
        figure('Position',[97   22   1216   776]); %make initial size large
        plot(u2,v2,'black.','markersize',(0.5)); axis equal; xlabel('position');ylabel('velocity');title('trajectory')
         xlabel('d\Theta_1dt','fontSize',15)
        ylabel('d\Theta_2dt','fontSize',15)
        title('Poincare Map: d\Theta_1dt vs d\Theta_2dt ','fontSize',20)
        h = legend('Point plotted every 0.005 seconds');
        set(h,'FontSize',11);
        grid minor
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---------------------------------------------------------------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %proof of concept to make sure the ode45 solution also works 
    case 8 %Solutions using ode45
        %%%PROOF OF CONCEPT FOR ODE45 SOLVER%%%
        figure('Position',[343    33   768   760]);
for pn = 1:1:length(t)
    clf
    plot(xx1(1:pn),yy1(1:pn),'b','linewidth',1.5)
    hold on
    plot(xx2(1:pn),yy2(1:pn),'g-o','linewidth',1.5)
    grid
    axis([-(l1+l2) l1+l2 -(l1+l2) l1+l2]);
    h=gca; 
    get(h,'fontSize');
    set(h,'fontSize',12)
    xlabel('X Position','fontSize',12);
    ylabel('Y Position','fontSize',12);
    title('4th Order Runga Kutta: X and Y Position','fontsize',14)
    pause(0.03)
end
    legend('Mass 1','Mass 2')
        
end %ends switch statement

%FUNCTION USED IN ODE45 SOLVER IN LINE 121
function [dydt] = solution(tt, yy)

% These conditions will be overwritten, but the differential equations below
% need to act on variables so we need to declare them in the function as
% well. 
l1=1;%placeholder
l2=1;%placeholder
m1=1;%placeholder
m2=1;%placeholder
g=9.8;%placeholder

aa = (m1+m2)*l1 ;
bb = m2*l2*cos(yy(1)-yy(3)) ;
cc = m2*l1*cos(yy(1)-yy(3)) ;
dd = m2*l2 ;
ee = -m2*l2*yy(4)* yy(4)*sin(yy(1)-yy(3))-g*(m1+m2)*sin(yy(1)) ;
ff = m2*l1*yy(2)*yy(2)*sin(yy(1)-yy(3))-m2*g*sin(yy(3)) ;

dydt(1) = yy(2);
dydt(3)= yy(4) ;
dydt(2)= (ee*dd-bb*ff)/(aa*dd-cc*bb) ;
dydt(4)= (aa*ff-cc*ee)/(aa*dd-cc*bb) ;
dydt=dydt';

end

%References: (1) http://scienceworld.wolfram.com/physics/DoublePendulum.html
% (2) https://www.mathworks.com/help/matlab/matlab_prog/create-functions-in-files.html
% (3) https://www.myphysicslab.com/pendulum/double-pendulum/double-pendulum-en.html
% (4) http://www.physicsandbox.com/projects/double-pendulum.html

