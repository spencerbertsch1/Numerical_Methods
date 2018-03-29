%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL EXAM: PY331
%~~~~~~~~~~~~~~~~~~~
% Program #1: Lorentz Attractor
%
%
%
% Author: Spencer Bertsch
% Date May 9, 2017 at 6:00pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all variables
clear
clc
%%%%%%%%%%%%%%%%%%%%

%% User Input
disp('Enter "1" to see X, Y, and Z components of the lorentz attractor plotted seperately against time.')
disp('Enter "2" to see a three dimensional representation of the attractor in real time.')
disp('Enter "3" to see a three dimensional representation with a different color scheme.')

problemnumber = input('Enter a number: ');

%% Variables
% Declare Variables 
sigma=10;
beta=8/3;
rho=28;

%Initial Conditions
x0 = 1; 
y0 = 1; 
z0 = 1; 

%time vector
dt = 0.005; %sec (appropriate dt = 0.005 sec)
t = 0:dt:70; %seconds 

%% Numerical Solution Euler Chromer Method 

x = zeros(1, length(t));     % Initialize the X vector
x(1) = x0;      % Set initial x

y = zeros(1,length(t)); % Initializes the Y vector
y(1) = y0; % Set initial y

z= zeros(1,length(t)); % Initializes the Z vector
z(1) = z0; % Set initial z

%Euler Method solution to system of differential equations
for i=2:length(t) 
    dxdt = sigma * (y(i-1) - x(i-1)); 
    x(i) = x(i-1) + dt * (dxdt); %Equation 1
    dydt = x(i-1) * (rho - z(i-1)) -y(i-1); 
    y(i) = y(i-1) + dt * (dydt); %Equation 2
    dzdt = (x(i-1) * y(i-1)) -beta*z(i-1); 
    z(i) = z(i-1) + dt * (dzdt); %Equation 3
end

%% Plots

switch problemnumber  
    case 1 %represents the first case in the switch statement
%% 2D PLOTTING %%% 
figure('Position',[204    52   907   745]);
subplot(3,1,1)
plot(t,x,'r','linewidth',1)
title('X Component of Attractor')
xlabel('Time (sec)')
ylabel('X Position')
grid

subplot(3,1,2)
plot(t,y,'g','linewidth',1)
title('Y Component of Attractor')
xlabel('Time (sec)')
ylabel('Y Position')
grid

subplot(3,1,3)
plot(t,z,'b','linewidth',1)
title('Z Component of Attractor')
xlabel('Time (sec)')
ylabel('Z Position')
grid

    case 2 %second plot
%%% 3D Plotting %%%
figure('Position',[204    52   907   745]); %make initial size large
for k = 1:35:length(t)
    plot3(x(1:k),y(1:k),z(1:k),'b','linewidth',0.8)
    title([num2str(t(k)), ' Seconds']);
    axis([ -24.9823  30.4832  -24.9823  30.4832  -10.9823  60.4832 ])
    ylabel('Y Position')
    zlabel('Z Position')
    set(gca,'fontsize',20)
    grid on
    
    %%% rotating view angle #1
     campos([(-155.5363 + k/20) (-389.7880 - k/20) (314.5237 + k/60)]) 
     camtarget([5 0 30])
      
    drawnow 

end

    case 3 %third plot
        %%% different color scheme just for fun %%% 
figure('Position',[204    52   907   745]); %make initial size large
for k = 1:35:length(t)
    plot3(x(1:k),y(1:k),z(1:k),'c','linewidth',0.95)
    title([num2str(t(k)), ' Seconds']);
    axis([ -24.9823  30.4832  -24.9823  30.4832  -10.9823  60.4832 ])
    xlabel('X Position')
    ylabel('Y Position')
    zlabel('Z Position')
    set(gca,'fontsize',20)
    colormap winter
    grid on 
    set(gca,'color','black') %<---vvv UNCOMMENT to make background black
    ax = gca;
    ax.GridColor = 'white';  % [R, G, B]
    
    %%% rotating view angle #1
     campos([(-155.5363 + k/50) (-389.7880 + k/8) (314.5237 - k/80)])
     camtarget([5 0 30])
    
    drawnow 

end

end %ends switch statement

%Work Cited 
% https://en.wikipedia.org/wiki/Lorenz_system