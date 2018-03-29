%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL EXAM: PY331
%~~~~~~~~~~~~~~~~~~~
% Program #3: Fields--- Particles interacting with E and B fields
%
%
%
% Author: Spencer R Bertsch
% Date: May 9, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all variables
clear
clc
%%%%%%%%%%%%%%%%%%%%



%% Analytical Solution
    
%Negatively charged particle motion in perpenicular magentic field
%and electric fields

dt = 0.01;  
tA = 0:dt:30;
a = 2; %arbitrary constant 
b = 2;  

%Analytic Solution for Cycloid Trajectory
%st = sin(t); %Helix
%ct = cos(t); %Helix
yAnalytical = (a*tA - b*sin(tA));
zAnalytical = (a - b*cos(tA));

%% Euler Method

% time vector 
dt = 0.01;  
t = 0:dt:70;

%Numerical Solution for Cycloid Trajectory
E = 2; %electric field strength
B = 4; %magnetic field strength
omega = 0.5; 

Ey0 = 0; 
Ez0 = 0; 

Ey = zeros(1,length(t)); % Initializes the delta x vector
Ey(1) = Ey0; % Set initial delta x to 0

Ez= zeros(1,length(t)); % Initializes the delta x vector
Ez(1) = Ez0; % Set initial delta x to 0

for i=2:length(t)

    dydt = 4*((E/B) * (1- cos(omega*t(i)))); 
    Ey(i) = Ey(i-1) + dt * (dydt);
    dzdt = 4*(E/B) * sin(omega*(t(i))); 
    Ez(i) = Ez(i-1) + dt * (dzdt);

end

%% Euler Chromer Method
 
ECy0 = 0; 
ECz0 = 0; 

ECy = zeros(1,length(t)); % Initializes the delta x vector
ECy(1) = ECy0; % Set initial delta x to 0

ECz= zeros(1,length(t)); % Initializes the delta x vector
ECz(1) = ECz0; % Set initial delta x to 0

for i=2:length(t)

    dydt = 4*((E/B) * (1- cos(omega*t(i)))); 
    ECy(i) = ECy(i-1) + dt * (dydt);
    dzdt = 4*(E/B) * sin(omega*(t(i-1))); 
    ECz(i) = ECz(i-1) + dt * (dzdt);

end

%% User Input
disp('Enter "1" to see an electric dipole.')
disp('Enter "2" to see the analytic solution.')
disp('Enter "3" to see the numerical solution using the Euler Method')
disp('Enter "4" to see the numerical solution using the Euler Chromer Method')
disp('Enter "5" to see cycloid motion of the particle')
problemnumber = input('Enter a number: ');
switch problemnumber 
%% Section One 
    case 1 %represents the first case in the switch statement

        %%% Positive and Negative charge in 2D Plane %%% 
        
[X,Y] = meshgrid(-2:0.2:2);
Z = X .* exp(-X.^2 - Y.^2);
figure('Position',[204    52   907   745]); %make initial size of figure large

[C,H] = contour (X, Y, Z, 25);
set (H, 'LineWidth', 1.25); %make lines thicker 

[U,V] = gradient(Z,0.2,0.2);
hold on
 lh=quiver(X,Y,U,V);
 set(lh,'linewidth',0.7);
 set(lh,'color',[0,0,1]);

title('Electric Dipole');
set(findall(gca, 'Type', 'Line'),'LineWidth',12);
set(gca,'fontsize',20)
%set(gca,'color','black')
hold off

%% PLOTS 
case 2

%%% ***** Plot of Analytical Solution *** %%%

figure('Position',[204    52   907   745]); %make initial size large
for k = 1:60:length(tA)
    clf
    plot3(yAnalytical(1:k),zAnalytical(1:k),tA(1:k),'b','linewidth',1.5)
    hold on 
    plot3(yAnalytical(k),zAnalytical(k),tA(k),'c.','markersize',(20)); %<--- Uncomment to see
%     electron moving through path
    title([num2str(tA(k)), ' Milliseconds']);
    %Uncomment next line to hold axis static 1
    %axis([ min([st(:);ct(:);t(:)]) max([st(:);ct(:);t(:)]) min([st(:);ct(:);t(:)]) max([st(:);ct(:);t(:)]) min([st(:);ct(:);t(:)]) max([st(:);ct(:);t(:)]) ])
    axis([ 0 61.9761 0 5.9761 0 41.9761 ])
    xlabel('Time (ms)')
    ylabel('Y Position')
    zlabel('X Position')
    set(gca,'fontsize',20)
    grid on 
    hold on 
    
    drawnow 

end
disp('Analytical Solution: 3D motion of the electron traveling through') 
disp('perpendicular electric and magnetic fields')

    case 3 
%%% ***** Plot of Euler Method *** %%%

figure('Position',[204    52   907   745]); %3D image
for k = 1:100:length(t)
    clf
    plot3(Ey(1:k),Ez(1:k),t(1:k),'r','linewidth',1.5)
    hold on %Uncomment to see electron moving through path 
    plot3(Ey(k),Ez(k),t(k),'c.','markersize',(20));
    axis([ 0 161.9761 0 10.9761 0 91.9761 ])
    title([num2str(t(k)), ' Milliseconds']);
    xlabel('Time (ms)')
    ylabel('Y Position')
    zlabel('X Position')
    set(gca,'fontsize',20)
    grid

    drawnow 
   
end
disp('Numerical Solution using Euler Method: 3D motion of the electron ')
disp('traveling through perpendicular electric and magnetic fields')

    case 4
        %%% ***** Plot of Euler Chromer Method *** %%%
        
figure('Position',[204    52   907   745]); %3D image
for k = 1:100:length(t)
    clf
    plot3(ECy(1:k),ECz(1:k),t(1:k),'r','linewidth',1.5)
    hold on %Uncomment to see electron moving through path 
    plot3(ECy(k),ECz(k),t(k),'c.','markersize',(20));
    axis([ 0 161.9761 0 10.9761 0 91.9761 ])
    title([num2str(t(k)), ' Milliseconds']);
    xlabel('Time (ms)')
    ylabel('Y Position')
    zlabel('X Position')
    set(gca,'fontsize',20)
    grid

    drawnow 
   
end
disp('Numerical Solution using Euler Chromer Method: 3D motion of the electron ')
disp('traveling through perpendicular electric and magnetic fields')
        
    case 5
        %%% ***** Plot of Cycloid Numerical Solution in 2D *** %%%

%figure('Position',[204    52   907   745]); %3D image
figure('Position',[59  311  1297  277]); %2D cycloid
for k = 1:100:length(t)
    clf
    plot3(Ey(1:k),Ez(1:k),t(1:k),'r','linewidth',1.2)
    hold on 
    plot(Ey(k),Ez(k),'c.','markersize',(20));
%     hold on 
%     plot3(y(1:k),z(1:k),t(1:k),'g--','linewidth',2)
    title([num2str(t(k)), ' Milliseconds']);
    %Uncomment next line to hold axis static 1
    xlabel('Time (ms)')
    ylabel('Y Position')
    zlabel('X Position')
    set(gca,'fontsize',20)

    %axis([min([z(:);y(:)]) max([z(:);y(:)]) min([z(:);y(:)]) max([z(:);y(:)])])
    axis([-1.2384e-05 141.7318 -1.2384e-05 10.7318])
        %Camera Stuff
    campos([75 3.5000 916.0254])
    camtarget([75 3.5000 50])
    camup([0 1 0])
    
    %creates cycloid plot%%%%%%%%%%  
    pbaspect([1 0.1209 0.1209])
    %daspect([15.0000 1 10.0000])
    
    drawnow 
end
h = legend('Cycloid motion of the electron','Negatively charged electron');
set(h,'FontSize',12);
disp('Numerical Solution: Cycloid motion of the electron traveling through') 
disp('perpendicular electric and magnetic fields')

end % ends switch statement


% WORK CITED
% http://jwilson.coe.uga.edu/EMAT6680Fa2014/Gieseking/Exploration%2010/Parametric%20Equations.html
% http://physicstasks.eu/402/the-motion-of-a-charged-particle-in-homogeneous-perpendicular-electric-and-magnetic-fields
