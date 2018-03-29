%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL EXAM: PY331
%~~~~~~~~~~~~~~~~~~~
% Program #4: Van der Pol Oscillator
%
%
% Author: Spencer Bertsch
% Date 9, 2017 at 6:00pm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is another system of differential equations I found and I wanted to
% test out different methods of integration and compare their results. 
% It's a brief example, but it shows the similarities and differences
% between the schemes. 

clc
clear

dt = 0.1; 
t = 0:dt:18; 

y10 = 2; 
y20 = 0; 

prev1_0 = 2; 
prev2_0 = -0.2; 

Ey1 = zeros(1, length(t));     % Initialize the velocity vector
Ey1(1) = y10;      % Set initial velocity
Ey2 = zeros(1,length(t)); % Initializes the delta x vector
Ey2(1) = y20; % Set initial delta x to 0

ECy1 = zeros(1, length(t));     % Initialize the velocity vector
ECy1(1) = y10;      % Set initial velocity
ECy2 = zeros(1,length(t)); % Initializes the delta x vector
ECy2(1) = y20; % Set initial delta x to 0

Hy1 = zeros(1, length(t));     % Initialize the velocity vector
Hy1(1) = y10;      % Set initial velocity
Hy2 = zeros(1,length(t)); % Initializes the delta x vector
Hy2(1) = y20; % Set initial delta x to 0

prev1 = zeros(1,length(t)); % Initialize previous value vector
prev2 = zeros(1,length(t)); % Initialize previous value vector

prev2(1) = prev2_0; % Set initial value based on euler solution
prev1(1) = prev1_0; % Set initial value based on euler solution

new_position1 = zeros(1, length(t));     % Initialize the velocity vector
new_position1(1) = 1.98;      % Set initial velocity
new_position2 = zeros(1,length(t)); % Initializes the delta x vector
new_position2(1) = -2; % Set initial delta x to 0

% Heun Method -it took a very long time to figure out how to apply the heun
% method to a system of two coupled ODEs 
for i=2:length(t)
    
    dy1dt = (Hy2(i-1));%
    
    new_position1(i) = Hy1(i-1) + dt * (dy1dt); 
    
    new_slope1 = new_position2(i-1); %new slope that will be average with the old slope
    
    slope_now1 = ((dy1dt + new_slope1)/2); 
    
    Hy1(i) = Hy1(i-1) + dt * slope_now1; %new position for y1
    
    
    dy2dt = (1-Hy1(i-1)^2)*Hy2(i-1)-Hy1(i-1);%
    
    new_position2(i) = Hy2(i-1) + (dt * dy2dt); 
    
    new_slope2 = (1 - new_position1(i)^2) * new_position2(i) - new_position1(i); %new slope that will be average with the old slope
    
    slope_now2 = ((dy2dt + new_slope2)/2); 
    
    Hy2(i) = Hy2(i-1) + dt * slope_now2; %new position for y2

end

% Verlet Method 
% No dvdt, so no acceleration. 

% Euler Method
for i=2:length(t)

    dy2dt = (1-Ey1(i-1)^2)*Ey2(i-1)-Ey1(i-1);
    Ey2(i) = Ey2(i-1) + dt * (dy2dt);
    
    dy1dt = (Ey2(i));
    Ey1(i) = Ey1(i-1) + dt * (dy1dt);
    
end

%Euler Chromer Method
for i=2:length(t)

    dy2dt = (1-ECy1(i-1)^2)*ECy2(i-1)-ECy1(i-1);
    ECy2(i) = ECy2(i-1) + dt * (dy2dt);
    
    dy1dt = (ECy2(i-1));
    ECy1(i) = ECy1(i-1) + dt * (dy1dt);

end

%ODE45 Solver / ode23 also works well here 
t_start = 0; 
t_stop = 15; 
[tt,y] = ode45(@function1,[t_start t_stop],[y10; y20]);
Y1 = y(:,1); %first vector represents y1 
Y2 = y(:,2); %second vector represents y2

%%% Error Plot 
EError = zeros(length(t),1); 
ECError = zeros(length(t),1);
HError = zeros(length(t),1);
for i=1:length(t)
    EError(i) = ((Y1(i) - Ey1(i))/Y1(i));
    ECError(i) = ((Y1(i) - ECy1(i))/Y1(i));
    HError(i) = ((Y1(i) - Hy1(i))/Y1(i));
end 

EError1 = EError(1:50); 
ECError1 = ECError(1:50);
HError1 = HError(1:50);

EError2 = EError(1:120); 
ECError2 = ECError(1:120);
HError2 = HError(1:120);

disp('Enter "1" to see the Euler solution.')
disp('Enter "2" to see the Euler Chromer solution.')
disp('Enter "3" to see the Heun solution.')
disp('Enter "4" to see the 4th order Runga Kutta solution.')
disp('Enter "5" to see the Error against the 4th order R.C solution over 50 time steps.')
disp('Enter "6" to see the Error against the 4th order R.C solution over 120 time steps.')

problemnumber = input('Enter a number: ');

switch problemnumber 
    case 1 %Euler
figure('Position',[204    52   907   745]);
plot(t,Ey1,'g-o','LineWidth',0.8)
hold on 
plot(t,Ey2,'b-o','LineWidth',0.8)
title('Van der Pol Equation- Euler Solutions for y1 and y2','fontSize',15);
xlabel('Time (sec)','fontSize',15);
ylabel('Solution','fontSize',15);
h = legend('y1','y2');
set(h,'FontSize',13);
grid on 

case 2 %Euler Chromer
figure('Position',[204    52   907   745]);
plot(t,ECy1,'g-o','LineWidth',0.8)
hold on 
plot(t,ECy2,'b-o','LineWidth',0.8)
title('Van der Pol Equation- Euler Chromer Solutions for y1 and y2','fontSize',15);
xlabel('Time (sec)','fontSize',15);
ylabel('Solution','fontSize',15);
h = legend('y1','y2');
set(h,'FontSize',13);
grid on  

case 3 %Heun
figure('Position',[204    52   907   745]);
plot(t,Hy1,'g-o','LineWidth',0.8)
hold on 
plot(t,Hy2,'b-o','LineWidth',0.8)
title('Van der Pol Equation- Heun Solutions for y1 and y2','fontSize',15);
xlabel('Time (sec)','fontSize',15);
ylabel('Solution','fontSize',15);
h = legend('y1','y2');
set(h,'FontSize',13);
grid on 

case 4 %ode45- 4th order Runga Kutta
figure('Position',[204    52   907   745]);
plot(tt,Y1,'g-o','LineWidth',0.8)
hold on 
plot(tt,Y2,'b-o','LineWidth',0.8)
title('Van der Pol Equation- 4th order Runga Kutta Solutions for y1 and y2','fontSize',15);
xlabel('Time (sec)','fontSize',15);
ylabel('Solution','fontSize',15);
h = legend('y1','y2');
set(h,'FontSize',13);
grid on 

case 5 %error compared to 4th order runga kutta
t = t(1:50); 
figure('Position',[204    52   907   745]);
plot(t,EError1,'*','MarkerSize',14)
hold on 
plot(t,ECError1,'^','MarkerSize',14)
hold on 
plot(t,HError1,'o','MarkerSize',14)
title('Error of Integration methods against 4th order Runga Kutta Solution','fontSize',15);
xlabel('Time (sec)','fontSize',15);
ylabel('Solution','fontSize',15);
h = legend('Euler Error','Euler Chromer Error','Heun Error');
set(h,'FontSize',13,'Location','northwest');
grid on

case 6 %error compared to 4th order runga kutta--- longer time interval 
t = t(1:120); 
figure('Position',[204    52   907   745]);
plot(t,EError2,'*','MarkerSize',14)
hold on 
plot(t,ECError2,'^','MarkerSize',14)
hold on 
plot(t,HError2,'o','MarkerSize',14)
title('Error of Integration methods against 4th order Runga Kutta Solution','fontSize',15);
xlabel('Time (sec)','fontSize',15);
ylabel('Solution','fontSize',15);
h = legend('Euler Error','Euler Chromer Error','Heun Error');
set(h,'FontSize',13,'Location','northwest');
title(h,'Error over Time')
grid on
    
end %ends switch statement 
%functions
function output = function1(~,y)
output = [y(2); (1-y(1)^2)*y(2)-y(1)];
end

%references: 
%(1) https://en.wikipedia.org/wiki/Van_der_Pol_oscillator
%(2) http://mathworld.wolfram.com/vanderPolEquation.html
