%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Method of Characteristics, Minimum Length Nozzle            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DESCRIPTION: Script generates a minimum-length nozzle contour based on the method of
%  characeteristics. 
%
%  CREDITS: MOC logic and details are outlined in Zucrow and Hoffman (1967) and
%  Anderson (2003), if you want more details. 
%
%  AUTHOR: Ethan Labianca-Campbell, MSAA @ Purdue University
%  https://www.linkedin.com/in/ethan-labianca-campbell/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gamma: specific heat ratio (Cp/Cv of gas)
% Me: Design-required exit Mach number of the nozzle
% n: number of "characteristics" 


clear all
clc
tic
%% Gas properties %%
%gamma=input('Enter gamma : ');
gamma = 1.2;
%% Design parameters %%
%Me=input('Enter exit Mach no. : ');
Me = 6.5;
theta_max=PMFan(Me,gamma)*180/(2*pi);
D=input('Enter your throat diameter (in any unit)'); % Non-Dimensional y co-ordinate of throat wall ()
%% Incident expansion wave conditions %%
n=input('Enter number of characteristics lines (greater than 2) emanating from sharp corner throat : ');
theta_0=theta_max/n ; % splits to n angles
%% Characteristic parameter solver
[v,Cp,Cm,theta]=MOC_2D(theta_max,theta_0,n);


%% COMPUTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mach number and Mach angle at each node
node=0.5*n*(4+n-1);
M=zeros(1,node);
mu=zeros(1,node);
for i=1:node
    M(i)=flowprandtlmeyer(gamma,v(i),'nu');
    mu(i)=Mu(M(i));
end
%% Grid generator
figure(1)
i=1;
x=zeros(1,node);
y=zeros(1,node);
wall=theta_max;

%% Plot characteristic lines within nozzle
while (i<=n+1)
    if i==1
        x(i)=-D/(tand(theta(i)-mu(i)));
        y(i)=0;
        plot([0 x(i)],[D 0]);
        hold on
    elseif i==n+1
            x(i)=(y(i-1)-D-x(i-1)*tand((theta(i-1)+theta(i)+mu(i-1)+mu(i))*0.5))/(tand(0.5*(wall+theta(i)))-tand((theta(i-1)+theta(i)+mu(i-1)+mu(i))*0.5));
            y(i)=D+x(i)*tand(0.5*(wall+theta(i)));
            plot([x(i-1) x(i)],[y(i-1) y(i)]);
            hold on
            plot([0 x(i)],[D y(i)]);
            hold on
    else
            x(i)=(D-y(i-1)+x(i-1)*tand(0.5*(mu(i-1)+theta(i-1)+mu(i)+theta(i))))/(tand(0.5*(mu(i-1)+theta(i-1)+mu(i)+theta(i)))-tand(theta(i)-mu(i)));
            y(i)=tand(theta(i)-mu(i))*x(i)+D;
            plot([x(i-1) x(i)],[y(i-1) y(i)]);
            hold on
            plot([0 x(i)],[D y(i)]);
            hold  on
    end
    i=i+1;
    hold on
end
h=i;
k=0;
i=h;
for j=1:n-1
    while (i<=h+n-k-1)
        if (i==h)
            x(i)=x(i-n+k)-y(i-n+k)/(tand(0.5*(theta(i-n+k)+theta(i)-mu(i-n+k)-mu(i))));
            y(i)=0;
            plot([x(i-n+k) x(i)],[y(i-n+k) y(i)]);
            hold on
        else if (i==h+n-k-1)
                x(i)=(x(i-n+k)*tand(0.5*(theta(i-n+k)+theta(i)))-y(i-n+k)+y(i-1)-x(i-1)*tand((theta(i-1)+theta(i)+mu(i-1)+mu(i))*0.5))/(tand(0.5*(theta(i-n+k)+theta(i)))-tand((theta(i-1)+theta(i)+mu(i-1)+mu(i))*0.5));
                y(i)=y(i-n+k)+(x(i)-x(i-n+k))*tand(0.5*(theta(i-n+k)+theta(i)));
                plot([x(i-1) x(i)],[y(i-1) y(i)]);
                hold on
                plot([x(i-n+k) x(i)],[y(i-n+k) y(i)]);
                hold on
            else
                s1= tand(0.5*(theta(i)+theta(i-1)+mu(i)+mu(i-1)));
                s2= tand(0.5*(theta(i)+theta(i-n+k)-mu(i)-mu(i-n+k)));
                x(i)=(y(i-n+k)-y(i-1)+s1*x(i-1)-s2*x(i-n+k))/(s1-s2);
                y(i)=y(i-1)+(x(i)-x(i-1))*s1;
                plot([x(i-1) x(i)],[y(i-1) y(i)]);
                hold on
                plot([x(i-n+k) x(i)],[y(i-n+k) y(i)]);
                hold on
            end
        end
        i=i+1;
    end
    k=k+1;
    h=i;
    i=h;
    hold on
end
title(sprintf('Characteristic lines for Mach=%d and gamma=%d',Me,gamma))
xlabel('x/x0');
ylabel('y/y0');
axis equal
xlim([0 x(node)+0.5])
ylim([0 y(node)+0.5])
%% Nozzle co-ordinates (x_wall,y_wall)
x_wall=zeros(1,n+1);
y_wall=zeros(1,n+1);
x_wall(1)=0;
y_wall(1)=1;
i=2;
j=n+1;
while (i<=n+1)
    x_wall(i)=x(j); % x_wall outputs 
    y_wall(i)=y(j); % y_wall outputs
    j=j+(n-i+2);
    i=i+1;
end

%plot contour of inner nozzle shape
figure(2)
plot(x_wall,y_wall)

% output information to command
length = x_wall(end);
R_e = y_wall(end);
A_e = pi*R_e^2;
toc
fprintf('\n**OUTPUTS**\n')
fprintf('Length from throat is = %.2f \n',length)
fprintf('Radius of exit from is = %.2f \n',R_e)
fprintf('Area of exit is = %.2f \n',A_e)



%% USER DEFINED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRANDTL-MEYER FAN
function v=PMFan(M,gamma)
    % Prandtl-Meyer expansion fan angle; generates the angle by which sonic
    % flow must turn to reach given mach number. 
    v=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(M^2-1)))-atan(sqrt(M^2-1));
end

%% MACH ANGLE
function mu=Mu(M)
    % computes mach wave angle
    mu=asind(1/M);
end

%% 2D MOC LOGIC 
function [v,Cp,Cm,theta]=MOC_2D(theta_max,theta_0,n)
% v - Prandtl-Meyer function
% Cp - Left running characteristic constant (C+)
% Cm - Right running characteristic constant (C-)
% theta - Flow angle relative to horizontal

dtheta=(theta_max-theta_0)/(n-1); % angle step
%7+1,6+1,5+1,4+1,3+1,2+1,1+1
node=0.5*n*(4+n-1); % num nodes
theta=zeros(1,node); 
v=zeros(1,node);
Cp=zeros(1,node);
Cm=zeros(1,node);

for i=1:n
    theta(i)=theta_0+(i-1)*dtheta;
    v(i)=theta(i);
    Cp(i)=theta(i)-v(i);
    Cm(i)=theta(i)+v(i);
end

i=n+1;
theta(i)=theta(i-1);
v(i)=v(i-1);
Cp(i)=Cp(i-1);
Cm(i)=Cm(i-1);
p=2;
q=n+2;
for k=1:n-1
    j=p;
    h=q;
    theta(h)=0;
    Cm(h)=Cm(j);
    v(h)=Cm(j)-theta(h);
    Cp(h)=theta(h)-v(h);
    j=j+1;
    for i=h+1:n-p+q
        Cm(i)=Cm(j);
        Cp(i)=Cp(i-1);
        theta(i)=0.5*(Cp(i)+Cm(i));
        v(i)=0.5*(Cm(i)-Cp(i));
        j=j+1;
    end
    if i==n-p+q
        h=i+1;
    else
        h=h+1;
    end
    theta(h)=theta(h-1);
    v(h)=v(h-1);
    Cp(h)=Cp(h-1);
    Cm(h)=Cm(h-1);
    p=p+1;
    q=h+1;
end
end