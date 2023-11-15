function [sys,x0,str,ts]=plant(t,x,u,flag)
switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1,
    sys=mdlDerivatives(t,x,u);
case 3,
    sys=mdlOutputs(t,x,u);
case {2, 4, 9 }
    sys = [];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 4;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 4;
sizes.NumInputs      = 4;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;
sys=simsizes(sizes);
x0=[-2 0 -2 0];
str=[];
ts=[];
function sys=mdlDerivatives(t,x,u) 
l1=1; % length of link 1
l2=1; % length of link 2
m1=1; % mass of link 1
m2=10; % mass of link 2

dq = [x(2);x(4)]; % actual angular velocity
torque=[u(3) u(4)]'; % optimal applied torque

% inertia matrix of manipulator nominal model
M0(1,1)=(m1+m2)*l1^2;
M0(1,2)=m2*l1*l2*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)));
M0(2,1)=M0(1,2);
M0(2,2)=m2*l2^2;

% uncertainty of M due to change of load
delta_M(1,1) = -5;
delta_M(1,2) = -5 * (sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3)));
delta_M(2,1)=delta_M(1,2);
delta_M(2,2) = -5;

% inertia matrix of manipulator practical model
M = M0 + delta_M;

% Coriolis/centrifugal force matrix of manipulator nominal model
C0(1,1)=0;
C0(1,2)= - m2*l1*l2*(cos(x(1))-sin(x(3))-sin(x(1))*cos(x(3)))*x(4);
C0(2,1)= - m2*l1*l2*(cos(x(1))-sin(x(3))-sin(x(1))*cos(x(3)))*x(2);
C0(2,2)=0;

% uncertainty of C due to change of inertia delta_M and friction
delta_C(1,1) = 0;
delta_C(1,2) = ( cos(x(1))*sin(x(2)) - sin(x(1))*cos(x(2)) ) * 5 * x(2)^2;
delta_C(2,1) = ( cos(x(1))*sin(x(2)) - sin(x(1))*cos(x(2)) ) * 5 * x(4)^2;
delta_C(2,2) = 0;

% Coriolis/centrifugal force matrix of manipulator practical model
C = C0 + delta_C;

% gravitational matrix of manipulator nominal model
g = 9.794;
G0(1)=-(m1+m2)*l1*g*sin(x(1));
G0(2)=-m2*l2*g*sin(x(3));
G0 = G0';

% perturbation of gravitational force due to configuration and change of manipulator total mass
delta_G(1)=5*g*sin(x(1));
delta_G(2)=5*g*sin(x(3));
delta_G = delta_G';

% gravitational matrix of manipulator practical model
G = G0 + delta_G;

% external perturbation is square wave with amplitude 2 and frequency 0.27 or 0.2
f1_disturbance = 0.27;
external_disturbance1 = 2*sign(sin(2*pi*f1_disturbance*t));
f2_disturbance = 0.2;
external_disturbance2 = 2*sign(sin(2*pi*f2_disturbance*t));
external_disturbance = [external_disturbance1;external_disturbance2];

S=inv(M)*(torque + external_disturbance - C*dq -G);

sys(1)=x(2); % dq1
sys(2)=S(1); % ddq1
sys(3)=x(4); % dq2
sys(4)=S(2); % ddq2
function sys=mdlOutputs(t,x,u)
sys(1)=x(1); % q1
sys(2)=x(2); % dq1
sys(3)=x(3); % q2
sys(4)=x(4); % dq2