function [sys,x0,str,ts] = spacemodel(t,x,u,flag)
switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1,
    sys=mdlDerivatives(t,x,u);
case 3,
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 1;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 4;
sizes.NumInputs      = 10;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = 0;
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)
sys=0;

function sys=mdlOutputs(t,x,u)
qd1=u(1);
dqd1=u(2);
ddqd1=u(3);
qd2=u(4);
dqd2=u(5);
ddqd2=u(6);
dqd=[dqd1 dqd2]';
ddqd=[ddqd1 ddqd2]';

q1=u(7); % actual angular displacement of link 1
q2=u(9); % actual angular displacement of link 2
q=[q1 q2]';
dq1=u(8); % actual angular velocity of link 1
dq2=u(10); % actual angular velocity of link 2
dq=[dq1 dq2]';

e1=q1-qd1; % angular displacement tracking error of link 1
e2=q2-qd2; % angular displacement tracking error of link 2
e=[e1 e2]'; % second term of state tracking error
de1=dq1-dqd1; % angular velocity tracking error of link 1
de2=dq2-dqd2; % angular velocity tracking error of link 2
de=[de1 de2]'; % first term of state tracking error
error_state = [de;e]; % state tracking error, defined in Eq. 7

l1=1; % length of link 1
l2=1; % length of link 2
m1=1; % mass of link 1
m2=10; % mass of link 2

% inertia matrix of manipulator nominal model
M0(1,1)=(m1+m2)*l1^2;
M0(1,2)=m2*l1*l2*(sin(u(7))*sin(u(9))+cos(u(7))*cos(u(9)));
M0(2,1)=M0(1,2);
M0(2,2)=m2*l2^2;

% Coriolis/centrifugal force matrix of manipulator nominal model
C0(1,1)=0;
C0(1,2)= - m2*l1*l2*(cos(u(7))*sin(u(9))-sin(u(7))*cos(u(9)))*u(10);
C0(2,1)= - m2*l1*l2*(cos(u(7))*sin(u(9))-sin(u(7))*cos(u(9)))*u(8);
C0(2,2)=0;

% gravitational matrix of manipulator nominal model
g = 9.794; % gravitational acceleration
G0(1)=-(m1+m2)*l1*g*sin(u(7));
G0(2)=-m2*l2*g*sin(u(9));
G0 = G0';

% robotic H_infty control design
% choose desired disturbance attenuation level gamma and positive definite weight matrix R
scenario = 'case4';
switch scenario
    case 'case1'
        gamma = inf;
        beta = 0.09;
    case 'case2'
        gamma = 1;
        beta = 0.04;
    case 'case3'
        gamma = 0.2;
        beta = 0.01;
    case 'case4'
        gamma = 0.05;
        beta = 0.002;
    otherwise
        disp('gamma does not have the value');
end

% calculate the Cholesky factorization
R = beta * eye(2);
R1 = 1/sqrt(1/beta - 1/(gamma^2)) * eye(2);

% select positive definite symmetric weight matrix Q, defined in Eq. 47
Q1 = eye(2);
Q2 = eye(2);
Q12 = zeros(2, 2);
Q21 = Q12;
Q = [Q1'*Q1 Q12; Q21 Q2'*Q2];

% solution matrix of Riccati-like algebraic equation (43), defined in Eq. 10
T11 = R1' * Q1;
T12 = R1' * Q2;
T21 = zeros(2);
T22 = eye(2);
T0 = [T11, T12;T21, T22];

B1 = eye(2);
B2 = zeros(2);
B = [B1;B2]; % input matrix in the state space method, defined in Eq. 8

% corresponding optimal applied torque with nominal parameter matrix,defined in Eq. 57
control_optimal = - inv(R) * B' * T0 * error_state; % state feedback control law design, defined in Eq. 58
torque_optimal = M0 * (ddqd - inv(T11) * T12 * de  - inv(T11) * inv(M0) * (C0 * B' * T0 * error_state - control_optimal) ) + C0 * dq + G0;

sys(1)=control_optimal(1);
sys(2)=control_optimal(2);
sys(3)=torque_optimal(1);
sys(4)=torque_optimal(2);