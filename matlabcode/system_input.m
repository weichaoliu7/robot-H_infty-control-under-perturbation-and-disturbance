% model of two link robot manipulator
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
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [0 0];
function sys=mdlOutputs(t,x,u)
qd1=1.5*cos(t)+0.5; % desired angular displacement of link 1
dqd1=-1.5*sin(t); % desired angular velocity of link 1
ddqd1=-1.5*cos(t); % desired angular acceleration of link 1
qd2=cos(t)+1; % desired angular displacement of link 2
dqd2=-sin(t); % desired angular velocity of link 2
ddqd2=-cos(t); % desired angular acceleration of link 2

sys(1)=qd1;
sys(2)=dqd1;
sys(3)=ddqd1;
sys(4)=qd2;
sys(5)=dqd2;
sys(6)=ddqd2;