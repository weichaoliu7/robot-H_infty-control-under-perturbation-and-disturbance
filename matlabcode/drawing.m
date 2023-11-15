close all;

figure(1);
subplot(211);
plot(t,qd(:,1),'r',t,y(:,1),'b');
xlabel('time(s)');ylabel('position tracking of link 1');
subplot(212);
plot(t,qd(:,2),'r',t,y(:,3),'b');
xlabel('time(s)');ylabel('position tracking of link 2');

figure(2);
subplot(211);
plot(t,dqd(:,1),'r',t,y(:,2),'b');
xlabel('time(s)');ylabel('velocity tracking of link 1');
subplot(212);
plot(t,dqd(:,2),'r',t,y(:,4),'b');
xlabel('time(s)');ylabel('velocity tracking of link 2');

figure(3);
subplot(211);
plot(t,y(:,1)-qd(:,1),'r');
xlabel('time(s)');ylabel('position tracking error of link 1');
subplot(212);
plot(t,y(:,3)-qd(:,2),'r');
xlabel('time(s)');ylabel('position tracking error of link 2');

figure(4);
subplot(211);
plot(t,y(:,2)-dqd(:,1),'r');
xlabel('time(s)');ylabel('velocity tracking error of link 1');
subplot(212);
plot(t,y(:,4)-dqd(:,2),'r');
xlabel('time(s)');ylabel('velocity tracking error of link 2');

figure(5);
subplot(211);
plot(t,u(:,1),'r');
xlabel('time(s)');ylabel('optimal control of Link1');
subplot(212);
plot(t,u(:,2),'r');
xlabel('time(s)');ylabel('optimal control of Link2');

figure(6);
subplot(211);
plot(t,tol(:,1),'r');
xlabel('time(s)');ylabel('corresponding optimal applied torque of Link1');
subplot(212);
plot(t,tol(:,2),'r');
xlabel('time(s)');ylabel('corresponding optimal applied torque of Link2');
