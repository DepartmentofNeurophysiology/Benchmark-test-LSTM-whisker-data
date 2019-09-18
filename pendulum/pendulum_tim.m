function dydt = pendulum_tim(t,theta,u,ut,cw0,w02,n)

gwn=n*randn(1);
u = interp1(ut,u,t);

dydt = [theta(2); -w02*sin(theta(1))-cw0*theta(2)+gwn+u];