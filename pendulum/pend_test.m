cw0=0.1;
w02=10; 
n=30;
tspan=[0 15];
init=[pi,0];
ut=linspace(tspan(1),tspan(2));
u=zeros(length(ut),1);%0.5*sin(sqrt(w02)*ut);
sol = ode45(@(t,theta) pendulum_tim(t,theta,u,ut,cw0,w02,n), tspan,init);
%%
t = linspace(tspan(1),tspan(2),400);
theta= deval(sol,t);
%%
theta=theta';
figure(1)
plot(t,theta(:,1))
title('Solution of pendulum');
xlabel('Time t');
ylabel('Theta');
legend('theta') 

y=-cos(theta(:,1));
x=sin(theta(:,1));

figure(2)
plot(t,x)
xlabel('Time t');
ylabel('x coordinate');
figure(3)
plot(t,y)
xlabel('Time t');
ylabel('y coordinate');

%%
L=0.5;
O=[0 0];
axis(gca,'equal')
axis([-0.7 0.7 -0.7 0.7])
p=tspan(2)/length(t);       % Time between each 'frame' 
grid on;
for i=1:length(t)
    P=L*[ x(i) y(i)];
    O_circ= viscircles(O,0.01);
    pend=line([O(1) P(1)],[O(2) P(2)]);
    ball=viscircles(P,0.05);
    pause(p);
    
    if i<length(t)
        delete(pend);
        delete(ball);
        delete(O_circ);
    end
end

    
    
