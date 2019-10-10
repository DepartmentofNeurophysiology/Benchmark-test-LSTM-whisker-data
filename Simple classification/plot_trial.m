function   plot_trial( trial_num, dat)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n_t=trial_num; % trial number
y_t = zeros(length(dat(n_t).touch_times));
y_p = zeros(length(dat(n_t).pole_times));
y_l = zeros(length(dat(n_t).lick_times));
figure(1)
subplot(2,1,1);
plot(dat(n_t).timeVec, dat(n_t).thetaVec)  % Plots theta
hold on
plot(dat(n_t).pole_times,y_p,'r*')          % Plots pole times  
if dat(n_t).touch_times ~= 0
hold on
plot(dat(n_t).touch_times,y_t,'gx')         % Plots touch times
end
if dat(n_t).lick_times ~= 0
hold on
plot(dat(n_t).lick_times,y_l,'mo')          % Plots lick times
end
hold off
xlabel('time')
ylabel('whisking angle')
legend('whisking','pole in reach')
subplot(2,1,2)
plot(dat(n_t).timeVec, dat(n_t).kappaVec)   % Plots curvature
hold on
plot(dat(n_t).pole_times,y_p,'r*')          % Plots pole times
if dat(n_t).touch_times ~= 0
hold on
plot(dat(n_t).touch_times,y_t,'gx')         % Plots touch times
end
if dat(n_t).lick_times ~= 0
hold on
plot(dat(n_t).lick_times,y_l,'mo')          % Plots lick times
end
hold off
xlabel('time')
ylabel('curvature')
legend('curvature','pole in reach')



end

