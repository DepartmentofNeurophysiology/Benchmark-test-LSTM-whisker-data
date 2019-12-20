function   plot_trial( trial_num, dat)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n_t = trial_num; % trial number
y_t = zeros(length(dat(n_t).touch_times));
y_p = zeros(length(dat(n_t).pole_times));
y_l = zeros(length(dat(n_t).lick_times));

figure(1)
subplot(2,1,1);
plot(dat(n_t).timeVec, dat(n_t).thetaVec,'Linewidth', 1.2 )  % Plots theta
hold all
plot(dat(n_t).pole_times,y_p,'r*','Linewidth', 1.5)          % Plots pole times  
if dat(n_t).touch_times ~= 0
plot(dat(n_t).touch_times,y_t,'gx')         % Plots touch times
end
% if dat(n_t).lick_times ~= 0
% plot(dat(n_t).lick_times,y_l,'mo')          % Plots lick times
% end
hold off
legend({'whisking','pole in reach'})
xlabel('time')
ylabel('Whisking angle')



subplot(2,1,2)
plot(dat(n_t).timeVec, dat(n_t).kappaVec,'Linewidth', 1.2)   % Plots curvature
hold on
plot(dat(n_t).pole_times,y_p,'r*','Linewidth', 1.5)          % Plots pole times
if dat(n_t).touch_times ~= 0
hold on
plot(dat(n_t).touch_times,y_t,'gx')         % Plots touch times
end
% if dat(n_t).lick_times ~= 0
% hold on
% plot(dat(n_t).lick_times,y_l,'mo')          % Plots lick times
% end
hold off
xlabel('time')
ylabel('Curvature')
legend({'whisking','pole in reach'})
ylim([-0.01 0.005])


end