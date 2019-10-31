function [ Xnz_c, Xnz_a ] = standard_nozero( tr1,dat)
%% This function takes one trial, standardises it and removes the zeros returning the whisker angle and curvature
% Select trial 't1'

% Whisker curvature
noNAN = dat(tr1).kappaVec;
noNAN(isnan(noNAN)) = 0;
trial_1 = (noNAN-mean(noNAN))/std(noNAN);

pole_1 = dat('t1').pole_times;
time_1 = dat('t1').timeVec;

X_1 = zeros(time_1(end),1);

i=1;
for t=1:length(X_1)
    X_1(t) = trial_1(i);
    if mod(t,2) == 0
        i = i + 1;
    end
    if i == length(trial_1)
        break
    end
end

% Only input between poles
for t=1:length(X_1)
    if t < pole_1(1) || t > pole_1(2)
        X_1(t) = 0;
    end
end
% Remove zeros
Xnz_c = [];
for t=1:length(X_1)
    if X_1(t) ~= 0
        Xnz_c = [Xnz_c, X_1(t)];
    end
end

% Whisker angle
noNAN = dat(tr1).thetaVec;
noNAN(isnan(noNAN)) = 0;
trial_1 = (noNAN-mean(noNAN))/std(noNAN);

pole_1 = dat('t1').pole_times;
time_1 = dat('t1').timeVec;
%

X_1 = zeros(time_1(end),1);

i=1;
for t=1:length(X_1)
    X_1(t) = trial_1(i);
    if mod(t,2) == 0
        i = i + 1;
    end
    if i == length(trial_1)
        break
    end
end

% Only input between poles
for t=1:length(X_1)
    if t < pole_1(1) || t > pole_1(2)
        X_1(t) = 0;
    end
end
% Remove zeros
Xnz_a = [];
for t=1:length(X_1)
    if X_1(t) ~= 0
        Xnz_a = [Xnz_a, X_1(t)];
    end
end
end

