function [conv_whisk, conv_curve] = convolve_kernel_whisk_curve( KernelStruct, trials, dat)
% Function convolves whisker and curvature arrays with kernels, for correct
% trials with licking. Returns this as two strcutures with each entry in
% the structure 200 (neurons) time arrays for each trial.

la = sum(KernelStruct{1,1}.kerneltime<0);
lc = sum(KernelStruct{1,1}.kerneltime>0);
[Nkernel, Ndimk] = size(KernelStruct{1,1}.Kernels);
ConvTrace = cell(Nkernel, 1);
trial=1;
for i=1:trials
    if dat(i).correct == 1 
        if dat(i).lick_times ~= 0
            trial = trial + 1;
            for nk = 1: Nkernel
                conv_whisk(trial).ConvTrace{nk,1} = convolve_kernel_acausal( dat(i).thetaVec, KernelStruct{1,1}.Kernels{nk, 1},la,lc);
                conv_curve(trial).ConvTrace{nk,1} = convolve_kernel_acausal( dat(i).kappaVec, KernelStruct{1,1}.Kernels{nk, 2},la,lc);
            end
        end
    end
end

end

