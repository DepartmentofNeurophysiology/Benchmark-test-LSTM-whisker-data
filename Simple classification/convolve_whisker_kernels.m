function [ ConvTrace ] = convolve_whisker_kernels( WhiskerStruct, KernelStruct )
%% Convolve recordings with relevant kernel
la = sum(KernelStruct{1,1}.kerneltime<0);
lc = sum(KernelStruct{1,1}.kerneltime>0);
[Ndimw, Ntrace] = size(WhiskerStruct.Recording);
[Nkernel, Ndimk] = size(KernelStruct{1,1}.Kernels);
binsize_kernels = KernelStruct{1,1}.kerneltime(2)-KernelStruct{1,1}.kerneltime(1);
ConvTrace = cell(Nkernel, Ntrace);

for nk = 1:Nkernel
    for nt = 1:Ntrace
        ConvTrace{nk,nt} = nan*ones(Ndimw, length(WhiskerStruct.Recording{nd,nt}));
        for nd = 1:Ndimw  
            ConvTrace{nk,nt}(nd,:) = convolve_kernel_acausal( WhiskerStruct.Recording{nd,nt}, KernelStruct{1,1}.Kernels{nk, nd}, la, lc);
        end
    end
end


end

