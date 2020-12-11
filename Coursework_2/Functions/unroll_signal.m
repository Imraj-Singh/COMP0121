function s_filter = unroll_signal(Time, N, signal, Mx, My, Mz)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
s_filter = [];


for i=1:length(N)
    s_filter = [s_filter, ones(1,N(i))*signal(i)];
    % Miss out the first element of the signal block so the k is kept
    % whilst ensuring 2^n points sampled for the fft
    if i>1
        s_filter(sum(N(1:(i-1)))+1)=0;
    end
end

end

