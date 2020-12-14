function s_filter = unroll_signal(N, signal)
%UNROLL_SIGNAL by Imraj Singh
% Inputs
%   M, signal
% Outputs
%   s_filter

% Initialise array
s_filter = [];


for i=1:length(N)
    % Create the binary filter to pick out sample values
    s_filter = [s_filter, ones(1,N(i))*signal(i)];
    
    % Miss out the last element of the signal block so the k is kept
    % whilst ensuring 2^n points sampled for the fft
    if i>1
        s_filter(sum(N(1:i)))=0;
    end
    
end
end

