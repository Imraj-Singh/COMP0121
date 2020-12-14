function Time = unroll_time(N,t)
%UNROLL_TIME by Imraj Singh
% Inputs
%   N,t
% Outputs
%   Time

% Create the time vector
Time = linspace(0, 0, sum(N));
Time(1:N(1)) = linspace(0, t(1), N(1));
for i = 2:length(N)
   Time(sum(N(1:i-1))+1:sum(N(1:i))) = linspace(sum(t(1:(i-1))),sum(t(1:i)),N(i));
end
end

