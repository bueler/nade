function Z = mysum(N)
% MYSUM Approximates an infinite series by its Nth partial sum:
%    sum_{n=1}^N sin(n) / (n^3 + 1)
% Example:
%   >> mysum(100)

nn = 1:N;
Z = sum(sin(nn) ./ (nn.^3 + 1));
