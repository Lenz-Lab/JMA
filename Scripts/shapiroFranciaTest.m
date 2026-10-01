function [is_normal, pvalue, W] = shapiroFranciaTest(x, alpha)
% Shapiro-Francia normality test.
% Same calculation as row 8 of normalitytest.m, pulled out on its own so
% JMA_03 does not compute the other nine tests it never uses.
%
% Inputs:
%   x         : data vector (at least 5 values)
%   alpha     : significance level (default 0.05, matching normalitytest.m)
%
% Outputs:
%   is_normal : 1 if p > alpha (fail to reject normality), otherwise 0
%   pvalue    : p-value of the test
%   W         : Shapiro-Francia test statistic
%
% Reference:
%   Oner, M., & Deveci Kocakoc, I. (2017). JMASM 49: A Compilation of Some
%   Popular Goodness of Fit Tests for Normal Distribution: Their Algorithms
%   and MATLAB Codes (MATLAB). Journal of Modern Applied Statistical
%   Methods, 16(2), 30.

if nargin < 2
    alpha = 0.05;
end

x = x(:)';
n = length(x);
y = sort(x);

mi      = norminv(((1:n)-0.375)/(n+0.25));
weights = mi./sqrt(mi*mi');
W       = sum(y.*weights)^2/sum((y-mean(y)).^2);

u1      = log(log(n))-log(n);
u2      = log(log(n))+2/log(n);
mu      = -1.2725+1.0521*u1;
sigma   = 1.0308-0.26758*u2;
pvalue  = 1-normcdf((log(1-W)-mu)/sigma,0,1);

is_normal = double(pvalue > alpha);
end
