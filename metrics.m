function [P, I, R, lags] = metrics(X, varargin)
narginchk(1,2);
[corrType] = varargin{1};
[R, lags]  = correlation(X, corrType);
P        = psl(R); 
I        = isl(R);
end