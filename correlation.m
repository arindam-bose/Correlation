function [C, lags] = correlation(X, varargin)
% CORRELATION calculates periodic and aperiodic cross-correlation function estimates.
%
% Inputs:
%       X        : a set of M sequence of length N vectors (N>1).
%       CORRTYPE : 'a' - aperiodic cross-correlation (default)
%                  'p' - periodic cross-correlation
%
% Outputs:
%       C        : M x M x 2*N-1 matrix of cross-correlation values
%                  M x M is the size of cross-correaltion matrix for each lag -N+1 <= l <= N-1
%       k        : corresponding lags
%
% Usage:
%       R        = correlation(X, 'a');
%       C        = correlation(X, 'p');
%       [R,k]    = correlation(X, 'a');
%       [C,k]    = correlation(X, 'p');
%
% Arindam Bose
% Fall 2017

narginchk(1,3);
[corrType] = parseinput(varargin{:});

[N,M] = size(X);

if (strcmp(corrType,'a') || strcmp(corrType,'A'))
    C = zeros(M,M,2*N-1);
    for i = 1 : M
        for j = 1 : M
            for k = 1 : 2*N-1
                C(i,j,k) = sum(X(:,i) .* shift(conj(X(:,j)), k-N));
            end
        end
    end
else                          % corrType = 'p' || corrType = 'P'
    C = zeros(M,M,2*N-1);
    for i = 1 : M
        for j = 1 : M
%             for k = 1 : 2*N-1
%                 C(i,j,k) = sum(X(:,i) .* circshift(conj(X(:,j)), k-N));
%             end
            temp = ifft(fft(X(:,i)) .* conj(fft(X(:,j))));
            C(i,j,:) = [flipud(conj(temp(2:end))); temp];
        end
    end
end
lags = (-N+1) : (N-1);


%--------------------------------------------------------------------------
function [corrType] = parseinput(varargin)
%   Parse the input and determine optional parameters:
%
%   Outputs:
%      corrType  - string with the type of correlation wanted

% Set some defaults:
% Assume aperiodic autocorrelation until mentioned otherwise

corrType = 'a';

switch nargin
   case 1                                 
      if ischar(varargin{1})
         corrType = varargin{1};           
      else                                 % Not recognized
         error('Input argument is not recognized.');
      end
end

if ( ~strcmp(corrType,'a') && ~strcmp(corrType,'A') && ~strcmp(corrType,'p') && ~strcmp(corrType,'P') )
   fprintf('Unknown corrType.  corrType set to ''a''.\n');
   corrType = 'a';
end


%--------------------------------------------------------------------------
function y = shift(x, k)
%   shift the input array to the right or left with k lags
%
%   Inputs:
%       x : the 1-D array
%       k : the amount of shift 
%           k >= 0 : right shift
%           k <  0 : left shift 
%   Outputs:
%      y : the shifted array padded with zero

y = zeros(size(x));
if k >= 0
	y(k+1:end) = x(1:end-k);
elseif k <= -1
	y(1:end+k) = x(-k+1:end);
end
