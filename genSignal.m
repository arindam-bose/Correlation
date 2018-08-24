function [X] = genSignal(ch_X, varargin)
% GENSIGNAL generates a set of M sequences of length N according to the following choices
% 
% Inputs:
%       ch_X :       1: random complex sequence a + ib (a,b are in [-1,1])
%                    2: binary sequence [1,-1]
%                    3: all ones 
%                    4: finite alphabate: rth roots of unity
%                    5: single r-th root ZadoffChu sequence
%              default: random real sequence [0,1]
%       M : number of sequences in the family
%       N : length of each sequence
%       r : r-th root of unity (optional, needed only when ch_X is 4)
%
% Outputs:
%       X : a set of M sequences of length N
%
% Usage:
%       X = genSignal(ch_X, M, N);
%       X = genSignal(ch_X, M, N, r);  % r needed only when ch_X is 4
%
% Arindam Bose
% Fall 2017

narginchk(3,4);
[M, N, r] = parseinput(ch_X, varargin{:});

if ch_X == 1             % random complex sequence a + ib (a,b are in [-1,1])
    X = (-1 + 2*rand(N,M)) + 1i*(-1 + 2*rand(N,M));
elseif ch_X == 2         % binary sequence [1,-1]
    X = randi([0, 1], N, M);
    X(X == 0) = -1;
elseif ch_X == 3         % all ones
    X = ones(N, M);
elseif ch_X == 4         % finite alphabate: rth roots of unity
    p = [1, zeros(1,r-1), -1];
    z = roots(p);
    X = z(randi(length(z), N, M));
elseif ch_X == 5         % single r-th root ZadoffChu sequence
    if M > 1
        error('M cannot be greater than 1');
    end
    X = lteZadoffChuSeq(r,N);
else                     % random real sequence [0,1]
    X = rand(N, M);
end


%--------------------------------------------------------------------------
function [M, N, r] = parseinput(ch_X, varargin)
%   Parse the input and determine optional parameters:
%
%   Outputs:
%      M : number of sequence
%      N : length of each sequence
%      r : r-th root of unity (needed only when ch_X is 4)

r = 1;
switch nargin
   case 3                              %
      if ch_X == 4
          error('Too few arguments, mention the value of r');
      else
          if isa(varargin{1}, 'double') && isa(varargin{2}, 'double') && ...
                  isscalar(varargin{1}) && isscalar(varargin{2}) %
             M   = varargin{1};
             N   = varargin{2};
          else                                 % Not recognized
             error('M and N should be double scalar');
          end
      end
   case 4
       if isa(varargin{1}, 'double') && isa(varargin{2}, 'double') && isa(varargin{3}, 'double') && ...
              isscalar(varargin{1}) && isscalar(varargin{2}) && isscalar(varargin{3}) %
         M   = varargin{1};
         N   = varargin{2}; 
         r   = varargin{3};
      else                                 % Not recognized
         error('M, N and r should be double scalar');
      end
end