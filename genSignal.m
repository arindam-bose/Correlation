function [X] = genSignal(ch_X, params)
% GENSIGNAL generates a set of M sequences of length N according to the following choices
% 
% Inputs:
%       ch_X :       1: random complex sequence a + ib (a,b are in [-1,1])
%                    2: binary sequence [1,-1]
%                    3: all ones 
%                    4: finite alphabate: rth roots of unity
%                    5: single r-th root ZadoffChu sequence
%                    6: chirp signal
%              default: random real sequence [0,1]
%       params.M      : number of sequences in the family
%       params.N      : length of each sequence
%       params.R      : R-th root of unity (needed only when ch_X is 4 and 5)
%
% Outputs:
%       X             : a set of M sequences of length N
%
% Usage:
%       X = genSignal(ch_X, params);
%
% Arindam Bose
% Fall 2017

narginchk(2,2);

emsg = '';
if isfield(params,'M')
    M = params.M;
end
if isfield(params,'N')
    N = params.N;
end
if isfield(params,'R')
    R = params.R;
end

if ch_X == 1                                         % random complex sequence a + ib (a,b are in [-1,1])
    if exist('M', 'var') && exist('N', 'var')
        X = (-1 + 2*rand(N,M)) + 1i*(-1 + 2*rand(N,M));
    else
        emsg = 'variables M,N are required in params';
    end
elseif ch_X == 2                                     % binary sequence [1,-1]
    if exist('M', 'var') && exist('N', 'var')
        X = randi([0, 1], N, M);
        X(X == 0) = -1;
    else
        emsg = 'variables M,N are required in params';
    end
elseif ch_X == 3                                     % all ones
    if exist('M', 'var') && exist('N', 'var')
        X = ones(N, M);
    else
        emsg = 'variables M,N are required in params';
    end
elseif ch_X == 4                                     % finite alphabate: R-th roots of unity
    if exist('M', 'var') && exist('N', 'var') && exist('R', 'var')
        p = [1, zeros(1,R-1), -1];
        z = roots(p);
        X = z(randi(length(z), N, M));
    else
        emsg = 'variables M,N,R are required in params';
    end
elseif ch_X == 5                                     % single R-th root ZadoffChu sequence
    if exist('M', 'var') && exist('N', 'var') && exist('R', 'var')
        if M > 1
            error('M cannot be greater than 1');
        end
        X = lteZadoffChuSeq(R,N);
    else
        emsg = 'variables M,N,R are required in params';
    end
elseif ch_X == 6                                    % chirp signal with bandwidth B, N = BT
    if exist('M', 'var') && exist('N', 'var')
        if M > 1
            error('M cannot be greater than 1');
        end
        n = 1:1:N;
        X = exp(1i*pi.*n.^2./N);
        X = X';
    else
        emsg = 'variables M,N,b are required in params';
    end
else                                                % random real sequence [0,1]
    disp('ch_X unrecognized, switching to default');
    if exist('M', 'var') && exist('N', 'var')
        X = rand(N, M);
    else
        emsg = 'variables M,N are required in params';
    end
end
error(emsg);
