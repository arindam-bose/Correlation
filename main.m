%%
clc;
clear all;
close all;

%% Single sequence generation
% X: a set of M sequences of length N
M         = 1;
N         = 139;
r         = 25;           % r-th root of unity

idx1      = 1;           % which pair you want to see the correlations
idx2      = 1;

normed    = 1;           % 0: un-normalized axis, 
                         % 1: normalized axis
ch_X      = 5;           % 1: random complex sequence a + ib (a,b are in [-1,1])
                         % 2: binary sequence [1,-1]
                         % 3: all ones 
                         % 4: finite alphabate: rth roots of unity
                         % 5: single r-th root ZadoffChu sequence

if (idx1 > M || idx2 > M)
    error('sequence indices cannot be larger than M');
end

X = genSignal(ch_X, M, N, r);

% cross-correlations
R    = correlation(X, 'a');
C    = correlation(X, 'p');

% plots
k = (-N+1) : (N-1);
if normed == 1
    temp_d   = (norm(X(:,idx1))*norm(X(:,idx2)));
    temp_str = 'Normalized amplitude'; 
elseif normed == 0
    temp_d   = 1;
    temp_str = 'Un-normalized amplitude'; 
end

figure(1); 
subplot(2,2,1); 
                plot(k, abs(squeeze(R(idx1,idx2,:))) / temp_d);
                title(sprintf('Aperiodic cross-correlation X(:,%d) and X(:,%d)', idx1, idx2));
                xlabel('k'); ylabel(temp_str); grid on;
subplot(2,2,2); 
                plot(k, abs(squeeze(R(idx2,idx1,:)))/ temp_d);
                title(sprintf('Aperiodic cross-correlation X(:,%d) and X(:,%d)', idx2, idx1));
                xlabel('k'); ylabel(temp_str); grid on;
subplot(2,2,3); 
                plot(k, abs(squeeze(C(idx1,idx2,:))) / temp_d);
                title(sprintf('Periodic cross-correlation X(:,%d) and X(:,%d)', idx1, idx2)); 
                xlabel('k'); ylabel(temp_str); grid on;
subplot(2,2,4); 
                plot(k, abs(squeeze(C(idx2,idx1,:))) / temp_d);
                title(sprintf('Periodic cross-correlation X(:,%d) and X(:,%d)', idx2, idx1));
                xlabel('k'); ylabel(temp_str); grid on;

%% psl and isl against sequence length
% X: a set of M sequences of maximum lengths N_max
M         = 5;
N_max     = 1000;
r         = 5;           % r-th root of unity
ch_X      = 1;           % 1: random requence
                         % 2: binary sequence
                         % 3: all ones 
                         % 4: finite alphabate: roots of unity
                         % 5: single r-th root ZadoffChu sequence

Apsl      = zeros(1, N_max); 
Aisl      = zeros(1, N_max);
Ppsl      = zeros(1, N_max);
Pisl      = zeros(1, N_max);
welchisl  = zeros(1, N_max);
welchApsl = zeros(1, N_max);
welchPpsl = zeros(1, N_max);

for n = 2:N_max
    X            = genSignal(ch_X, M, n, r);
    R            = correlation(X, 'a');
    C            = correlation(X, 'p');
    Apsl(n)      = psl(R); 
    Aisl(n)      = isl(R);
    Ppsl(n)      = psl(C); 
    Pisl(n)      = isl(C);
    sigma        = norm(X,'fro')^2/M;
    welchApsl(n) = sigma * sqrt((M-1)/(2*M*n-M-1));
    welchPpsl(n) = sigma * sqrt((M-1)/(M*n-1));
    welchisl(n)  = sigma^2 * M *(M-1);
end

% plots
figure(2); 
subplot(2,2,1); 
                plot(Apsl);  hold on;
                plot(welchApsl);  
                title('Aperiodic PSL growth'); legend('PSL', 'Welch bound');
                xlabel('Sequence length'); grid on; hold off;
subplot(2,2,2); 
                plot(Aisl); hold on;
                plot(welchisl);
                title('Aperiodic ISL growth'); legend('ISL', 'Welch bound');
                xlabel('Sequence length'); grid on;  hold off;
subplot(2,2,3); 
                plot(Ppsl); hold on;
                plot(welchPpsl);
                title('Periodic PSL growth'); legend('PSL', 'Welch bound');
                xlabel('Sequence length'); grid on;  hold off;
subplot(2,2,4); 
                plot(Pisl); hold on;
                plot(welchisl);
                title('Periodic ISL growth'); legend('ISL', 'Welch bound');
                xlabel('Sequence length'); grid on;  hold off;

                
%% frequency and phase distribution of a signle sequence
N    = 1001;
r    = 25;               % r-th root of unity

ch_x = 5;                % 1: random complex sequence a + ib (a,b are in [-1,1])
                         % 2: binary sequence [1,-1]
                         % 3: all ones 
                         % 4: finite alphabate: rth roots of unity
                         % 5: single r-th root ZadoffChu sequence
                         
x    = genSignal(ch_x, 1, N, r);
X    = fftshift(fft(x,N));

% plots
figure(3);
subplot(2,1,1); plot(20*log10(abs(X))); title('Magnitude of the frequency spectrum');
                xlabel('Frequency (N)'); ylabel('magnitude (dB)'); grid on;
subplot(2,1,2); plot(angle(X)/pi); title('Phase of the frequency spectrum');
                xlabel('Frequency (N)'); ylabel('phase / \pi'); grid on;            