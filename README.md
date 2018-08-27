# What is this package?
This package calculates the periodic and aperiodic cross-correlation matrices of a set of sequences X of M-members, each of length N.
It returns a 3D matrix of size M x M x 2*(N-1): one M x M matrix for each lag from -(N-1) to (N-1).
It also calculates the peak sidelobe level (PSL) and integrated sidelobe level (ISL) of the set of sequences X.

# Funtion descriptions
# 1. CORRELATION
CORRELATION calculates periodic and aperiodic cross-correlation function estimates.
1. Inputs:
    X: a set of M sequence of length N vectors (N>1).
    CORRTYPE: 'a' - aperiodic cross-correlation (default), 'p' - periodic cross-correlation

2. Outputs:
    C: M x M x 2*N-1 matrix of cross-correlation values, M x M is the size of cross-correaltion matrix for each lag -N+1 <= l <= N-1
    k: corresponding lags

3. Usage:
R        = correlation(X, 'a');
C        = correlation(X, 'p');
[R,k]    = correlation(X, 'a');
[C,k]    = correlation(X, 'p');

# 2. PSL
PSL calculates the peak sidelobe level for periodic or aperiodic cross-correlation.

1. Inputs:
       R : M x M x 2*N-1 matrix of cross-correlation values.

2. Outputs:
       P : max([R(u,v,k)_{u~=v,k}] U [R(u,v,k)_{u=v,k~=0}])

3. Usage:
       P = psl( R );
       
# 3. ISL   
ISL calculates the integrated sidelobe level for periodic or aperiodic cross-correlation.

1. Inputs:
       R : M x M x 2*N-1 matrix of cross-correlation values.

2. Outputs:
       I : sum([R(u,v,k)_{u~=v,k}]^2 + [R(u,v,k)_{u=v,k~=0}]^2)

3. Usage:
       I = isl( R );
