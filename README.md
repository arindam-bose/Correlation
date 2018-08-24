# Correlation
CORRELATION calculates periodic and aperiodic cross-correlation function estimates.

Inputs:
       X        : a set of M sequence of length N vectors (N>1).\\
       CORRTYPE : 'a' - aperiodic cross-correlation (default)
                  'p' - periodic cross-correlation

Outputs:
       C        : M x M x 2*N-1 matrix of cross-correlation values
                  M x M is the size of cross-correaltion matrix for each lag -N+1 <= l <= N-1

Usage:
       R        = correlation(X, 'a');
       C        = correlation(X, 'p');

# PSL
PSL calculates the peak sidelobe level for periodic or aperiodic cross-correlation.

Inputs:
       R : M x M x 2*N-1 matrix of cross-correlation values.

Outputs:
       P : max([R(u,v,k)_{u~=v,k}] U [R(u,v,k)_{u=v,k~=0}])

Usage:
       P = psl(R);
       
# ISL   
ISL calculates the integrated sidelobe level for periodic or aperiodic cross-correlation.

Inputs:
       R : M x M x 2*N-1 matrix of cross-correlation values.

Outputs:
       I : sum([R(u,v,k)_{u~=v,k}]^2 + [R(u,v,k)_{u=v,k~=0}]^2)

Usage:
       I = isl(R);
