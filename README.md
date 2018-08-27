# What is this package?
This package calculates the periodic and aperiodic cross-correlation matrices of a set of sequences X of M-members, each of length N.<br/>
It returns a 3D matrix of size M x M x 2*(N-1): one M x M matrix for each lag from -(N-1) to (N-1).<br/>
It also calculates the peak sidelobe level (PSL) and integrated sidelobe level (ISL) of the set of sequences X.<br/>

# Funtion descriptions
# 1. CORRELATION
CORRELATION calculates periodic and aperiodic cross-correlation function estimates.<br/>
1. Inputs:<br/>
    + X: a set of M sequence of length N vectors (N>1).<br/>
    + CORRTYPE: 'a' - aperiodic cross-correlation (default), 'p' - periodic cross-correlation<br/>

2. Outputs:<br/>
    + C: M x M x 2*N-1 matrix of cross-correlation values, M x M is the size of cross-correaltion matrix for each lag -N+1 <= l <= N-1<br/>
    + k: corresponding lags<br/>

3. Usage:<br/>
    + R        = correlation(X, 'a');<br/>
    + C        = correlation(X, 'p');<br/>
    + [R,k]    = correlation(X, 'a');<br/>
    + [C,k]    = correlation(X, 'p');<br/>

# 2. PSL
PSL calculates the peak sidelobe level for periodic or aperiodic cross-correlation.<br/>

1. Inputs:<br/>
    + R : M x M x 2*N-1 matrix of cross-correlation values.<br/>

2. Outputs:<br/>
    + P : max([R(u,v,k)_{u~=v,k}] U [R(u,v,k)_{u=v,k~=0}])<br/>

3. Usage:<br/>
    + P = psl( R );<br/>
       
# 3. ISL   
ISL calculates the integrated sidelobe level for periodic or aperiodic cross-correlation.<br/>

1. Inputs:<br/>
    + R : M x M x 2*N-1 matrix of cross-correlation values.<br/>

2. Outputs:<br/>
    + I : sum([R(u,v,k)_{u~=v,k}]^2 + [R(u,v,k)_{u=v,k~=0}]^2)<br/>

3. Usage:<br/>
    + I = isl( R );<br/>
