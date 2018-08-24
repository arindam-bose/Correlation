function [P] = psl(R)
% PSL calculates the peak sidelobe level for periodic or aperiodic cross-correlation.
%
% Inputs:
%       R : M x M x 2*N-1 matrix of cross-correlation values.
%
% Outputs:
%       P : max([R(u,v,k)_{u~=v,k}] U [R(u,v,k)_{u=v,k~=0}])
%
% Usage:
%       P = psl(R);
%
% Arindam Bose
% Fall 2017

Ld   = size(R,3);
N    = (Ld+1)/2;
Ra   = abs(R);

maxL = max(max(max (Ra(:,:, 1:N-1))));
maxM = max(max(tril(Ra(:,:,  N),-1)));

P    = max(maxL, maxM);