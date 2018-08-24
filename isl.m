function [I] = isl(R)
% ISL calculates the integrated sidelobe level for periodic or aperiodic cross-correlation.
%
% Inputs:
%       R : M x M x 2*N-1 matrix of cross-correlation values.
%
% Outputs:
%       I : sum([R(u,v,k)_{u~=v,k}]^2 + [R(u,v,k)_{u=v,k~=0}]^2)
%
% Usage:
%       I = isl(R);
%
% Arindam Bose
% Fall 2017

Ld   = size(R,3);
N    = (Ld+1)/2;
Ra   = abs(R);

sumL = sum(sum(sum (Ra(:,:, 1:N-1).^2)));
sumM = sum(sum(tril(Ra(:,:,  N),-1).^2));

I    = 2*(sumL + sumM);