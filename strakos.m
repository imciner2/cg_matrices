function A = strakos(N, p, lambda_min, kappa, seed)
%STRAKOS Generate a matrix with eigenvalues distributed in an interval
%
% Generate a matrix of size n by n that has a condition number of kappa
% and eigenvalues in the range [lambda_min, kappa*lambda_min]. The eigenvalues
% are distributed in the range according to the parameter p, where p=1
% creates a uniform distribution and p=0 clusters all eigenvalues at lambda_min.
% Values in between create eigenvalues distributed across the range, with a cluster
% at lambda_min.
%
% Parameters:
%    N - Size of the matrix
%    p - Clustering parameter (must be in [0,1])
%    lambda_min - Minimum eigenvalue of the matrix
%    kappa - Condition number of the matrix
%
% Based on the test matrix from:
% Z. Strakos, ‘On the Real Convergence Rate of the Conjugate Gradient Method’,
% Linear Algebra and its Applications, vol. 154–156, pp. 535–549, 1991,
% doi: 10.1016/0024-3795(91)90393-B.
%

% Created by: Ian McInerney
% Created on: August 17, 2022
% SPDX-License-Identifier: MIT

lambda_max = kappa * lambda_min;

i = 2:1:N-1;

% Generate the eigenvalues
lambda = zeros(N, 1);
lambda(1) = lambda_min;
lambda(end) = lambda_max;
lambda(i) = lambda_min + (i - 1)./(N-1).*(lambda_max - lambda_min).*p.^(N-i);

% Generate a set of basis vectors, from https://nhigham.com/2020/04/22/what-is-a-random-orthogonal-matrix/
stream = RandStream( 'mt19937ar', 'Seed', seed );
[Q, R] = qr( rand( stream, N ) );
Q = Q*diag( sign( diag(R ) ) );

A = Q*diag( lambda )*inv(Q);

end
