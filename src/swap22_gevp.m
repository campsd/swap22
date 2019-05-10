function [ A, B, Q, Z ] = swap22_gevp( A, B, qrup )
% Swap 2x2 blocks in generalized real Schur form.
%
% (A,B) is expected to have complex-conjuagete eigenvalues.
% Input
% -----
%   A  :  4x4 block upper triangular having 2x2 blocks
%   B  :  4x4 block upper triangular having 2x2 blocks
%   qrup :  boolean, if true, the transformations are computed via a QR
%         decomposition, otherwise via an SVD.
%
% Output
% ------
%   A  :  4x4 block upper triangular with swapped eigenvalues
%   B  :  4x4 block upper triangular with swapped eigenvalues
%   Q  :  equivalence transformation
%   Z  :  equivalence transformation
%
% Software based on: "Swapping 2x2 blocks in the Schur and generalized
%                     Schur form" - D. Camps, N. Mastronardi, R. Vandebril,
%                     and P. Van Dooren.
warning('off','MATLAB:nearlySingularMatrix');
I = eye(2); 
M = [kron(I,A(1:2,1:2)) kron(transpose(A(3:4,3:4)),I); ...
     kron(I,B(1:2,1:2)) kron(transpose(B(3:4,3:4)),I)];
x = M\[-A(1:2,3); -A(1:2,4); -B(1:2,3); -B(1:2,4)];
X = reshape(x(1:4),2,2);
Y = reshape(x(5:end),2,2);
% COMPUTE ORTHONORMAL BASIS FOR LEFT AND RIGHT DEFLATING SUBSPACES
if qrup == 1 % Use QR decomposition
  L = [eye(2,2) Y];
  R = [X; eye(2,2)];
  [Q,~] = qr(L(2:-1:1,4:-1:1)'); Q = Q(end:-1:1, end:-1:1)' ;
  [Z,~] = qr(R);
else % Use SVD
  [UX,SX,VX] = svd(X);SX=diag(SX);
  [COSX,SINX] = cosine_sv(SX);
  [UY,SY,VY] = svd(Y);SY=diag(SY);
  [COSY,SINY] = cosine_sv(SY);
  Q = [-VY*diag(SINY)*UY' VY*diag(COSY)*VY'; ...
       UY*diag(COSY)*UY' UY*diag(SINY)*VY'] ;
  Z = [UX*diag(SINX)*VX' UX*diag(COSX)*UX' ; ...
       VX*diag(COSX)*VX' -VX*diag(SINX)*UX'] ;                                
end  
A = Q*A*Z;
B = Q*B*Z;
err=[norm(A(3:4,1:2))/norm(A),norm(B(3:4,1:2))/norm(B)];
if (err(1) > eps || err(2) > eps)
  [A,B,Qup,Zup] = refine_swap_gevp(A,B);
  Q = Qup * Q;
  Z = Z * Zup;
end
A(3:4,1:2) = 0; B(3:4,1:2) = 0;
warning('on','MATLAB:nearlySingularMatrix');  
end
