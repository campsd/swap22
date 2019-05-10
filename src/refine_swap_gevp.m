function [A,B,Qupt,Zupt] = refine_swap_gevp(A, B)
% REFINE_SWAP_GEVP -- for use with swap22_gevp
% A,B is at this stage nearly block upper triangular. We compute a
% transformation to reduce the norm of the off-diagonal blocks.
%
% Software based on: "Swapping 2x2 blocks in the Schur and generalized
%                     Schur form" - D. Camps, N. Mastronardi, R. Vandebril,
%                     and P. Van Dooren.
I = eye(2,2);
err = [1 1]; % initialize to 1, only called when error is too large
Qupt = eye(4,4);
Zupt = eye(4,4);
maxit = 20; nit = 0;
while (err(1) > eps || err(2) > eps) && (nit < maxit)
  nit = nit + 1;
  M = [kron(I,A(3:4,3:4)) -kron(transpose(A(1:2,1:2)),I); ...
       kron(I,B(3:4,3:4)) -kron(transpose(B(1:2,1:2)),I)];
  x = M\[A(3:4,1); A(3:4,2); B(3:4,1); B(3:4,2)];
  X = reshape(x(1:4),2,2);
  Y = reshape(x(5:8),2,2);
  [UX,SX,VX] = svd(X);SX=diag(SX);
  [COSX,SINX] = cosine_sv(SX);
  [UY,SY,VY] = svd(Y);SY=diag(SY);
  [COSY,SINY] = cosine_sv(SY);
  Qup = [VY*diag(COSY)*VY' VY*diag(SINY)*UY' ; ...
        -UY*diag(SINY)*VY' UY*diag(COSY)*UY'] ;
  Zup = [VX*diag(COSX)*VX' VX*diag(SINX)*UX' ; ...
        -UX*diag(SINX)*VX' UX*diag(COSX)*UX'] ;                                    
  A = Qup.'*A*Zup;
  B = Qup.'*B*Zup;
  err=[norm(A(3:4,1:2))/norm(A),norm(B(3:4,1:2))/norm(B)];
  Qupt = Qup.' * Qupt;
  Zupt = Zupt * Zup;
end

end