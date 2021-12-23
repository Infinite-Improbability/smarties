function X = invertLUcol(B)
%% Invert matrix using LU
% This performs the matrix inversion using LU with partial pivoting of columns
% Note that in Matlab, this is equivalent for near-singular matrices to
% X = transpose(transpose(B)\I), but not to X = I/B which does
% pivoting on rows of B
% It is also equivalent to X=inv(B.').', but not to X=inv(B).
% The code below is more explicit and gives the same results in both Matlab and Octave
    
   [L,U,P] = lu(B.'); % LU decomposition with partial pivoting LU = PB.'
    optsLT.LT=true;
    optsUT.UT=true;
    Y=linsolve(L,P,optsLT); % use lower triangular solver
    X=linsolve(U,Y,optsUT).'; % use upper triangular solver

end