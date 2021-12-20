function [Q] = exportQmatrix(stQ, complete, invert, out, format)
%% exportQmatrix
% Reshaping to long format and exporting Q-matrix entries to a text file
%
% PARAMETERS:
% - stQ: structure containing Q-matrix elements, as returned by the program
% - complete: [logical] compute negative m's
% - invert: [logical] export Q^-1 matrix
% - out: optional output filename
% - format: string format for text output
%
% RETURNS: Q-matrix elements and indices consolidated in a single matrix
% The output format consists of 8 columns:
% s sp m mp n np Qr Qi
% whereby
% * m  1st m-index
% * mp 2nd m-index (identical, due to rotational symmetry)
% * n  1st n-index
% * np 2nd n-index
% * s  1st block index (electric/magnetic)
% * sp 2nd block index (electric/magnetic)
% * Qr real (Q_s,sp,m,mp,n,np)
% * Qi imag (Q_ssp,m,mp,n,np)
%
% Dependency:
% combine_oeeo

% The input stQ structure is a cell array where each cell is associated
% with a single m-index. Each cell contains Qeo and Qoe structs (in
% addition to P structs, which we shall ignore). 
%
% The Q (eo or oe) structs contain a matrix in four blocks
%
%     Q11(n,np) | Q12(n,np)
%     ----------+----------
%     Q21(n,np) | Q22(n,np)
%
% and vectors ind1 and ind2 indicating which row and column indices are
% contained in each block.
%
% To export the stQ structure we flatten everything out, first for each m,
% then merge all the per-m matrices together.
% For a given m, we start by merging the Qeo and Qoe matrices into one
% using combine_oeeo.
% Then we 
end