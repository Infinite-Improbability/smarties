function [stCoa, CstTRa] = slvForTCoated(stParamsCore, stParamsCoat, stOptions)
  %% slvForT
% Calculates the T-matrix and orientation-averaged properties for a coated
% particle
% Input:
%       - stParamsCore and stParamsCoat: structs for core and coating respectively.
%              The following parameters should be defined in stParams:
%              - a: semi-axis along x,y
%              - c: semi-axis along z
%              - k1: wavevector in embedding medium (of refractive index nM) (k1=2*pi*nM/lambda)
%                    For the core this is k0 * nM, where k0 is the
%                    wavevector incident on the coating. Likewise, s for
%                    the coating is nM for the core.
%              - s: relative refractive index of particle (s=n_Particle / nM)
%              TODO: Modify so user doesn't have to calculate s for core?
%              - N: number of multipoles for T-matrix
%              - nNbTheta: number of thetas for quadratures
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
%       - stGeometry (optional): struct with geometry
%       TODO: Restore geometry support
%
% Output:
%       - stCoa: struct with orientation averaged cross-sections
%                Cext, Csca, Cabs
%       - CstTRa: cell {1 x M} of structs defining the T (and possibly R)
%                matrices
%
% Dependency: 
% rvhGetAverageCrossSections, rvhGetSymmetricMat, rvhGetTRfromPQ,
% rvhTruncateMatrices, slvGetOptionsFromStruct, sphCalculatePQ,
% sphEstimateDelta, sphEstimateNB, sphMakeGeometry

% TODO: Adjust N estimation for coated spheroids

% Sanity check
if stParamsCore.k ~= stParamsCoat.k * stParamsCoat.s
    warning("slvForTCoated: stParamsCore.k does not equal stParamsCoat.k * stParamsCoat.s")
end

% Get T1 (T-matrix for core in medium matching coating)
[stCoaCore, CstTRaCore] = slvForT(stParamsCore, stOptions, stGeometry);

% Get P2, Q2
PQ = getPQ(stParamsCoat, false);

% Get PP2, QQ2 by using hankel functions in place of bessel functions
PPQQ = getPQ(stParamsCoat, true);

% Combine to get T-matrix for coated particle


end

function CstPQa = getPQ(stParams, coated)
%% getPQ
% TODO: A lot of this is duplicated from slvForT. Can it be cleaned up?
% Also this func needs proper documentation

stk1s.k1 = stParams.k1;
stk1s.s = stParams.s;

N = stParams.N;
nNbTheta = stParams.nNbTheta;

[~,Delta,NB,absmvec,~, bOutput] = slvGetOptionsFromStruct(stParams,stOptions);

stk1s.bOutput=bOutput;

% Make structure describing spheroidal geometry and quadrature points for
% numerical integrations
stGeometry = sphMakeGeometry(nNbTheta, a, c);

if Delta<0 % then need to estimate Delta
    [Delta, T2211err]= sphEstimateDelta(stGeometry, stk1s);
    if isnan(Delta)
        disp ('ERROR: Delta could not be found. Results are likely to be non-converged. Try choosing Delta manually instead.');
        return;
    end
    disp (['Delta estimated to \Delta=', int2str(Delta),' with relative error in T_{11}^{22,m=1} of ',  num2str(T2211err)]);
end

NQ = N+Delta;% NQ>=N: Maximum multipole order for computing P and Q matrices

CstPQa =  sphCalculatePQ(NQ, absmvec, stGeometry, stk1s, NB, coated);
end

function cellMat = cstMultiply(A, B, typeA, typeB)
%% cstMultiply
% Multiplies two of the matrixes stored in smarties cell structure, as used
% in the variables named Cst__.
% type: type of matrix used, e.g. P or Q

if size(A) ~= size(B)
    warning('slvForTCoated: Attempting to multiply mismatched cells.')
end

cellMat = cellfun(cellMultiply, A, B, typeA, typeB);


end

function cellProd = cellMultiply(cellA, cellB, typeA, typeB)
%% cellMultiply
% Multiples two cells of the Cst__ matrices

typeStr = append('st4M', typeA);
if any(strcmp(cellA.CsMatList, typeStr))
    cellA.type = typeStr;
else
    warning('slvForTCoated: cellA missing matching type of matrix')
end
matricesA = 

typeStr = append('st4M', typeB);
if any(strcmp(cellB.CsMatList, typeStr))
    cellB.type = typeStr;
else
    warning('slvForTCoated: cellB missing matching type of matrix')
end



end