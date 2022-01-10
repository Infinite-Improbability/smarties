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
%              - N: number of multipoles for T-matrix. Should probably
%              match for in and out?
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
if stParamsCore.k1 ~= stParamsCoat.k1 * stParamsCoat.s
    warning("slvForTCoated: stParamsCore.k does not equal stParamsCoat.k * stParamsCoat.s")
end

% Expand some options we'll need later
% TODO: Figure out core or coat here
[bGetR,Delta,~,absmvec,bGetSymmetricT, ~] = slvGetOptionsFromStruct(stParamsCoat,stOptions);
N = stParamsCoat.N;
k1 = stParamsCoat.k1;

%% Get T1 (T-matrix for core in medium matching coating)
[~, TRCore] = slvForT(stParamsCore, stOptions);

%% Get P2, Q2
PQ = getPQ(stParamsCoat, stOptions, false);

%% Get PP2, QQ2 by using hankel functions in place of bessel functions
PPQQ = getPQ(stParamsCoat, stOptions, true);

%% Combine to get T-matrix for coated particle
% We'll get a variable with the right structure we can overwrite.
PQCombined = PQ;

% Then find the combined matrix for each m
for m = 1:length(absmvec)
    suffixes = ["eo" "oe"];
    for sufIndex = 1:2
        suffix = char(suffixes(sufIndex));
        T = TRCore{m}.(['st4MT', suffix]);
        P = PQ{m}.(['st4MP', suffix]);
        Q = PQ{m}.(['st4MQ', suffix]);
        PP = PPQQ{m}.(['st4MP', suffix]);
        QQ = PPQQ{m}.(['st4MQ', suffix]);
        
        [P11, P12, P21, P22] = combineMatrixStructures(P, PP, T);
        [Q11, Q12, Q21, Q22] = combineMatrixStructures(Q, QQ, T);
        
        PQCombined{m}.(['st4MP', suffix]).M11 = P11;
        PQCombined{m}.(['st4MP', suffix]).M12 = P12;
        PQCombined{m}.(['st4MP', suffix]).M21 = P21;
        PQCombined{m}.(['st4MP', suffix]).M22 = P22;
        
        PQCombined{m}.(['st4MQ', suffix]).M11 = Q11;
        PQCombined{m}.(['st4MQ', suffix]).M12 = Q12;
        PQCombined{m}.(['st4MQ', suffix]).M21 = Q21;
        PQCombined{m}.(['st4MQ', suffix]).M22 = Q22;
    end
end

% Get T (and possibly R)
CstTRa = rvhGetTRfromPQ(PQCombined,bGetR);

% If needed, discard higher order multipoles
% (which are affected by the finite size of P and Q)
% N+Delta>=N: Maximum multipole order for computing P and Q matrices
if (N+Delta)>N
    CstTRa = rvhTruncateMatrices(CstTRa, N);
 end
% T and R matrices now go up to N multipoles

% If required, symmetrize the T-matrix
if bGetSymmetricT
    CstTRa = rvhGetSymmetricMat(CstTRa, {'st4MT'});
end

% Calculate the (Ext, Abs, Sca) cross-sections for orientation-averaged excitation
stCoa = rvhGetAverageCrossSections(k1, CstTRa);


end



function CstPQa = getPQ(stParams, stOptions, coated)
%% getPQ
% TODO: A lot of this is duplicated from slvForT. Can it be cleaned up?
% Also this func needs proper documentation

a = stParams.a;
c = stParams.c;

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


function [D11, D12, D21, D22] = combineMatrixStructures(A,B,C)
    %% combineMatrixStructures(A,B)
    % Combines two matrices A and B partitioned into four equal sized
    % blocks (M_row,column) to produce
    % D = A + B * C
    % with the same structure
    
    A11 = A.M11;
    A12 = A.M12;
    A21 = A.M21;
    A22 = A.M22;
    
    B11 = B.M11;
    B12 = B.M12;
    B21 = B.M21;
    B22 = B.M22;
    
    C11 = C.M11;
    C12 = C.M12;
    C21 = C.M21;
    C22 = C.M22;
    
    D11 = A11 + (B11*C11 + B12*C21);
    D12 = A12 + (B11*C12 + B12*C22);
    D21 = A21 + (B21*C11 + B22*C21);
    D22 = A22 + (B21*C12 + B22*C22);
    
end