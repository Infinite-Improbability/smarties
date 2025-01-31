function [N, nNbTheta, err] = sphEstimateNandNT(stParams, stOptions, maxAcc, stParamsCore)
%% sphEstimateNandNT
% Estimates the required number of multipoles N and quadrature points nNbTheta
%
% sphEstimateNandNT(stGeometry, stParams, NQmax, acc)
% combines calls to sphEstimateNbTheta and sphEstimateN to obtain estimates
% of these two parameters.
%
% Input:
%       - stParams: Structure containing simulation parameters.
%                   If stParams.N is defined, then the
%                   convergence test only checks for values larger than
%                   this. If working with a coated spheroid this should be
%                   the parameters for the coat.
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
%       - maxAcc: (optional) desired relative accuracy (default 1e-20). If
%               this cannot be reached, the best possible accuracy is
%               returned.
%       - stParamsCore: (otpional) struct with parameter for core, same
%       format as stParams.
%
% Dependency:
% sphEstimateN, sphEstimateNbTheta, coaEstimateN

%set default parameters
coated = false;
if nargin < 3
    maxAcc = 1e-20;
elseif nargin == 4
    coated = true;
end

warning('off');
% First estimate the number of quadrature point for a small value
% of N (which is much faster)
N=5;
stParams.N=N;
stParams.nNbTheta=0;
[nNbTheta,errTh]=sphEstimateNbTheta(stParams,stOptions,maxAcc);


if isnan(nNbTheta) % not converged
    N=NaN;
    err=NaN;
else
    % Then estimate the number of required N for this nNbTheta
    stParams.N=1;
    stParams.nNbTheta=nNbTheta+50; % We add 50 here to avoid problems with small values (occurring when h close to 1)
    if coated
        stParamsCore.nNbTheta = nNbTheta + 50;
        [N, errN] = coaEstimateN(stParamsCore, stParams, stOptions, maxAcc);
    else
        [N,errN] = sphEstimateN(stParams,stOptions,maxAcc);
    end
    if ~isnan(N)
        err=max(errN,errTh);
        % Finally see whether we need more theta's for this N
        stParams.N=N;
        stParams.nNbTheta=nNbTheta;
        [nNbTheta,errTh]=sphEstimateNbTheta(stParams,stOptions,maxAcc);
        if isnan(nNbTheta) % not converged
            N=NaN;
            err=NaN;
        else
            err=max(errTh,err);
        end
    else
        err=NaN;
    end
end
warning('on');
