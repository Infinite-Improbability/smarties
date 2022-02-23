function [Nreturn, err] = coaEstimateN(stParamsCore, stParamsCoat, stOptions, maxAcc)
%% sphEstimateN
% Estimates the required number of multipoles N for convergence of physical properties
%
% coaEstimateN(stParams, stOptions, maxAcc)
% finds an estimate for N by studying the convergence of the
% the orientation-averaged extinction cross-section calculated for m=0 and
% m=1 only to speed up test.
% Also returns the estimated converged precision.
%
% Input:
%       - stParams:   Structure containing simulation parameters.
%                   If stParams.N is defined, then the
%                   convergence test only checks for values larger than
%                   this.
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
%       - maxAcc: desired relative accuracy (default 1e-20). If this cannot
%               be reached, the best possible accuracy is returned.
%
% Dependency:
% sphEstimateN, sphMakeGeometry, slvForTCoated


%set default parameters
if nargin < 3
    maxAcc = 1e-20;
end
minAcc = 1e-3;
maxAcc = min(maxAcc,minAcc);

absmvec = [0,1]; % only m=0 and 1

% This works on only one wavelength, so we choose the largest k1 * s
% as representative of the worst case
% Find max and min relative refractive index
[~,ind] = max(abs(stParamsCoat.k1 .* stParamsCoat.s));
stParamsCoat.s = stParamsCoat.s(ind); stParamsCore.s = stParamsCore.s(ind); 
stParamsCoat.k1 = stParamsCoat.k1(ind); stParamsCore.k1 = stParamsCore.k1(ind);

% Calculate parameters for core
% We can do this easily by using the existing functions for homogenous spheroids
[coreN, coreErr] = sphEstimateN(stParamsCore, stOptions, maxAcc);
stParamsCore.N = coreN;
Nmin = coreN;

if isfield(stParamsCoat,'N') % if defined, only look for larger values
    Nmin=max(Nmin,stParamsCoat.N);
end

NQarr=[21,51,101,151,201]; % Arrays of NQs that will be used for tests
NQarr=NQarr(NQarr>Nmin);
nqmax=length(NQarr);

% Number of points used for linear fit to check if error has reached a
% plateau
NforConv = 4;
Afit = [ones(NforConv,1), (1:NforConv).'];

warning('off', 'SMARTIES:missingm'); % suppress warnings in rvhGetAverageCrossSections
for nq=1:nqmax % loop on NQ
    NQ=NQarr(nq);

    % Store errors
    Ntest=Nmin:2:NQ; % Test by steps of 2
    nnmax=length(Ntest);
    Qerr = zeros(nnmax,1);
    Q=1e100;

    for nn = 1:nnmax % Loop over truncation
        % Truncate P,Q to N and get corresponding T
        N=Ntest(nn);

        % Update supplied options with tested params
        stOptions.absmvec = absmvec;
        stParamsCoat.N = N;

        [stCoa, ~] = slvForTCoated(stParamsCore, stParamsCoat, stOptions);
        Qnew=stCoa.Cext;

        % Relative error
        Qerr(nn) = abs (Q./Qnew-1);
        if nn>1
%            fprintf('nn=%d, N= %d, err= %e\n',nn,N,Qerr(nn));
            if Qerr(nn)<maxAcc % then desired convergence was achieved
                Nreturn = Ntest(nn-1);
                err = Qerr(nn);
                warning('on', 'SMARTIES:missingm'); % reactivate warnings in rvhGetAverageCrossSections
                return;
            end
        end

        % Now check if plateau has been reached
        % Minimum requirement for convergence testing
        if nn > NforConv && max(Qerr((nn-NforConv+1):nn))<minAcc
            % Then test if the last NforConv errors have flattened out by a
            % linear regression to the last NforConv points
            coef=Afit \ log10(abs(Qerr((nn-NforConv+1):nn)));
            slope = coef(2);
            % if slope<0, then still converging
            if slope>=0 % this means no longer better over NforConv steps
                Nreturn = Ntest(nn-NforConv);
                err = mean(Qerr((nn-NforConv+1):nn));
                warning('on', 'SMARTIES:missingm'); % reactivate warnings in rvhGetAverageCrossSections
                return;
            end
        end
        % Prepapre for next step
        Q=Qnew;

    end
    % Convergence was not found, continue loop with a larger NQ
    Nmin=NQ-2*(NforConv-1);
end
%disp('Convergence for N was not found..')
Nreturn= NaN;
err=NaN;
warning('on', 'SMARTIES:missingm'); % reactivate warnings in rvhGetAverageCrossSections
