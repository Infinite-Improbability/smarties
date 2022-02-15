function [PQcells, PPQQcells] = coaCalculatePQ(nMax, absmvec, stGeometry, stParams)
%coaCalculateTMatrix - Description
%
% Syntax: T = coaCalculateTMatrix(input)
%
% Long description

% TODO: Use arguments block

s = stParams.s; % relative refractive index
absmvec = transpose(absmvec); % transposing the vector so we can easily
% iterate over it in the for loop argument

% Gaussian quadrature is provided as geometry

% Calculate n and n(n+1)
nVec = 1:(nMax+1);
nVecProd = nVec .* (nVec+1); % n(n+1)

gm = sqrt((2 * nVec + 1) ./ nVecProd); % An, but where smarties uses
% sqrt((2n+1)/(2n(n+1))), this is sqrt((2n+1)/(n(n+1)))
% i.e. An = gm / sqrt(2)

PQcells = cell(1, length(absmvec));
PPQQcells = cell(1, length(absmvec));

% LISA loops from 0 to stParams.n, but we have absmvec for that
% m is the multipole we are computing at.
% TODO: Vectorise more? Both this loop and interior ones
for m = absmvec
    mMin = max(1, m); % Minimum value of i,j for M_ij to be non-zero
    n = nMax - mMin + 1; % size of the matrices for a given m (since n,k>=m)

    i = mMin:nMax;
    rel = zeros(1, nMax);
    rel(i) = i-mMin+1; % the assignment to a specific slice of rel is important
    % rel gives the position of a submatrix index in the larger matrix
    % the code can probably be reworked to eliminate it
    
    % Now we calculate P, Q matrices for current m
    % The PP and QQ matrices use hankel functions in place of bessel funcs
    % with the arguement s*radius
    % M_ij; i = row index; j = column index
    P = zeros(2*nMax); PP = zeros(2*nMax);
    Q = zeros(2*nMax); QQ = zeros(2*nMax);

    % Loop over theta
    for thetaIndex = 1:length(stGeometry.theta)
        theta = stGeometry.theta(thetaIndex);

        % Don't need to integrate from pi/2 to pi due to symmetry
        % Just double the integral.
        % Since we only go from 0 to pi/2 this leaves Q and P at half their
        % 'correct' values. But since we then multiply T = -P * Q^-1 this
        % cancels out.
        if theta > (pi / 2)
            continue % TODO: Is it safe to break out of the loop entirely?
        end
        % sphMakeGeometry skips over generating these points entirely

        sinTheta = sin(theta);
        sin2Theta = sinTheta ^ 2;
        cosTheta = cos(theta);
        weightedSin = stGeometry.wTheta(thetaIndex) * sinTheta;
        radius = stGeometry.r(thetaIndex); % r(theta)
        dRadius = stGeometry.drdt(thetaIndex); % dr/dtheta
        drSinTheta = dRadius * sinTheta;
        
        % We use Wigner functions in place of Legendre funcs since they
        % don't explode for high n,m
        wig = wigner(sinTheta, cosTheta, m, nMax);

        % spherical bessel functions
        % all are first kind
        % TODO: Consider Riccati-Bessel funcs
        hankel = sqrt(pi/(2*radius)) * besselh((0:nMax)+0.5, radius); % hn(radius)
        bessel = sqrt(pi/(2*s*radius)) * besselj((0:nMax)+0.5, s*radius); % jn(s*radius)
        hankel2 = sqrt(pi/(2*s*radius)) * besselh((0:nMax)+0.5, s*radius); % hn(s*radius)
    
        % Setup a delta function
        % Delta(i)= cos(theta) * Delta(i-1) - i * sin(theta)^2 * Wigner(i)
        % for m=0; when m ~= 0, it is slightly different
        % Also note that i is an integer, not sqrt(-1)
        delta = zeros(1, nMax);
        if m == 0
            delta(1) = -sin2Theta;
            delta(2) = -3*sin2Theta*cosTheta;
            for k = 3:nMax
                delta(k) = delta(k-1)*cosTheta - nVec(k)*sin2Theta*wig(k);
            end
        else
            for k = mMin:nMax
                delta(k) = nVec(k)*cosTheta*wig(k+1) - wig(k)*sqrt((nVec(k)^2 - nVec(m)^2));
            end
        end
    
        % See that P,Q = 0 when (i or j are less than m) and m>1
        % So the loop goes from mMin to nMax, mMin being the minimum value
        % of i,j so that we only calculate non-zero matrix elements
        iz = zeros(1, nMax); jz = zeros(1, nMax);
        ia = zeros(1, nMax); ja = zeros(1, nMax);
        ib = zeros(1, nMax); ka = zeros(1, nMax);
        wig2 = zeros(1, nMax); cb = zeros(1, nMax);
    
        iz(i) = hankel(i+1) .* weightedSin;
        ia(i) = hankel(i) ./ hankel(i+1);
        jz(i) = real(iz(i));
        ja(i) = real(hankel(i)) ./ real(hankel(i+1));
        ib(i) = s .* bessel(i) ./ bessel(i+1);
        ka(i) = ja(i) - ia(i);
        wig2(i) = wig(i+1) .* nVecProd(i) .* drSinTheta;
        cb(i) = s .* hankel2(i) ./ hankel2(i+1);

        for j = i % remember i is a vector
            % Dividing matrices in four block matrices
            %   M11 | M12
            % ------+------
            %   M21 | M22
            % Each are nMax * nMax
            % TODO: Split the P, Q, PP and QQ matrices into 4 (sub) matrices each

            % M12 and M21 first
            % If M=0, II and III quadrant elements equal zero
            if m ~= 0
                % For plane-symmetric particles and i+j=even, matrix elements are zero
                % i = nMin:nMax so if i(1)+j=nMin+j is even, i(1:2:end)+j will be even
                if mod(mMin+j, 2) ~= 0
                    ii = i(1:2:end);
                else
                    ii = i(2:2:end);
                end

                ipj = nVec(ii) .* nVec(j) ./ radius;
                i1 = wig(ii+1) .* wig(j+1) .* drSinTheta;
                gm1 = gm(ii) .* gm(j) / sin2Theta / 2;
                if mod(m,2) ~= 0
                    gm1 = -gm1;
                end
                
                % TODO: Do the nVec calls change anything over using i,j
                % directly?
                za = radius + ia(ii) .* (radius .* ib(j) - nVec(j)) - nVec(ii) .* ib(j) + ipj;
                zb = (delta(ii) .* wig(j+1) + delta(j) .* wig(ii+1)) * radius;
                zc = i1 .* (nVecProd(ii) .* ib(j) + nVecProd(j) .* ia(ii) - ipj .* (nVec(ii) + nVec(j) + 2));
                zr = zc + ka(ii) .* nVecProd(j) .* i1;
                x1 = nVec(m) .* gm1 .* bessel(j+1);

                % M21 calculations
                Q(n+rel(ii),rel(j)) = Q(n+rel(ii),rel(j)) - (iz(ii) .* x1 .* (za .* zb + zc)).';
                zax = za + ka(ii) * (radius .* ib(j) - nVec(j));
                P(n+rel(ii),rel(j)) = P(n+rel(ii),rel(j)) - (jz(ii) .* x1 .* (zax .* zb + zr)).';

                % M12 calculations
                za = za + radius .* (s^2 - 1);
                Q(rel(ii),n+rel(j)) = Q(rel(ii),n+rel(j)) - (iz(ii) .* x1 .* (za .* zb + zc) ./ s).';
                zax = za + ka(ii) .* (radius * ib(j) - nVec(j));
                P(rel(ii),n+rel(j)) = P(rel(ii),n+rel(j)) - (jz(ii) .* x1 .* (zax .* zb + zr) ./ s).';

                % Additional matrices required by coat, PP and QQ
                % Calculate using h(kr) in place of the j(kr) bessel function
                za = radius * (1 + ia(ii) .* cb(j)) - nVec(j) .* ia(ii) - nVec(ii) .* cb(j) + ipj;
                zc = zc + (cb(j) - ib(j)) .* nVecProd(ii) .* i1;
                zr = zc + ka(ii) .* nVecProd(j) .* i1;
                x1 = nVec(m) .* gm1 .* hankel2(j+1);
                % M21
                QQ(n+rel(ii),rel(j)) = QQ(n+rel(ii),rel(j)) - (iz(ii) .* x1 .* (za .* zb + zc)).';
                zax = za + ka(ii) .* (radius .* cb(j) - nVec(j));
                PP(n+rel(ii),rel(j)) = PP(n+rel(ii),rel(j)) - (jz(ii) .* x1 .* (zax .* zb + zc)).';
                % M12
                za = radius .* (s^2 + ia(ii) .* cb(j)) - nVec(j) .* ia(ii) - nVec(ii) .* cb(j) + ipj;
                QQ(rel(ii),n+rel(j)) = QQ(rel(ii),n+rel(j)) - (iz(ii) .* x1 .* (za .* zb + zc) ./ s).';
                zax = za + ka(ii) .* (radius .* cb(j) - nVec(j));
                PP(rel(ii),n+rel(j)) = PP(rel(ii),n+rel(j)) - (jz(ii) .* x1 .* (zax .* zb + zr) ./ s).';
            end

            % M11, M22 calculations
            % For plane-symmetric particles and i+j=odd, matrix elements are zero
            if mod(mMin+j, 2) == 0
                ii = i(1:2:end);
            else
                ii = i(2:2:end);
            end

            gm1 = gm(ii) .* gm(j) / sin2Theta / 2;
            if mod(m,2) ~= 0
                gm1 = -gm1;
            end
            gm2 = complex(0, gm1);

            zd = radius .* (ib(j) - s^2 .* ia(ii)) + s^2 .* nVec(ii) - nVec(j);
            ze = delta(ii) * delta(j);
            if m ~= 0
                ze = ze + (nVecProd(m) - nVec(m)) * wig(ii+1) * wig(j+1);
            end
            ze = ze .* radius;
            zfx = wig2(j) .* delta(ii);
            zfy = wig2(ii) .* delta(j);
            x1 = gm2 .* bessel(j+1);

            % M22
            zf = zfx - s^2 .* zfy;
            Q(n+rel(ii),n+rel(j)) = Q(n+rel(ii),n+rel(j)) + (iz(ii) .* x1 .* (zd .* ze + zf) ./ s).';
            zdx = zd - (ja(ii) - ia(ii)) .* radius .* s^2;
            P(n+rel(ii),n+rel(j)) = P(n+rel(ii),n+rel(j)) + (jz(ii) .* x1 .* (zdx .* ze + zf) ./ s).';
            % M11
            zd = radius .* (ib(j) - ia(ii)) + nVec(ii) - nVec(j);
            zf = zfx - zfy;
            Q(rel(ii),rel(j)) = Q(rel(ii),rel(j)) + (iz(ii) .* x1 .* (zd .* ze + zf)).';
            zdx = zd - (ja(ii) - ia(ii)) * radius;
            P(rel(ii),rel(j)) = P(rel(ii),rel(j)) + (jz(ii) .* x1 .* (zdx .* ze + zf)).';

%             if m == 0 && rel(j) == 1
%                 fprintf('m = %g, Theta Index = %g, Q(29,1) = (%g, %g)\n', m, thetaIndex, real(Q(29,1)), imag(Q(29,1)));
%             end

            % PP and QQ again
            zd = radius .* (cb(j) - s^2 .* ia(ii)) + s^2 .* nVec(ii) - nVec(j);
			zf = zfx - s^2 .* zfy;
			x1 = gm2 .* hankel2(j+1);
	        QQ(n+rel(ii),n+rel(j)) = QQ(n+rel(ii),n+rel(j)) + (iz(ii) .* x1 .* (zd .* ze + zf) ./ s).';
			zdx = zd - (ja(ii) - ia(ii)) .* radius .* s^2;
	        PP(n+rel(ii),n+rel(j)) = PP(n+rel(ii),n+rel(j)) + (jz(ii) .* x1 .* (zdx .* ze + zf) ./ s).';
			zd = radius .* (cb(j) - ia(ii)) + nVec(ii) - nVec(j);
			zf = zfx - zfy;
	        QQ(rel(ii),rel(j)) = QQ(rel(ii),rel(j)) + (iz(ii) .* x1 .* (zd .* ze + zf)).';
			zdx = zd - (ja(ii) - ia(ii)) .* radius;
	        PP(rel(ii),rel(j)) = PP(rel(ii),rel(j)) + (jz(ii) .* x1 .* (zdx .* ze + zf)).';
            
            % testing
            if m == 0
                % disp([theta j Q(60,60)]);
            end
        
        end % end of the loop over j
        

    end % end of the loop over theta

    % Save it in the standard smarties format
    PQcells{1, m+1} = cellExport(P, Q, m, nMax);
    PPQQcells{1, m+1} = cellExport(PP, QQ, m, nMax);


end % end of the loop over m

end


function wig = wigner(sinTheta, cosTheta, m, nMax)
%% wigner
%
% Returns Wigner d-matrix elements for d_{m,0}^{n}(theta)
%
% Inputs:
%   - sinTheta: sin(theta) for a single theta
%   - cosTheta: cos(theta)
%   - m: indice corresponding to projected angular momentum
%   - nMax: maximum multipole order
%
% Outputs:
%   - wig: [nMax+1 x 1] wig(n) = d_{m,0}^{n-1}(theta) for n=1:nMax

% d_{0,0}^{n}(theta) = P_l(cos(theta))

wig = zeros(1, nMax+1);

if m == 0
    wig(1) = 1; % d_{0,0}^{0}(theta) = 1
    wig(2) = cosTheta; % d_{0,0}^{1}(theta) = cos(theta)
    for n = 2:nMax
        wig(n+1) = (wig(n) * (2*n-1) * cosTheta - (n-1) * wig(n-1)) / n;
    end
else
    wig(m+1) = 1;
    if mod(m,2) ~= 0
        wig(m+1) = -wig(m+1);
    end
    for n = 0:(m-1)
        wig(m+1) = (wig(m+1) * sinTheta * sqrt((n+1) * (m+n+1))) / (2 * (n+1));
    end
    for n = (m+1):nMax
        wig(n+1) = (((2*n-1) * cosTheta * wig(n)) - sqrt((n-1)^2 - m^2) * wig(n-1)) / sqrt(n^2 - m^2);
    end
end
    
end


function out = cellExport(P, Q, m, nMax)
    % TODO: doc string

    nMin = max(m, 1);
    n = nMax - nMin + 1;
    i = nMin:nMax;
    rel = i-nMin+1; % unlike creation of rel in main function we don't want a specific slice

    P11 = P(rel, rel);
    P12 = P(rel, rel+n);
    P21 = P(rel+n, rel);
    P22 = P(rel+n, rel+n);

    Q11 = Q(rel, rel);
    Q12 = Q(rel, rel+n);
    Q21 = Q(rel+n, rel);
    Q22 = Q(rel+n, rel+n);

    % This code was taken from sphCalculatePQ
    % The code below uses the symmetry and the even-odd formulation
    evenodd=mod(nMin,2);
    % if nMin is even then the first index (1) is even so need to swap even-odd
    inde=(1+evenodd):2:n;
    indo=(2-evenodd):2:n;

    out.st4MQeo.M12 = Q12(inde,indo);
    out.st4MQeo.M21 = Q21(indo,inde);
    out.st4MQeo.M11 = Q11(inde,inde);
    out.st4MQeo.M22 = Q22(indo,indo);
    out.st4MQeo.m = m;
    out.st4MQeo.ind1 = inde;
    out.st4MQeo.ind2 = indo;

    out.st4MQoe.M12 = Q12(indo,inde);
    out.st4MQoe.M21 = Q21(inde,indo);
    out.st4MQoe.M11 = Q11(indo,indo);
    out.st4MQoe.M22 = Q22(inde,inde);
    out.st4MQoe.m = m;
    out.st4MQoe.ind1 = indo;
    out.st4MQoe.ind2 = inde;

    out.st4MPeo.M12 = P12(inde,indo);
    out.st4MPeo.M21 = P21(indo,inde);
    out.st4MPeo.M11 = P11(inde,inde);
    out.st4MPeo.M22 = P22(indo,indo);
    out.st4MPeo.m = m;
    out.st4MPeo.ind1 = inde;
    out.st4MPeo.ind2 = indo;

    out.st4MPoe.M12 = P12(indo,inde);
    out.st4MPoe.M21 = P21(inde,indo);
    out.st4MPoe.M11 = P11(indo,indo);
    out.st4MPoe.M22 = P22(inde,inde);
    out.st4MPoe.m = m;
    out.st4MPoe.ind1 = indo;
    out.st4MPoe.ind2 = inde;
    
    out.CsMatList={'st4MQ', 'st4MP'};

end