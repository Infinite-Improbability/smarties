function [Q] = exportQmatrix(stQ, complete, invert, filename, format)
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
% s sp n np m mp Qr Qi
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

%% Parse input parameters
% TODO: Decide between positional or keyword args

arguments
    stQ cell
    complete logical = false % For T matrix this defaults to true
    invert logical = false
    filename string = []
    format string = '%d %d %d %d % d % d %.15g %.15g\n'
end

%% Conversion Method
% The input stQ structure is a cell array where each cell is associated
% with a single m-index. Each cell contains Qeo and Qoe structs (in
% addition to P structs, which we shall ignore). 
%
% The Q (eo or oe) structs contain a matrix in four blocks
%
%     M11(n,np) | M12(n,np)
%     ----------+----------
%     M21(n,np) | M22(n,np)
%
% and vectors ind1 and ind2 indicating which row and column indices are
% contained in each block.
%
% To export the stQ structure we flatten everything out, first for each m,
% then merge all the per-m matrices together.
% For a given m, we start by merging the Qeo and Qoe matrices into one
% using combine_oeeo.
% Then we convert into a matrix with columns
%   s sp n np m mp Qr Qi
% Once we have iterated over all m, we stack the matrices and sort the
% rows.

%% Extract Q-matrix values for each m and list with corresponding indices
% Loop over positive m-values (following storage as independent cells)
mMax = length(stQ);
tmp = cell(1,mMax);

for i = 1:mMax
    % i holds the current value of m in the loop + 1
    % the plus one is because cells are indexed from 1, while m starts from
    % 0.
    
    if invert
        % TODO: Verify this produces the correct Q^-1 matrix
        suffixes = {'eo', 'oe'};
       
        for j = 1:2
            idStr = ['st4MQ', suffixes{j}];
            [stQ{i}.([idStr]).M11,...
                stQ{i}.([idStr]).M12,...
                stQ{i}.([idStr]).M21,...
                stQ{i}.([idStr]).M22...
            ] = blockInvert(stQ{i}.([idStr]).M11,...
                stQ{i}.([idStr]).M12,...
                stQ{i}.([idStr]).M21,...
                stQ{i}.([idStr]).M22);
        end
       
    end
    
    % Get all values for the current m
    [M, n_vec] = combine_oeeo(stQ{i}, 'st4MQ');
    
    % Create array of indices
    [n,s, np,sp] = ndgrid(n_vec,1:2, n_vec,1:2);
    
    % since m-cell storage was
    %
    %  M11(n,np) | M12(n,np)
    %  --------- + ---------
    %  M21(n,np) | M22(n,np)
    %
    % to grab values by columns, want to vary n, then s, then np, then sp
    
    % reshape as vectors
    vec_t = M(:);
    vec_n = n(:); vec_np = np(:);
    vec_s = s(:); vec_sp = sp(:);
    % m = mp is constant
    vec_m = 0*vec_s + (i - 1);
    vec_mp = vec_m; % axisym
    
    % strip analytical zeros due to eo-oe symmetry
    ids = ~isnan(vec_t);
    
    % save to cell as matrix of format
    % s sp n np m mp Qr Qi
    tmp{i} = [vec_s(ids) vec_sp(ids) vec_n(ids) vec_np(ids) vec_m(ids) vec_mp(ids) real(vec_t(ids)) imag(vec_t(ids))];
end

%% Combine all values for positive m's
% We now have an array of the form
% [vec_s vec_sp vec_n vec_np vec_m vec_mp q_real q_imag]
Q = cell2mat(tmp.');

%% Strip remaining analytic zeroes for m=0
% 1. check components of Q for nonzero values
% 2. if both components are nonzero keep it
% TODO: Shouldn't we keep it if *either* component is nonzero?
ids = all(Q(:,[7 8]), 2); % ids = all(Q(:,[7 8]) ~= 0, 2) seems to have a redundant comparison
Q = Q(ids,:);

%% Add values for negative m's if requested
% TODO: Verify symmetry rules hold for Q matrix
% Warning in place until then

if complete
    warning('Export of negative m values relies on an unconfirmed symmetry relationship.');
    
    ids = Q(:,5) ~= 0; % don't duplicate m=0 cases
    
    vec_m = -Q(ids,5);
    vec_mp = vec_m; % axisym
    % duplicate n and s indices (for m!=0)
    vec_s =  Q(ids,1);
    vec_sp = Q(ids,2);
    vec_n =  Q(ids,3);
    vec_np = Q(ids,4);

    % Q_{-m} = Q_{m} if s=sp, -Q_{m} otherwise
    sgn = (-1).^(Q(ids,1) + Q(ids,2));
    vect_r = sgn .* Q(ids,7);
    vect_i = sgn .* Q(ids,8);

    % update Q to include negative m's
    Q2 = [vec_s vec_sp vec_n vec_np vec_m vec_mp vect_r vect_i];
    Q = [Q ; Q2];

end
    
%% Reorder rows by increasing s sp n np m mp (slow to fast)
% i.e. start with Q11 block
% start with n=1, np=1
% vary m=mp=0:n
% then n=1, np=2, ...
% Note: careful with sortrows when Q contains complex numbers,
% it doesn't treat negative m's like we'd want...
% so now taking real part and imag separately...

% Q = sortrows(Q, [6 5 4 3 2 1]);

Q = sortrows(Q, 1:6);
% [s sp n np m mp Qr Qi]
% [vecs vecsp vecn vecnp vecm vecmp vectr vecti];

% write to a file
if(~isempty(filename))
    if strcmp(filename, 'stdout')
        fileID = 1;
    else
        fileID = fopen(filename, 'w');
    end
    
    if inverse
        fprintf(fileID, '%d elements of Q^(-1)-matrix\n', size(Q, 1));
    else
        fprintf(fileID, '%d elements of Q-matrix\n', size(Q, 1));
    end
    fprintf(fileID, '* s  1st block index (electric/magnetic)');
    fprintf(fileID, '* sp 2nd block index (electric/magnetic)');
    fprintf(fileID, '* n  1st n-index');
    fprintf(fileID, '* np 2nd n-index');
    fprintf(fileID, '* m  1st m-index');
    fprintf(fileID, '* mp 2nd m-index');
    fprintf(fileID, '* Qr real (Q_s,sp,m,mp,n,np)');
    fprintf(fileID, '* Qi imag (Q_ssp,m,mp,n,np)');
    fprintf(fileID, format, Q.');
    fprintf(fileID, 's sp n np m mp Qr Qi \n');

    if ~strcmp(filename, 'stdout')
        fclose(fileID);
    end
end

end

function [E, F, G, H] = blockInvert(A, B, C, D)
    %% Invert 2x2 block matrix
    % input matrix given as blocks arranged as so
    %  A | B   inversion     E | F
    % ---+--- ----------->  ---+---
    %  C | D                 G | H
    %
    % Dependency:
    % invertLUcol
    
    % TODO: Cache or memoize repeated operations?
    
    Ai = invertLUcol(A); % we assume A is nonsingular
    
    H = invertLUcol(D - C*Ai*B); % inverted Schur complement
    E = Ai + Ai * B * H * C * Ai;
    F = - Ai * B * H;
    G = - H * C * Ai;
end