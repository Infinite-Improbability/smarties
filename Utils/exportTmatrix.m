function [T] = exportTmatrix( stT, complete, filename, format )
%% exportTmatrix
% Reshaping to long format and exporting T-matrix entries to a text file
%
% PARAMETERS:
% - stT: structure containing T-matrix elements, as returned by the program
% - complete: [logical] compute negative m's
% - out: optional output filename
% - format: string format for text output
%
% RETURNS: T-matrix elements and indices consolidated in a single matrix
% The output format consists of 8 columns:
% s sp m mp n np Tr Ti
% whereby
% * m  1st m-index
% * mp 2nd m-index (identical, due to rotational symmetry)
% * n  1st n-index
% * np 2nd n-index
% * s  1st block index (electric/magnetic)
% * sp 2nd block index (electric/magnetic)
% * Tr real(T_sspmmpnnp)
% * Ti imag(T_sspmmpnnp)
%
% Dependency:
% combine_oeeo

if(nargin < 2)
    complete = true;
end
if(nargin < 3)
    filename = [];
end
if(nargin < 4)
%     format = '%d %d %d %d %d %d % d % d %.15g %.15g\n';
    format = '%d %d %d %d % d % d %.15g %.15g\n';
end

mMax = length(stT);

%% extract T-matrix values for each m and list with corresponding indices
% loop over positive m-values [following storage as independent cells]
tmp = cell(1,mMax); % temporarily store in a cell array
for i_m =  1:mMax

    % get all values for current m
    [M, nvec] = combine_oeeo(stT{i_m});
    % create array of indices
% seems to be incorrect order (transposed)
%     [np,sp, n,s] = ndgrid(nvec,1:2, nvec,1:2);
    [n,s, np,sp] = ndgrid(nvec,1:2, nvec,1:2);
    % since m-cell storage was
    %
    %  M11(n,np) | M12(n,np)
    %  --------- + ---------
    %  M21(n,np) | M22(n,np)
    %
    % to grab values by columns, want to vary n, then s, then np, then sp
    % want to vary
    % reshape as vectors
    vect = M(:);
    vecn = n(:); vecnp = np(:);
    vecs = s(:); vecsp = sp(:);
    % m = mp is constant
    vecm = 0*vecs + (i_m - 1);
    vecmp = vecm; % axisym
    % strip analytical zeros due to eo-oe symmetry
    ids = ~isnan(vect);
    tmp{i_m} = [vecs(ids) vecsp(ids) vecn(ids) vecnp(ids) vecm(ids) vecmp(ids) real(vect(ids)) imag(vect(ids))];

end

% combine all values for positive ms
% now have an array of the form
% vecs vecsp vecm vecmp vecn vecnp treal timag
T = cell2mat(tmp.');

% strip remaining analytic zeroes for m=0
ids = all(T(:,[7 8]) ~= 0, 2);
T = T(ids,:);

%% add values for negative m's if requested
if complete

    ids = T(:,5) ~= 0; % don't duplicate m=0 cases

    vecm = -T(ids,5);
    vecmp = vecm; % axisym
    % duplicate n and s indices (for m!=0)
    vecs =  T(ids,1);
    vecsp = T(ids,2);
    vecn =  T(ids,3);
    vecnp = T(ids,4);

    % T_{-m} = T_{m} if s=sp, -T_{m} otherwise
    sgn = (-1).^(T(ids,1) + T(ids,2));
    vectr = sgn .* T(ids,7);
    vecti = sgn .* T(ids,8);

    % update T to include negative m's
    T2 = [vecs vecsp vecn vecnp vecm vecmp vectr vecti];
    T = [T ; T2];

end


%% reorder rows by increasing s sp n np m mp (slow to fast)
% i.e. start with T11 block
% start with n=1, np=1
% vary m=mp=0:n
% then n=1, np=2, ...
% Note: careful with sortrows when T contains complex numbers,
% it doesn't treat negative m's like we'd want...
% so now taking real part and imag separately...

% T = sortrows(T, [6 5 4 3 2 1]);
% T = sortrows(T, 1:6);

T = sortrows(T, 1:6);
% [vecs vecsp vecn vecnp vecm vecmp vectr vecti];
vecp = p_index(T(:,3), T(:,5)); % n, m -> p 
% vecpp = p_index(T(:,4), T(:,6));
% T = [vecp vecpp T];


% write to a file
if(~isempty(filename))
    if strcmp(filename, 'stdout')
        fileID = 1;
    else
        fileID = fopen(filename, 'w');
    end
    fprintf(fileID, '%d elements of T-matrix\n', size(T, 1));
    fprintf(fileID, 's sp n np m mp Tr Ti \n');
    fprintf(fileID, format, T.');
    if ~strcmp(filename, 'stdout')
        fclose(fileID);
    end
end



end


function [out1] = p_index(in1,in2)
% (n,m) -> p = n * (n+1) + m
out1 = in1 .* (in1 + 1) + in2;
end
