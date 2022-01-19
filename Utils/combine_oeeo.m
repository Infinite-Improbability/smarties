function [M, nvec] = combine_oeeo(stMa, fieldname)
    %% Combine oe and eo matrices
    % Returns full matrix from a struct stored in block-rvh form
    % modified from rvhGetFullMatrix to have NaNs
    % 
    %
    % Input:
    %           - st4Ma: struct of matrix in block-rvh form
    %                   with a CsMatList field and at least two fields (ending in "eo" and "oe")
    %                   each is a struct with fields M11,M12,M21,M22,m,ind1,ind2
    % Output:
    %           - M: the full square matrix of size [2Nm x 2Nm] where Nm=N+1-m
    %                (or N if m=0)
    %                Note that Nm=length(ind1)+length(ind2)
    %           - nvec: [Nm x 1] the n (or k) - values each block of the matrix
    %                   corresponds to
    %
    % Dependency:
    % none
    
    if(nargin < 2)
        fieldname = 'st4MT';
    end
    
    
    % Get oe struct
    % interpolation of field name
    st4M = stMa.([fieldname 'oe']);
    % st4M = stMa.st4MToe;
    
    ind1=st4M.ind1;
    ind2=st4M.ind2;
    N1=length(ind1);
    N2=length(ind2);
    
    Nm=N1+N2;
    M=NaN(2*Nm); % fill with NAs to keep track of analytical zeros
    m=st4M.m;
    
    M(ind1,ind1) = st4M.M11;
    M(ind1,Nm+ind2) = st4M.M12;
    M(Nm+ind2,ind1) = st4M.M21;
    M(Nm+ind2,Nm+ind2) = st4M.M22;
    nvec = ( (max(m,1)): (Nm + max(m,1) - 1)) .';
    
    % Get eo struct
    st4M = stMa.([fieldname 'eo']);
    % st4M = stMa.st4MTeo;
    
    % Complete the matrix
    % (note that in principle, ind1 is same as ind2 before and vice versa)
    ind1=st4M.ind1;
    ind2=st4M.ind2;
    N1=length(ind1);
    N2=length(ind2);
    Nm=N1+N2;
    
    M(ind1,ind1) = st4M.M11;
    M(ind1,Nm+ind2) = st4M.M12;
    M(Nm+ind2,ind1) = st4M.M21;
    M(Nm+ind2,Nm+ind2) = st4M.M22;
    
    end