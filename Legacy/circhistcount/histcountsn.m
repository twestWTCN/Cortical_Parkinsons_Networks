%HISTNC n-dimensional histogram count.
%   N = HISTCOUNTSN(X,EDGES), for a set of vectors X, counts the vectors that
%   fall into each of the n-dimensional histogram bins specified by the
%   EDGES cell array. X is given as m x n matrix, with n the dimension and
%   m the number of vectors. EDGES is a cell array with n entries, each
%   cell specifies the edges of the bins along one dimension (which must be
%   monotonically non-decreasing values). If the histogram is
%   one-dimensional, EDGES can also be a vector. N a matrix of size
%   numel(EDGES{1}) x numel(EDGES{2}) x ... x numel(EDGES{n}) containing
%   the counts for each bin.
%
%   N will count a vector X(i,:) if EDGES{d}(k) <= X(i,d) < EDGES{d}(k+1)
%   for all d in [1,n]. The last bin in each dimension will count any
%   values of X(:,d) in EDGES{d}(end-1) <= X(i,d) <= EDGES{d}(end).
%   Values outside the values in EDGES are not counted.  Use -inf and inf
%   in EDGES to include all non-NaN values.
%
%   N = HISTCOUNTSN(X,EDGES,ISCIRCULAR) lets you specify circular bins
%   in a dimension. ISCIRCULAR is a 1 x n vector with false for
%   non-circular dimensions (default) and true for circular ones. The edges
%   have to be monotonically increasing in [0,2*pi). X may take on any
%   value in that dimension, but the values are assumed to be given in
%   radians.
%
%   N = HISTCOUNTSN(X,EDGES,ISCIRCULAR,RETURNMAT) return a histogram matrix
%   of size numel(EDGES{1}) x numel(EDGES{2}) x ... x numel(EDGES{n})
%   instead of a column vector if RETURNMAT is true (default) or all bins
%   in a single column vector if RETURNMAT is false. The second options
%   should be slightly faster.
%
%   [N,BIN] = HISTCOUNTSN(X,EDGES,...) also returns an index matrix BIN. If X is
%   a vector, N(K) = SUM(BIN==K). BIN is NaN for out of range values (note
%   that this is different from histc, where out of range values are 0).
%
%   [N,BIN,SUB] = HISTCOUNTSN(X,EDGES,...) also returns an subscript matrix SUB.
%   SUB is a m x n matrix containing the bin subscript for each vector in
%   X. SUB is NaN for out of range values.
%
%   See also HISTC,HISTCOUNTS,DISCRETIZE.

function [n,bin,sub] = histcountsn(x,edges,iscircular,returnmat)

    if ~iscell(edges)
        edges = {edges};
    end

    if size(x,2) ~= numel(edges)
        error('Must give one bin range for each dimension.');
    end
    
    if nargin < 3 || isempty(iscircular)
        iscircular = false(1,size(x,2));
    elseif numel(iscircular) ~= size(x,2)
        error('Option iscircular must contain one value for each dimension.');
    end
    
    if nargin < 4 || isempty(returnmat)
        returnmat = true;
    elseif ~isscalar(returnmat)
        error('Option returnmat must be a scalar.');
    end
    
    % get n-dimensional subscripts of the bins
    sub = cell(1,size(x,2));
    nbins = zeros(1,size(x,2));
    for i=1:size(x,2)
        
%         % version >= r2015a (maybe faster, but untested)
%         if iscircular(i)
%             % shift all values and all bin ranges so that the first bin
%             % range start is at 0
%             sub{i} = discretize(...
%                 mod(x(:,i)-edges{i}(1),2*pi),...
%                 [edges{i}-edges{i}(1),2*pi]); % newer (but untested)
%             nbins(i) = numel(edges{i});
%         else
%             sub{i} = discretize(x(:,i),edges{i}); % newer (but untested)
%             nbins(i) = numel(edges{i})-1;
%         end
        
        % version < r2015a (but still supported for the forseeable future)
        if iscircular(i)
            % shift all values and all bin ranges so that the first bin
            % range start is at 0
            [~,sub{i}] = histc(...
                mod(x(:,i)-edges{i}(1),2*pi),...
                [edges{i}-edges{i}(1),2*pi]);
            nbins(i) = numel(edges{i});
        else
            [~,sub{i}] = histc(x(:,i),edges{i});
            sub{i}(sub{i}==numel(edges{i})) = numel(edges{i})-1; % last bin also has values == edges{i}(end)
            nbins(i) = numel(edges{i})-1;
        end
        sub{i}(sub{i}==0) = nan; % replace 0's (not in a bin) with nan
        
    end
    
    % get linear indices of the bins
    bin = sub2ind(nbins,sub{:});
    
    % accumulate count into linear bins
    n = accumarray(bin(~isnan(bin)),1,[prod(nbins),1]);
    
    if nargout >= 3
        sub = [sub{:}];
    end
    
    if returnmat
        n = reshape(n,nbins);
        if nargout >= 2
            bin = reshape(bin,nbins);
        end
    end
end
