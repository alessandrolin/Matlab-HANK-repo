function Vq = nakeinterp1(X, V, Xq)
% function Vq = nakeinterp1(X,V,Xq)
% this is the control function for mex file of nakeinterp1()
% Inputs:
%     X  - double column vector
%     V  - double column vector
%     Xq - double column vector
% Outputs:
%     Vq - double column vector
%
% author: Yu Li
% adjusted JP 5/2020
% Date  : 20th Mar, 2014
% https://raw.githubusercontent.com/mxgiuliani00/M3O-Multi-Objective-Optimal-Operations/master/lib/nakeinterp1.m

if ( isa(X,'double')~=1 || isa(V,'double')~=1 || isa(Xq,'double')~=1 )
    error('input must be ''double'' type')
end

if ( isvector(X)~=1 || isvector(V)~=1 )
    error('input X Y must be a column vector ')
end

if ( isempty(Xq) == 1)
    error('''Xq'' is undefined');
    Vq = [];
end

if ( isvector(Xq)~=1 && isscalar(Xq)~=1)
    error('input Xq must be a column vector or scalar')
end

if issorted(X,'monotonic')
    % make interpolation
    Vq = nakeinterp1mx(X(:), V(:), Xq);
else
    % sort X in ascending order
    [X_sort, idx_x] = sort(X);
    V_sort = V(idx_x);
    
    % make interpolation
    Vq = nakeinterp1mx(X_sort(:), V_sort(:), Xq);
end
