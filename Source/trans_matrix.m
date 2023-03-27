function [L,indeces,values] = trans_matrix(A,grid,P,recover)
% Create sizes;
na = size(A,1);
nz = size(A,2);


L = sparse(size(grid,1)*size(P,1),size(grid,1)*size(P,1));
% 

% 1 - Find weights 
ind_down = ones(na,nz);
for ia = 1:na
    ind_down(A>=grid(ia))=ia;
end
ind_down = min(ind_down,na-1);
wabove_a = (A-grid(ind_down))./(grid(ind_down+1)- grid(ind_down));
wbelow_a= max(1-wabove_a,0);

for ii = 1:size(P,1)
    aux = sparse([1:na,1:na]',[ind_down(:,ii);ind_down(:,ii)+1],[wbelow_a(:,ii);wabove_a(:,ii)],na,na); % generate sparse matrix
    L((ii-1)*size(grid,1)+1:ii*size(grid,1),:) = kron(P(ii,:),aux);
end


% If we want to recover the matrix of indeces and values:
if exist('recover','var')
    [rows,columns,values] = find(L);
    if nargout>1
        indeces               = [rows,columns];
    end
else
    indeces = nan; values = nan;
end    


