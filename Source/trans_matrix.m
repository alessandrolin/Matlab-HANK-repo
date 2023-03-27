function [L,indeces,values] = trans_matrix(A,grid,P,recover)%#codegen
% Create sizes;
na = size(A,1);
nz = size(A,2);




% % Changes 1 - Sparse right away
% tic
% L2 = sparse(size(grid,1)*size(P,1),size(grid,1)*size(P,1));
% 
% % 
% 
% for ii = 1:size(P,1)
%     aux1 = (A(:,ii)>=grid').*(A(:,ii)<=([grid(2:end);grid(end)])').*(1-(A(:,ii)-grid')./([grid(2:end); grid(end)]'-grid')); 
%     aux2 = [aux1(:,1:end-1) zeros(size(grid,1),1)]+ [zeros(size(grid,1),1) (1-aux1(:,1:end-1)).*(aux1(:,1:end-1)>0) ];
%     aux2(:,end) = (A(:,ii)>=grid(end-1)).*(1-aux2(:,end-1));
%     aux2 = aux2./sum(aux2,2);
%     L2((ii-1)*size(grid,1)+1:ii*size(grid,1),:) = kron(P(ii,:),aux2);
% end
% toc



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



% 
% % Changes 3 - completely different assignment;
% tic
% L3 = sparse(size(grid,1)*size(P,1),size(grid,1)*size(P,1));
% 
% % 1 - Find weights 
% tic
% ind_down = ones(na,nz);
% for ia = 1:na
%     ind_down(A>=grid(ia))=ia;
% end
% 
% wabove_a = (A-grid(ind_down))./(grid(ind_down+1)- grid(ind_down));
% wbelow_a= 1-wabove_a;
% 
% 
% 
% %  tic
% % for ii = 1:size(P,1)
% %     aux1 = (A(:,ii)>=grid').*(A(:,ii)<=([grid(2:end);grid(end)])').*(1-(A(:,ii)-grid')./([grid(2:end); grid(end)]'-grid')); 
% %     aux2 = [aux1(:,1:end-1) zeros(size(grid,1),1)]+ [zeros(size(grid,1),1) (1-aux1(:,1:end-1)).*(aux1(:,1:end-1)>0) ];
% %     aux2(:,end) = (A(:,ii)>=grid(end-1)).*(1-aux2(:,end-1));
% %     aux2 = aux2./sum(aux2,2);
% % end
% % toc
% 
% % 2 - Assign index mapping, values
% step_up_z = kron(0:nz-1,na);
% sparser = cell(na*nz,1); % cell object temp used to collect indices and values to be filled in sparse mat
% for iz = 1:nz
%     for ia = 1:na
%     % Find corresponding state row
%     row     = (ia + (iz-1)*na)*ones(nz*2,1);
%     column  = [ind_down(ia,iz),ind_down(ia,iz)+1]'+step_up_z; 
%     column  = column(:);
%     values  = kron(P(iz,:),[wbelow_a(ia,iz),wabove_a(ia,iz)])';
%     subsparser   = horzcat(row,column,values);
%     % Concatenate Indeces with what was there before;
%     sparser{row(1)}  = subsparser;
%     
%     b=1;
%         
%         
%         
%     end
% end
% 
% %3 -  Fill
% sparserM = cell2mat(sparser);
% L3 = sparse(sparserM(:,1),sparserM(:,2),sparserM(:,3),na*nz,na*nz); % generate sparse matrix
% toc
% 


