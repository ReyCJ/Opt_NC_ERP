% ------------ Usinf diffiusion maps in spectral clustering ---------------
%   References:
%   Ronald R. Coifman and Stéphane Lafon, "Diffusion maps" 2006
%   implemented by: Reza Mahini 2016

function [D_ERP]=DSC_GAERP(x,k)

samples=size(x,1);

% ------------------- Similarity matrix for our data --------------------

W=zeros(length(x),length(x)); % weight matrix making



W1 = squareform(pdist(x));
a=max(max(W1))^2;
b=min(min(W1))^2;
sigma_test=(a- b)/(2*log(a/(b+eps)));  % sigma_test = sigma^2

[epsilon Li]=epsilon_Est(x,'euclidean',samples,'figures');
%----------------------- Making the Gaussian kernel ----------------------
epsilon_mtx=epsilon;

for i=1:length(x)
    
    for j=i+1:length(x)
        dist = norm(x(i,:) - x(j,:))^2;
        W(i,j) = exp(-dist/epsilon);
        W(j,i)=W(i,j);
    end
end

% ------------------------ Making Dgreee matrix ---------------------------

D = zeros(length(x),length(x));
for i=1:length(W)
    D(i,i)=sum(W(i,:));  % degree matrix
end

P=(D^-1)*W; % normalized weight matrix (W)

D(D == 0) = eps; % to prevent divided by zero
D=1./sqrt(D);

P1=D * P * D;

L=D * W * D;

%**************************************************************************

[U, Si, Va]=svd(L); % Singular value decomposition outputs by suing part...
% of "Si" we can decrease dimention of our data

P2=U*Si*Va;
V1=D * U;

psi=V1*Si; % this is the maps

K_lg_eig1=psi(:,1:3); % selecting d=3 dimention space which d<<p


%[V, E] = eig(L); % E =eigenvalue matrix V=eignvector matrix
%[V, E] = eigs(L,10,sigma); % this instruction just make 6 largest eigenvalue and eigenvectors


%---------------------- normalize k largest eigenvectors -----------------

%        K_lg_eig2=normr(K_lg_eig1); % normalize rows

for ii=1:length(K_lg_eig1)
    K_lg_eig(ii,:)=K_lg_eig1(ii,:)./(sqrt(sum(power(K_lg_eig1(ii,:),2))));
    
end
N_space=K_lg_eig();
% ---------------------- clustering based on k-means ---------------------

[c_idx,cen]=kmeans(K_lg_eig,k);
D_ERP=c_idx;

end




