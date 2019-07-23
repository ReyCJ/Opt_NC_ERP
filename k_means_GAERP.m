% ---------K-means Clustering for ERP data different methods--------------

function [K_ERP]=k_means_GAERP(x,k)


%         [c_idx,cen]=kmeans(x,K); % normal kmeans toolbox
%         [c_idx,cen]=kmeans(x,K,'distance','cityblock');
%         [c_idx,cen]=kmeans(x,K,'distance','cosine');
[c_idx,cen]=kmeans(x,k,'distance','correlation');
%         [c_idx,cen]=kmeans(x,K,'distance','Hamming');

K_ERP=c_idx;

end




