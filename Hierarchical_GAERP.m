%  --------------------Hirerichal Clustering ERP data ------------------
% The ERP data has well structure for clustering time sammples and

function [H_ERP]=Hierarchical_GAERP(x,k)

%  clustering with diferent similarity functions---------------------------

%         cl=clusterdata(x,'linkage','complete','distance','euclidean','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','seuclidean','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','cityblock','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','minkowski','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','chebychev','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','mahalanobis','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','cosine','maxclust',k);
cl=clusterdata(x,'linkage','complete','distance','correlation','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','spearman','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','hamming','maxclust',k);
%         cl=clusterdata(x,'linkage','complete','distance','jaccard','maxclust',k);
H_ERP=cl;

end

%----------------- The end of Hierarchical clustering ---------------------