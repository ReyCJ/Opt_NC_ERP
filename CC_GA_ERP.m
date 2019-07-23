% Cluster-based Similarity Partitioning Algorithm (CSPA) for Grand average
% on Stimuli for ERP data analysis.

%%
% * * _%By Reza Mahini May 2017_ * *

function [compGroup_CC]=CC_GA_ERP(inData,inDaGA,chanlocs,Subj,...
    Sa,St,G,start_ms,end_ms,twStart,twEnd,inLabel,a,b,filepath)


[v,w]=time_conv_ts(Sa,start_ms,end_ms,twStart,twEnd);



% -------------------------------------------------------------------------
% [selChan,ch_loc]=selChan_Satis(chanlocs);

tic

for k=a:b  % CC from 2 to 15 by default
    
   fig1=1;
% %    for g=1:G 
    f_result_K=[];
    Cl_pow_res=[];
    Rank_Res=[];
    STD_value=[];
    cluster_N=[];
    
% %     x=ERP_Schiz_GA;  % Fetching dataset
    x=inDaGA;  % Fetching dataset

        
    methods=[];
    
    methods=squeeze(inLabel(k).data); %******* important Sel. Dataset ********* 1050x5
    
    [p, q]=size(methods); % all method results in one single dataset...
    
   
    %% Preparing HyperGraph (binary) ------------------------------------------
  
    H=zeros(p,k*5);  % p=4800
    for i=1:size(methods,2) % =5
        for st=1:k % for all clusters
            %temp=[]; no need this
            temp=find(methods(:,i)==st);
            for j=1:size(temp,1)
                H(temp(j),k*i+st-k)=1; % Hypergraph creation
            end
        end
    end
    
    T=zeros(p); % p=1050
    HE=[]; % hyperedges for each mathod
    

    for i=1:size(methods,2) % number of methods =5
        HE(i,:,:)=H(:,k*i-(k-1):k*i); % separation of HyperEdges
        S(i,:,:)=squeeze(HE(i,:,:))*squeeze(HE(i,:,:))'; % Similarity matrix

        T=T+squeeze(S(i,:,:));

    end
    
    % The Combination of Similarity for all methods --------------------------
    i=i+1;
    L=T/5; % Averaging similarities same as S=(1/r)HH'

    %% Reclustering ------------------------------------------------------------
    
    
    %sel=input('Select method for clustering "K" for k-means and "H" for Hierarchical : ')
    
    clust_idx=clusterdata(L,'linkage','complete','distance','minkowski','maxclust',k);
    Cl_idx(:,k)=clust_idx;
  

%% ------------------- indexing -----------------------------------------
 for g=1:G   
    for i=1:k % number of clusters 
        cl_inf(i).data= find(clust_idx(g*Sa*St-(Sa*St-1):g*Sa*St)==i);
    end
    % index(1,:)=[1 c_idx(1)]; % may need this
    
    clust_idx1=clust_idx(g*Sa*St-(Sa*St-1):g*Sa*St); 
    
    index=[];
    index(1,:)=[1 clust_idx1(1)];
    j1=2;
    for j=1:Sa*St-1 % 2400-1
        if abs(clust_idx1(j)-clust_idx1(j+1))>=1
            index(j1,:)=[j+1 clust_idx1(j+1)]; % Start end of each clusters
            j1=j1+1;
        end
    end
    index(j1,:)=[j+1 clust_idx1(j)];
    
   
 % ------------------------------------------------------------------------

    fac=(end_ms-start_ms)/Sa;
    
        x1=squeeze(x(g*Sa*St-(Sa*St-1):g*Sa*St,:));
    for t=1:St
 %% ------------------ Auto Area Detection and results ------------------
    
    Clu_idx=squeeze(Cl_idx(g*Sa*St-(Sa*St-1):g*Sa*St,k,:));
 
    [CSPA_f_result,CSPA_Cl_pow_res,SelStd_val,innerCorr]=Comp_detect_ERP_CC_100s(Clu_idx,x1,chanlocs,k,Sa,St,v,w,g,filepath);
    [selected_TW(:,:,g),TWs_ms(:,:,g),selTWs_ms(:,:,g),sel_STD(:,g),sel_innerCorr(g).gr]=Sel_TW(CSPA_f_result(:,:,g),SelStd_val(:,:,g),innerCorr, v,w,St,g,Sa,start_ms,end_ms,filepath); % TWs selection algorithm
 
 end
 
    compGroup_CC.std(k).data=sel_STD; %  Sel_STD(st,g) , st =stimulus, g= group
    compGroup_CC.Corr(k).data=sel_innerCorr; % innerCorr 
    compGroup_CC.sel_TW(k).data=selected_TW;
    compGroup_CC.idx(k).data=Clu_idx;
    compGroup_CC.sel_TW_ms(k).data=TWs_ms;
    
    
end
end

