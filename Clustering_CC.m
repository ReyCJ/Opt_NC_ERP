% Clustering with multiple clustering methods ---------------------------


% By: Reza Mahini April 2018, Dalian University of Technology-ASAP_Lab,
% Email : r_mahini@mail.dlut.edu.cn, r.mahini@foxmail.com
% You welcome to use the toolbox and share your idea about it.

clc;
clear;
close all;

% loading data and channel loc -------------------------------------------

load chanlocs.mat;
load Data_LN_N.mat;

Data=squeeze(Data_noise(:,:,:,:,1)); % ERP data with low level noise
size(Data)



% Initializing-------------------------------------------------------------

G=1;
St=2;
Sa=150;
Subj=20;
start_ms=-100; % start of ephochs
end_ms=600; % ephochs end
twStart=[183 231]; % roughly selected measurment interval by experimenter (here we defined by software)
twEnd=[278 350];
Chan=65; % number of channels
Comp=2; % No. Components of interest

a=2; % NC_from
b=15; % NC_to

maxCount=5;


stimSet={'St1','St2'};
compSet={'N200','P300'};

% Data Preparing ---------------------------------------------------------

inData=DimPrep(Data,Chan,Sa,St,Subj,G);   % Channel x Sample x Stim x Subject x Group

[ERP_Subj,inDaGA_M1]=Data_Preparing(inData,Subj,St,Sa,G); % Concatinating the ERP dataset (M1)

x=inDaGA_M1;  % Fetching dataset



%% CC for GA ERP data ----------------------------------------------------

inData=DimPrep(Data,Chan,Sa,St,Subj,G);
% Channel x Sample x Stim x Subject x Group
[GA_ERP]=grand_ERP(inData,G);

[ERP_Subj,inDaGA_M1]=Data_Preparing(inData,Subj,St,Sa,G); % subjects data and grand average data modeled data

x=inDaGA_M1;  % Fetching dataset

tic

for count=1:maxCount
    
    [count]

%% ------------------------------------------------------------------------
% *********************   Clustering with 5 methods    *******************

for k=a:b
    disp('current number of clusters = ')
    [k];
    K_ERP(:,k)=k_means_GAERP(x,k);
    H_ERP(:,k)=Hierarchical_GAERP(x,k);
    F_ERP(:,k)=FCM_GAERP(x,k);
    S_ERP(:,k)=SOM_GAERP(x,k);
    D_ERP(:,k)=DSC_GAERP(x,k);
    
    
    inLabel(k).data=Label_comb_Sim(K_ERP(:,k),H_ERP(:,k),F_ERP(:,k),S_ERP(:,k),D_ERP(:,k));
end


    
    for k=a:b % CC from 2 to 15 % Set OPNC **************

        [k];
        f_result_K=[];
        Cl_pow_res=[];
        Rank_Res=[];
        STD_value=[];
        cluster_N=[];
        
        
        
        methods=[];
      
        
        methods=inLabel(k).data;
        
        
        [p, q]=size(methods); % all method results in one single dataset...
        
        
        %% Preparing HyperGraph (binary) ------------------------------------------
        
        
        H=zeros(p,k*5);  % p=300 number of samples in restructured data
        for i=1:size(methods,2) % =5 number of methods
            for l=1:k % for all clusters
                %temp=[]; no need this
                temp=find(methods(:,i)==l);
                for j=1:size(temp,1)
                    H(temp(j),k*i+l-k)=1; % Hypergraph creation
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
  
        
        clust_idx=clusterdata(L,'linkage','complete','distance','minkowski','maxclust',k);
        CC_idx(:,k)=clust_idx;
        CC_idx_1(:,k,count)=clust_idx;
  %%      
  
        for com=1:Comp
            
            [v,w]=time_conv_ts(Sa,start_ms,end_ms,twStart(com),twEnd(com));
            
            if com==1
                selChan={'P2','P6','PO4'};
                ch_loc=[49   56    62];
            else
                selChan={'CP2','CPz','Cz'};
                ch_loc=[42    58    65];
            end

  
            
            % ------------------ Auto Area Detection and results ------------------
            
            % x1=squeeze(x(2400-2399:g*2400,:));
            
            Clu_idx=clust_idx;
            
            g=1;
            
            [CSPA_f_result,comp_pow,SelStd_val,innerCorr,winnID]=Comp_detect_ERP_CC_100s(Clu_idx,x,chanlocs,k,Sa,St,v,w,g,com);
            
            
            [selected_TW,TWs_ms,selTWs_ms,sel_STD,sel_innerCorr.gr]=...
                Sel_TW_100s(CSPA_f_result,SelStd_val,innerCorr, v,w,St,g,Sa,start_ms,end_ms); % TWs selection algorithm

            % --------------------------------------------------------------------
            
            selTWs(com).data=selTWs_ms;
            
            disp(compSet{com})
            disp(selTWs(com).data)
            % %     %% ------------------ Auto Area Detection and results ------------------
            
            compGroup_CC(count).comp(com).std(k).data=sel_STD; %  Sel_STD(st,g) , st =stimulus, g= group
            compGroup_CC(count).comp(com).Corr(k).data=sel_innerCorr; % innerCorr
            compGroup_CC(count).comp(com).sel_TW(k).data=selected_TW;
            compGroup_CC(count).comp(com).idx=CC_idx;

        end
        
    end
end

save compGroup_CC.mat

toc
