% This procedure could find the optimum number of clusters using Consensus
% Clustering

% This version is the Demo version with the alreadey prepared results by
% consensus clustering on simulated ERP data ("Data_LN_N")

% By: Reza Mahini 2018 updated Apr 2019
% 
% Respect to Academic rules please cite the toolbox if you find it useful
% and use it in your work!
% Here is how to cite:
% 
% You are welcome to use the toolbox, we will be
% glad to hear your oponion about the toolbox.

% Cheers, 

% Reza Mahini

% Email : r.mahini@gmail.com, Dalian University of Technology, ASAP Group


clc;
clear;
close all;

% loading data and channel loc -------------------------------------------

load chanlocs.mat;
load Data_LN_N.mat;

% The dataset for optimal number of clusters anaelysis 

% % load compGroup_CC; % loading all results from CC (i.e. consensus clustering results...
% from 5 clustering methods)

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

Thr=0.95; % the threshold for inner-similarity evaluation #default =0.95
corThr=0.02; % Correlation stability parameter #default = 0.02

min_count=1;
Max_count=input('Number of iterative running the method (i.e. must bigger than 1, recommended 20 for Simulated ERP) ='); % maximum iterations




% Data Preparing ---------------------------------------------------------

inData=DimPrep(Data,Chan,Sa,St,Subj,G);   % Channel x Sample x Stim x Subject x Group

[ERP_Subj,inDaGA_M1]=Data_Preparing(inData,Subj,St,Sa,G); % Concatinating the ERP dataset (M1)

x=inDaGA_M1;  % Fetching dataset

compGroup_CC= CC_ERP_GA_Sim(x,a,b,chanlocs,Sa,St,Comp,Max_count,twStart,twEnd,start_ms,end_ms);




%%  Main procedure  -------------------------------------------------------------------

tic

for count=min_count:Max_count
    
    [count]
    
    for com=1:Comp
        
        compGroup_CC_inp=[];
        
        [v,w]=time_conv_ts(Sa,start_ms,end_ms,twStart(com),twEnd(com));
        
        % Setting interest channel locations for analyzing components
        if com==1
            selChan={'P2','P6','PO4'}; 
            ch_loc=[49   56    62];
        else
            selChan={'CP2','CPz','Cz'};
            ch_loc=[42    58    65];
        end
        
        
        % Select the component of interest dataset 
        
        compGroup_CC_inp=compGroup_CC(count).comp(com);
        
        
       
        % Optimal number of clusters detection ----------------------------
        
        [Opt_NC,Opt_TW,corr_opt_NC,std_Opt_NC,corr_all]=Opt_NC_Det_Sim(compGroup_CC_inp,a,b,St,G,com,Thr);
        
        % end
        
        
        % Collectin all useful information such as: inner-similarities, 
        
        for st=1:St
            
            List_OpNC(com).data(count,st)=Opt_NC(st); % OPTNC for
            %        List_OpNC(com).data(count,2)=Opt_NC.stim(2).data(1);
            List_Opt_TW(com).data(count,st,:)=Opt_TW(st,:);   % Opt_TW.stim(st).data(2:3);
            List_corr_opt_NC(com).data(count,st)=corr_opt_NC(st);
            List_std_Opt_NC(com).data(count,st)=std_Opt_NC(st);
            List_allCorr(com,count,st,:)=corr_all(st,:);
        end
        
        
        disp(List_OpNC(com).data); % Observing the selections
        
        
        % Statistical analysis for selected optimal NC --------------------
  
        [SPSS_tab_avg]=ERP_statTable2_100s(Opt_TW,inData,ch_loc,Subj,St,G);
        [ranova_tbl]=ranova_ERP_Sim(SPSS_tab_avg);
        
        stim_pvalue(count,com,:)=ranova_tbl.pValue(3);
        chan_pvalue(count,com,:)=ranova_tbl.pValue(5);
        intStCh_pvalue(count,com,:)=ranova_tbl.pValue(7);
        
        % -----------------------------------------------------------------
        
    end
    
end


%%  Selecting Overal_optimal number of clusters --------------------------

maxThr=Thr; % Optimal selecting parameter # #default = 0.95

g=1;
for c=1:Comp
    
    for st=1:St
        
        allcorrList=squeeze(List_allCorr(c,min_count:Max_count,st,:));
        corr_avg(2,:)=mean(allcorrList',2)';
        corr_avg(1,:)=[2:1:15];
        fi=0; % not detected
        for i=2:size(corr_avg,2)
            if (corr_avg(2,i)>=Thr) && (abs(corr_avg(2,i-1)-corr_avg(2,i))<corThr)
                OPNC(c).stim(st).data=[corr_avg(1,i),corr_avg(2,i)];
                fi=1;
                break
            end
        end
        
        corThr1=corThr;  % We need to be flaxable for finding second choice
        while fi==0 && corThr1<0.1 % maximum tolorence
            
            for i=2:size(corr_avg,2)
                if (corr_avg(2,i)>=maxThr) && (abs(corr_avg(2,i-1)-corr_avg(2,i))<corThr1)
                    OPNC(c).stim(st).data=[corr_avg(1,i),corr_avg(2,i)];
                    fi=1;
                    break
                end
            end
            
            maxThr=maxThr-0.0005; % decreasing the maximum threshold
            corThr1=corThr1+0.00001;
            
          
        end
        
      
    end
end

toc


%% Plot the results ------------------------------------------------------

% Correlation plot for selected TW/NC  -----------------------------------

component={'N2','P3'};


for c=1:Comp
    
    figure('Renderer', 'painters', 'Position', [10 10 600 600])
    for st=1:St
        
        subplot(2,1,st,'align');
        corrplot=squeeze(List_allCorr(c,min_count:Max_count,st,:));
        plot(corrplot'); % List_allCorr(com,count,st,:)
        hold on
        plot(mean(corrplot',2),'-o','LineWidth',2,'MarkerSize',5);
        
        plot([(OPNC(c).stim(st).data(1)-a)+1,(OPNC(c).stim(st).data(1)-a)+1],[0 ,1.1],'--k');
        plot([0,(b-a)+1],[Thr,Thr],'--r');
        
        hold off
        xlabel('Number of clusters #');
        ylabel('Correlation Coefficient');
        set(gca,'fontsize',12);
        title(['Correlation Coefficient for Sel.TW, ' ,component{c}, ', Stim ', int2str(st)]);
        xticks(1:1:14);
        xticklabels(2:1:15);
        
    end
end

%% Overlap test ------------------------------------------------------------

for count=min_count:Max_count
    for com=1:Comp
        [v,w]=time_conv_ts(Sa,start_ms,end_ms,twStart(com),twEnd(com));
        v=int32(v); w=int32(w);
        %     for g=1:G
        for st=1:St
            
            temp=squeeze(List_Opt_TW(com).data(count,st,:));
            
            TW=temp(2:3); % be carefull
            
            if (TW(1)-v<0 && TW(2)-v<0)||(TW(1)-w>0 && TW(2)-w>0) % outside condition
                List_Overlap(com).data(count,st)=0;
            else
                List_Overlap(com).data(count,st)=1;
            end
        end
    end
end

%% Correlation of OpNC ----------------------------------------------------

for c=1:Comp
    
    figure('Renderer', 'painters', 'Position', [10 10 700 600])
    for st=1:St
        
        subplot(2,1,st,'align');
        bar(List_corr_opt_NC(c).data(:,st));
        xlabel('Run no. #');
        ylabel('Correlation coefficient');
        set(gca,'fontsize',12);
        title(['Correlation coefficient Sel.TW ' ,component{c}, ' and Cond ', int2str(st)]);
        %         ylim([0 1.1])
    end
end


%% Plot for P-values from different factors -------------------------------

all_pvalue=[stim_pvalue,chan_pvalue,intStCh_pvalue];
sataSet={'Stimulus', 'Electrode', 'Stimulus x Electrode'};
component={'N2','P3'};

cu=[1 3 5;2 4 6];

for c=1:Comp
    figure('Renderer', 'painters', 'Position', [10 10 600 700])
    
    for fact=1:3
        
        subplot(3,1,fact,'align');
        bar(all_pvalue(:,cu(c,fact)));
        xlabel('Run no. #');
        ylabel('p-value #');
        set(gca,'fontsize',11);
        title(['Statistical analysis for ' , component{c}, ', ' ,sataSet{fact}]);
        %     ylim([0 .06]);
        
    end
    
end


%% Box plot --------------------------------------------------------------

figure('Renderer', 'painters', 'Position', [10 10 700 400]);

var2={'N2_Stim','P3_Stim','N2_Chan','P3_Chan','N2_Stim*Chan','P3_Stim*Chan'};
boxplot([stim_pvalue,chan_pvalue,intStCh_pvalue],var2);
%ylim([-0.01 0.1]);
% xlabel('Condition #');
ylabel(' p-value ratio %#');
set(gca,'fontsize',12);
title(['p-value for 100 times running procedure']);

% The END -----------------------------------------------------------------
