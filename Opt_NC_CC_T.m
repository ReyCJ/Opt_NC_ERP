% Optimal number of clusters determination for ERP data

% by Reza Mahini 2018 
% 1st update Apr 2019
% 2nd Update Nov 2020
% 3rd Update Jan 2022

% Please reference these studies if you used the method or found it useful:

%[1] Mahini R et al. (2020) Determination of the Time Window of Event-Related Potential
% Using Multiple-Set Consensus Clustering Front Neurosci 14 doi:10.3389/fnins.2020.521595

%[2] Mahini R et al. (2019) Optimal Number of Clusters by Measuring Similarity among
% Topographies for Spatio-temporal ERP Analysis arXiv preprint arXiv:191109415
% ### Under review in Brain Topography journal ### 

% Visit our Github @ https://github.com/remahini

% Contact us : r.mahini@gmail.com, remahini@student.jyu.fi

%------------------------------------------------------------------------

clc;
clear;
close all;

% loading data and channel loc -------------------------------------------

load chanlocs.mat;
load Data_noise_awgn;
chanLoc=chanlocs;

% noise levels : 1=-10dB, 2=-5dB, 3=0dB, 4=5dB, 5=10dB, 6=20dB

Data=squeeze(Data_noise_awgn(:,:,:,:,6));

size(Data) % chan x sam x subj x cond


% The information box for optimal number of clusters anaelysis -------
% Pre-saved"compGroup_CC_test.mat/compGroup_CC_test_new85.mat" are needed

% 100x results are in "compGroup_CC_test_new85.mat" file

load compGroup_CC_test_new85.mat;
% OR
% load compGroup_CC_test; % you need to run "ERP_CC_SIM.m" first. 


% Initializing-------------------------------------------------------------

G=1;
St=2;
Sa=300;
Subj=20;

twStart=[183 231]; % (ms) roughly selected measurment interval by experimenter
twEnd=[278 350]; % (ms)

startEph=-100; % (ms)
endEph=600; % (ms)

Chan=65; % number of electrodes
Comp=2; % No. Components of interest

prompt = {'Min number of clusters:','Max number of clusters:'};
dlgtitle = 'Input';
dims = [1 45];
definput = {'2','10'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
a=str2num(answer{1});
b=str2num(answer{2});

thr=0.95; % the threshold for inner-similarity evaluation
min_count=1;
Max_count=size(compGroup_CC,2); % in general 100; % maximum iterations


% Data Preparing ---------------------------------------------------------

inData=DimPrep(Data,Chan,Sa,St,Subj,G);   % Channel x Sample x Stim x Subject x Group

[ERP_Subj,inDaGA_M1]=Data_Preparing(inData,Subj,St,Sa,G); % Concatinating the ERP dataset (M1)

x=inDaGA_M1;  % Fetching dataset


% Main procedure  -------------------------------------------------------------------

tic

for count=min_count:Max_count
    
    [count]
    
    for com=1:Comp
        
        compGroup_CC_inp=[];
        
        [v,w]=time_conv_ts(Sa,twStart(com),twEnd(com),startEph,endEph);
        
        % Setting interest channel locations for analyzing components
        if com==1
            selChan={'P2','P6','PO4'};
            ch_loc=[49   56    62];
        else
            selChan={'CP2','CPz','Cz'};
            ch_loc=[42    58    65];
        end
        
        
        % Select the component of interest dataset
        
        compGroup_CC_inp=compGroup_CC(count).comp(com); % Warning ______       
        
        
        % Optimal number of clusters detection ----------------------------
        
        [Opt_NC,Opt_TW,corr_opt_NC,corr_all]=Opt_NC_Det_SIM(compGroup_CC_inp,a,b,St,G,com,thr);
        
        
        % Collectin all useful information such as: inner-similarities,
        
        for st=1:St
            
            List_OpNC(com).data(count,st)=Opt_NC(st); % OPTNC for
            List_Opt_TW(com).data(count,st,:)=Opt_TW(st,:);  
            List_corr_opt_NC(com).data(count,st)=corr_opt_NC(st);
            List_allCorr(com,count,st,:)=corr_all(st,:);
        end
        
        
        disp(List_OpNC(com).data); % Observing the selections
        
        
        % Statistical analysis for selected optimal NC --------------------
        
        [SPSS_tab_avg]=ERP_statTable2_SIM(Opt_TW,inData,ch_loc,Subj,St,G);
        [ranova_tbl]=ranova_ERP_100(SPSS_tab_avg);
        
        stim_pvalue(count,com,:)=ranova_tbl.pValue(3);
        chan_pvalue(count,com,:)=ranova_tbl.pValue(5);
        intStCh_pvalue(count,com,:)=ranova_tbl.pValue(7);
        
        % -----------------------------------------------------------------
        
    end
    
end

%% Ranking between NCs options

for k=a:b
    for comp=1:Comp
        for st=1:St
            Rank(comp).stim(st).data(k)=length(find(List_OpNC(comp).data(:,st)==k));
        end
    end
end

%%  Selecting Overal_optimal number of clusters --------------------------

corThr=0.03; % Correlation stability parameter #default = 0.02
Thr=0.947;
maxThr=Thr; % Optimal selecting parameter # #default = 0.95

g=1;
for c=1:comp
    
    for st=1:St
        
        allcorrList=squeeze(List_allCorr(c,min_count:Max_count,st,:));
        corr_avg(2,:)=mean(allcorrList',2)';
        corr_avg(1,:)=[a:1:b];
        fi=0; % not detected
        for i=2:size(corr_avg,2)
            if (corr_avg(2,i)>=Thr) && (abs(corr_avg(2,i-1)-corr_avg(2,i))<corThr)&&...
                    (abs(corr_avg(2,i+1)-corr_avg(2,i))<corThr)
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
            
            %             corThr1=corThr1+0.002; % less strict in stability
            % %
            % %         OPNC(g).stim(st).data=[0,0];
        end
        
        %%%
        
        %         id=OPNC(g).stim(st).data(c,:); % OPNC ID
        %          [mival,~]=min(corr_avg(2,:));
    end
end

toc


%% Plot the results ------------------------------------------------------

% Correlation plot for selected TW/NC  -----------------------------------

component={'N2','P3'};

cond=1;
for c=1:comp
    
    figure('Renderer', 'painters', 'Position', [10 10 550 600])
    for st=1:St
        
        subplot(2,1,st,'align');
        corrplot=squeeze(List_allCorr(c,min_count:Max_count,st,:));
        
        meanCorr_cond(:,cond)=mean(corrplot',2);
        cond=cond+1;
        
        plot(corrplot'); % List_allCorr(com,count,st,:)
        hold on
        P1=plot(mean(corrplot',2),'-o','LineWidth',1.2,'MarkerSize',3);
        
        plot([(OPNC(c).stim(st).data(1)-a)+1,(OPNC(c).stim(st).data(1)-a)+1],[0 ,1.15],'--k');
        plot([0,(b-a)+1],[thr,thr],'--k');
        legend([P1],{'Mean CC'},'FontSize',10,'Location','southeast');
        hold off
        xlabel('Number of clusters #');
        ylabel('Inner-similarity');
        
        set(gca,'fontsize',11);
        title([component{c}, ' Component', ', Cond ', int2str(st)]);
        xticks(a-1:1:b-1);
        xticklabels(a:1:b);
        
    end
end

%% Optimal number of clusters

corr_avgmean(1,:)=[a:1:b];
corr_avgmean(2,:)=mean(meanCorr_cond,2);
fi=0; % not detected
for i=2:size(corr_avgmean,2)
    if (corr_avgmean(2,i)>=Thr) && (abs(corr_avgmean(2,i-1)-corr_avgmean(2,i))<corThr)&&...
            (abs(corr_avgmean(2,i+1)-corr_avgmean(2,i))<corThr)
        OPNC_all=[corr_avgmean(1,i),corr_avgmean(2,i)];
        fi=1;
        break
    end
end

corThr1=corThr;  % We need to be flaxable for finding second choice
while fi==0 && corThr1<0.1 % maximum tolorence
    
    for i=2:size(corr_avg,2)
        if (corr_avgmean(2,i)>=maxThr) && (abs(corr_avgmean(2,i-1)-corr_avgmean(2,i))<corThr1)
            OPNC_all=[corr_avgmean(1,i),corr_avgmean(2,i)];
            fi=1;
            break
        end
    end
    
    maxThr=maxThr-0.0005; % decreasing the maximum threshold
    corThr1=corThr1+0.00001;
    
end

figure('Renderer', 'painters', 'Position', [10 10 550 290])

plot(meanCorr_cond); % List_allCorr(com,count,st,:)
hold on
plot(mean(meanCorr_cond,2),'-o','LineWidth',1.2,'MarkerSize',3,'Color','#D95319');


plot([OPNC_all(1)-1,OPNC_all(1)-1],[0 ,1],'--k');
plot([1,(b-a)+1],[thr,thr],'--k');


xlabel('Number of clusters #');
ylabel('Inner-similarity');
set(gca,'fontsize',11);
title('Mean of all conditions');
%         xticks(1:1:14);
%         xticklabels(2:1:15);

tickValues =1:(b-a+1) ;
set(gca,'XTick',tickValues);
label_value=a:b;
set(gca,'XtickLabels',label_value);
ylim([0 1.15])
legend({'N2-Cond1','N2-Cond2','P3-Cond1','P3-Cond1','Mean-all-conds'},'FontSize',10,'Location','southeast')



%% The accuracy of selected TWs compare with original TW (pecentage cover)


GT=[201,261,201,261;266,357,266,364];

for com=1:comp
    
    for st=1:St
        for k=a:b
            for count=1:Max_count
                TW(com,st,k,count,1:2)=compGroup_CC(count).comp(com).sel_TW_ms(k).data(st,2:3);
                TW_diff(com,st,k,count,1)=abs(TW(com,st,k,count,1)-GT(com,2*st-1));
                TW_diff(com,st,k,count,2)=abs(TW(com,st,k,count,2)-GT(com,2*st));
                TW_dur(com,st,k,count)= abs(abs(TW(com,st,k,count,1)-TW(com,st,k,count,2))...
                    -abs(GT(com,2*st-1)-GT(com,2*st)));
            end
        end
    end
end


size(TW_diff)     % 2     2    15   100     2

% average TW

for com=1:comp
    for st=1:St
        for k=a:b
            %           squeeze(TW(com,st,k,:,1:2))
            AvgTW(com,st,k,1:2)= mean(squeeze(TW(com,st,k,:,1:2)));
            SDTW(com,st,k,1:2)= std(squeeze(TW(com,st,k,:,1:2)));
            TW_Duration(com,st,k)=mean(abs(TW(com,st,k,:,1)-TW(com,st,k,:,2)));
        end
    end
end


%% Correlation of OpNC ----------------------------------------------------

for c=1:comp
    
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

disp('st_pv_N2  st_pv_P3  ch_pv_N2  ch_pv_P3  intStCh_pv_N2  intStCh_pv_P3')
all_pvalue=[stim_pvalue,chan_pvalue,intStCh_pvalue];
disp(all_pvalue)


%% Overlap test ------------------------------------------------------------

% % for count=min_count:Max_count
% %     for com=1:2
% %         [v,w]=time_conv_ts(Sa,start_ms,end_ms,twStart(com),twEnd(com));
% %         v=int32(v); w=int32(w);
% %         %     for g=1:G
% %         for st=1:St
% %
% %             temp=squeeze(List_Opt_TW(com).data(count,st,:));
% %
% %             TW=temp(2:3); % be carefull
% %
% %             if (TW(1)-v<0 && TW(2)-v<0)||(TW(1)-w>0 && TW(2)-w>0) % outside condition
% %                 List_Overlap(com).data(count,st)=0;
% %             else
% %                 List_Overlap(com).data(count,st)=1;
% %             end
% %         end
% %     end
% % end

%% Statistical power analysis for within factors (Stim x Chann) -----------

function [ranova_tbl]=ranova_ERP_100(SPSS_tab_avg)

Rdata=SPSS_tab_avg(:,2:7);

Group={'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1'};

varNames={'Group','St1_ch1','St1_ch2','St1_ch3','St2_ch1','St2_ch2','St2_ch3'};

tbl = table(Group,Rdata(:,1),Rdata(:,2),Rdata(:,3),Rdata(:,4),...
    Rdata(:,5),Rdata(:,6),'VariableNames',varNames);

factNames = {'Stim','Chan'};

within_R = table({'St1';'St1';'St1';'St2';'St2';'St2'},{'Ch1';'Ch2';'Ch3';'Ch1';'Ch2';'Ch3'},'VariableNames',factNames);

rm = fitrm(tbl,'St1_ch1-St2_ch3~1','WithinDesign',within_R);

[ranova_tbl] = ranova(rm, 'WithinModel','Stim*Chan');
end


% The END -----------------------------------------------------------------
