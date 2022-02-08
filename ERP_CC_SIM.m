% Consensus clustering using Cluster-based Similarity Partitioning Algorithm (CSPA) for Grand average
% on Stimuli for grand average ERP data analysis.

% History :(
% First version By Reza Mahini May 2017
% updated at Aug 2019  Stabilization
% 1st Update Apr 2019
% 2nd Update Nov 2020
% 3rd Update Jan 2022

% Visit our Github @ https://github.com/remahini
% Contact us : r.mahini@gmail.com, remahini@student.jyu.fi

% ***** NOTE *****: Running this code provides a long loop (many repeats)
% and it might take long depend on your Data and machine hardware
% Aprox:1h~2d :|

% Please reference these studies if you used the method or found it useful:

%[1] Mahini R et al. (2020) Determination of the Time Window of Event-Related Potential
% Using Multiple-Set Consensus Clustering Front Neurosci 14 doi:10.3389/fnins.2020.521595

%[2] Mahini R et al. (2019) Optimal Number of Clusters by Measuring Similarity among
% Topographies for Spatio-temporal ERP Analysis arXiv preprint arXiv:191109415
% ### Under review in Brain Topography journal ### 

% -------------------------------------------------------------------------


clc;
clear;
close all;
delete K_lg_eig_M2.mat % to make sure the dataset is replaced with new

oldpath = path;
path(oldpath,'ClsMth&Funcs')

% ---------------------------- Input ERP ----------------------------------

load chanlocs;
load Data;

chanLoc=chanlocs;

size(Data) % chan x sam x subj x cond

% ------------------------------------------------------------------------
%  **************************** Initializing ******************************

G=1;
St=2;
Sa=300; % upsampled from 150 to 300
Subj=10; % limitted data due to volume
startEph=-100;
endEph=600;
Comp=2;
Chan=65;

% Processing parameters -------------------------------------------------

Maxcount=input('how many repeat do you need? (def= 100)'); % repetitions
if isempty(Maxcount)
   Maxcount=100;
end

prompt = {'Min number of clusters:','Max number of clusters:'};
dlgtitle = 'Input';
dims = [1 45];
definput = {'2','10'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
a=str2num(answer{1});
b=str2num(answer{2});

minSamThr=15;
InSim_Thr=input('The inner similarity threshold 0 <InSim_Thr <1 (def. = 0.95)?');
if isempty(InSim_Thr)
   InSim_Thr=0.90; % **** Inner-similarity threshold, it lets to include map with ...
end % higher than this threshold to join candidate map list

twStart=[190 256]; % (ms) [194 250];%  [194 240]; % twStart=[63 73]; % time sample(ts)
twEnd=[255 310]; % (ms) [250 350];% [278 385]; %  twEnd=[81 104];

stimSet={'Cond1','Cond2'};
compSet={'N2','P3'};
% [G,k,start_ms,end_ms,twStart,twEnd]=Init_ERP_CC();

% ------------------------------------------------------------------------
% ****************************** Data preparing  **************************

inData=DimPrep(Data,Chan,Sa,St,Subj,G);
% Channel x Sample x Stim x Subject x Group

[GA_ERP]=grand_ERP(inData,G);

[ERP_Subj,inDaGA_M1]=Data_Preparing(inData,Subj,St,Sa,G); % subjects data and grand average data modeled data

x=inDaGA_M1; % the grand average data ""samples x channels x group""

nSam=St*Sa; % special parameter

% Methods for generation step ---------------------------------------------

stb=input('Do we need stablilization? Enter 1=yes, 0=No (def)?');
if isempty(stb)
   stb=0; % **** Inner-similarity threshold, it lets to include map with ...
end % higher than this threshold to join candidate map list

clmethod={'KM','HC','FCM','SOM ','DFS', 'MKMS', 'AAHC', 'SPC', 'KMD', 'GMM'};

% **** Note: 
% We suggest to select group of polarity-variant (KM, FCM, etc.) OR
% polarity-invariant (MKM, AAHC) in once generation phase.

% Method codes: 1=K-means', 2= 'Hierachial', 3= 'FCM', 4 ='SOM ','5 = DFS',...
% 6= 'MKMS', 7= 'AAHC', 8= 'SPC', 9= 'KMD', 10='GMM'

M_list=[1 2 4 8 9]; % configuration of Consenus clustering (Recommended for test [1 2 8 9 10]);
Stb_list=[1 9]; % methods need stabilization
rep=[4 3]; % stabilization repeatation


% % % IndAnlz=input('Do you need individual subject analysis results (No =0, Yes=1) ? ');

plotonoff=input(' Do you need the component plot for each subject (No =0 (def), Yes=1) ?');
if isempty(plotonoff)
   plotonoff=0;
end


% the main procedure -----------------------------------------------------
tic
for count=1:Maxcount

   for k=a:b

      %  ------------------------------------------------------------------
      % ******************** Consensus Clustering ************************

      disp(['Counting  = ',num2str(count)]);
      disp(['current number of clusters = ',num2str(k)]);

      MethLabs_M1(count).ncl(k).labels=Labeling_AllCls_M1(inDaGA_M1,k,rep,stb,M_list,Stb_list);

      methods=[];

      methods=MethLabs_M1(count).ncl(k).labels.data;

      [p, q]=size(methods); % all method results in one single dataset...

      % Consenus finction ------------------------------------------------

      Clu_idx=CSPA(methods,k);


      %%  Time window determination ---------------------------------------

      for com=1:Comp

         [v,w]=time_conv_ts(Sa,twStart(com),twEnd(com),startEph,endEph);

         for g=1:G

            x1=squeeze(x(:,:,g));
            index=Indexing_M1(nSam,Clu_idx,k,St,g);

            CSPA_f_result=[];
            comp_pow=[];
            innerCorr=[];
            winnID=[];
            InnSimWcl=[];

            [CSPA_f_result,comp_pow,innerCorr,winnID,InnSimWcl]=...
               Comp_detect_ERP_CC_Upd(Clu_idx,x1,chanlocs,k,Sa,St,v,w,com,stimSet,compSet,InSim_Thr,minSamThr);

            [selected_TW,TWs_ms,selTWs_ms,sel_innerCorr,InnSim,selPower_amp]=...
               Sel_TW_Upd(CSPA_f_result,innerCorr,v,w,St,g,Sa,startEph,endEph,winnID,InnSimWcl,comp_pow); % TWs selection algorithm


            for st=1:St

               compCorr=[];
               compCorr=sel_innerCorr(st).data;
               n=size(compCorr,1);
               meanRow=sum(sum(compCorr,1))-n;
               inSim(st)=meanRow/(n^2-n);

               meantop_amp(st,:)=selPower_amp(st).data;

            end

            % Plot the component1 --------------------------------------------
            if plotonoff
               %
               if com==1
                  selChan1={'P6','PO4'}; %{'P2','P6','PO4'};
                  ch_loc1=[56    62]; %[49   56    62];
               else
                  selChan1={'CP2','Cz'};%{'CP2','CPz','Cz'};
                  ch_loc1=[42   65];% [42    58    65];
               end
         
               PlotAmp_M1_SIM(x1,index,Sa,startEph,endEph,St,ch_loc1,selChan1,com,compSet,stimSet);

               for st=1:St

                  WI=selected_TW(st,1);
                  figure('Renderer', 'painters', 'Position', [10 10 750 350])

                  subplot(1,2,1);
                  topoplot(squeeze(meantop_amp(st,:)),chanLoc)

                  title(['Topography map, ClustNo.', int2str(WI),', ', stimSet{st},]);
                  set(gca,'fontsize',12);
                  colorbar;
                  caxis([-2 2]);

                  subplot(1,2,2)
                  imagesc(sel_innerCorr(st).data);
                  title(['Samples Correlation']);
                  xlabel('Sample #');
                  ylabel('Sample #');
                  set(gca,'fontsize',12);
                  colorbar;
                  caxis([-1 1]);

               end

            end

         end

         if com==1
            selChan={'P2','P6','PO4'};
            ch_loc=[49   56    62];
         else
            selChan={'CP2','CPz','Cz'};
            ch_loc=[42    58    65];
         end


         selTWs(com).data=selTWs_ms;
         disp(compSet{com})
         disp(selTWs(com).data)

         compGroup_CC(count).comp(com).innSimm(k).data=inSim; %  access innSim(st,g), st =stimulus, g= group
         compGroup_CC(count).comp(com).sel_TW(k).data=selected_TW;
         compGroup_CC(count).comp(com).idx(k).data=Clu_idx;
         compGroup_CC(count).comp(com).sel_TW_ms(k).data=selTWs(com).data; % ms
         compGroup_CC(count).comp(com).meantop(k).data=meantop_amp; % access meantop_amp(st,:,g)


         [SPSS_tab_avg]=ERP_statTable2_100(selected_TW,inData,ch_loc,Subj,St,G);
         SPSStab_avg(count).comp(com).data=SPSS_tab_avg;


         [ranova_tbl]=ranova_ERP_100(SPSS_tab_avg);
         ranovatbl_all(count).comp(com).NC(k).data=ranova_tbl;


         stim_pvalue(count,com,k,:)=ranova_tbl.pValue(3);
         chan_pvalue(count,com,k,:)=ranova_tbl.pValue(5);
         intStCh_pvalue(count,com,k,:)=ranova_tbl.pValue(7);


      end

   end
end

save compGroup_CC_test.mat compGroup_CC;

all_pvalue=[stim_pvalue,chan_pvalue,intStCh_pvalue];
% disp(all_pvalue);

for count=1:Maxcount
   for com=1:Comp
      for g=1:G
         for st=1:St
            for k=a:b
               List_allCorr_CC(com,count,st,k-1,g)=compGroup_CC(count).comp(com).innSimm(k).data(st);
            end
         end
      end
   end
end

toc


%% Statistical power analysis for within factors (Stim x Chann) -----------


function [ranova_tbl]=ranova_ERP_100(SPSS_tab_avg)

Rdata=SPSS_tab_avg(:,2:7);

Group={'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1';'G1'};

varNames={'Group','St1_ch1','St1_ch2','St1_ch3','St2_ch1','St2_ch2','St2_ch3'};

tbl = table(Group,Rdata(:,1),Rdata(:,2),Rdata(:,3),Rdata(:,4),...
   Rdata(:,5),Rdata(:,6),'VariableNames',varNames);

factNames = {'Stim','Chan'};

within_R = table({'St1';'St1';'St1';'St2';'St2';'St2'},{'Ch1';'Ch2';'Ch3';'Ch1';'Ch2';'Ch3'},'VariableNames',factNames);

rm = fitrm(tbl,'St1_ch1-St2_ch3~1','WithinDesign',within_R);

[ranova_tbl] = ranova(rm, 'WithinModel','Stim*Chan');
end

