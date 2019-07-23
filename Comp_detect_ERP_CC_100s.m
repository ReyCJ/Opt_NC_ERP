% March 2017 by Reza Mahini
% Updated in march 2019, 

function [sel_info,comp_pow,SelStd_val,innerCorr,winnID]=Comp_detect_ERP_CC_100s(c_idx1,x,Chanloc,K,Sa,St,v,w,g,com,stimSet,compSet)

mclThr=2; % sample threshoold for merging clusters

% % stimSet={'St1','St2'};
% % % % groupSet={'RS','HC'};
% % compSet={'N2','P3'};

%chan=size(Chanloc,2);
chan=size(x,2); % check it pls


% for g=1:G
    
    for st=1:St
        
        analysis_area=c_idx1(st*Sa-(Sa-1):Sa*st,1); % Select area for detecting for single number of clusters--> c_idx1(st*Sa-(Sa-1):Sa*st,1);
        x1=x(st*Sa-(Sa-1):Sa*st,:);       
        cnt=0; %
        
        std_value=[];
        save_detec=[];
        std_v=[];
        s_detec=[];
        stdVal=[];
        Rank=[];
        
        
        % x11=squeeze(x1(:,:,g)); % samples x1 channels x1 group
        
        for f=1:K % Checking for each cluster in the area
            
            s_detec=[];
            innerSim=[];
            innerSim1=[];

            
            d1=[]; A=[]; b1=[]; % parameters adjusment
            
            b1=find(analysis_area()==f);
            
            c_info(f).data=b1;
            si(f)=size((b1),1);
            
            % Saving the areas information for each cluster ---------------
            j=1;
            s_detec(j).data=[]; % Saving the duration of each cluster.
            sc_mem=c_info(f).data;
            for i=2:length(c_info(f).data)
                if sc_mem(i)-sc_mem(i-1)<mclThr % 3 is threshold for merging clusters
                    s_detec(j).data(size(s_detec(j).data,2)+1)=sc_mem(i-1);
                else
                    s_detec(j).data(size(s_detec(j).data,2)+1)=sc_mem(i-1);
                    j=j+1;
                    s_detec(j).data=[];
                end
            end
            
            
            % Here we need to know which clusters are inthe range ---------
            
            % Test for cluster area is suitable or in range or not , to consider
            % which clusters are suitable for considering
            for t=1:size(s_detec,2)
                if isempty(s_detec(t).data)||(s_detec(t).data(1)-v<0 && s_detec(t).data(size(s_detec(t).data,2))-v<0)...
                        ||(s_detec(t).data(1)-w>0 && s_detec(t).data(size(s_detec(t).data,2))-w>0)
                    s_detec(t).data=[];
                end
            end
            
            % Acceptable areas for each cluster selection -----------------
            cnt=0; %
            thr=10; % min number of samples for selecting comp 1samp=4.5ms this ex1ample 

            for ii=1:length(s_detec) % for detected areas which one is >25
                if length(s_detec(ii).data)>=thr % size>25 (thereshold for acceptable area)
                    cnt=cnt+1;
                    save_detec(f).sel(cnt).data=s_detec(ii).data;
                end
                
            end
            if cnt==0
                save_detec(f).sel=[]; %
                
            end
            %%test to know if no comp is found
            flag=0;
            for r=1:size(save_detec,2)
                fl=any( structfun(@isempty, save_detec(r)));      % isempty(save_detec); % _______________Change____________________
                if fl~=1
                    flag=1;
                end
            end
            
            if flag==0 % if no first component , we have to decrease threshold
                while cnt==0 && thr>=4 % min threshold
                    for ii=1:length(s_detec) % for detected areas which one is >25
                        if length(s_detec(ii).data)>=thr% size>35 (thereshold for acceptable area)
                            cnt=cnt+1;
                            save_detec(f).sel(cnt).data=s_detec(ii).data;
                        end
                        
                    end
                    thr=thr-1;
                end
            end
            
            % Ranking the selected clusters ------------------------------------------
            %for i1=1:size(sve_detec,2) % 6 clusters maybe no need this ?
            emp=isempty(save_detec); % for the clusters are not satissfied with condition
            if emp~=1 % if there is any >25 memeber
                te=isempty(save_detec(f).sel);
                if te~=1
                    for i2=1:size(save_detec(f).sel,2) % area in each cluster maybe more than two one area has that condition
                        % calculating STD for the selected clusters -------------------------------
                        f1=save_detec(f).sel(i2).data;
                        innerSim = corr(x1(f1,:)'); % correlation in specific area
                        
                        
                        % innerSimilarity calculation
                        for t=1:size(innerSim,2)
                            A=innerSim(t,:);
                            item=A(t);
                            A(t)=[];
                            for jj=1:length(A)
                                d1(t,jj)=abs(item-A(jj)); 
                            end
                        end
                        d(f).dist(i2).data=mean(d1,1);
                        std_v(f).sel(i2)=std(d(f).dist(i2).data);
                        
                    end % the end of clusters with several areas ------------------
                    
                    clear min
                    std_value(f)=min(std_v(f).sel); % finding minimum if morethan one area >25 for one cluster
                else if te==1
                        std_value(f)=10 ; % Avoiding from select empty clusters (rejecting sparse clusters)
                    end
                end
            else if emp==1
                    std_value(f)=10; % Avoiding from select empty clusters (rejecting sparse clusters)
                end
            end
            
        end
        
        [stdVal,Rank]=sort(std_value);   % ,'descend'); 3 2 1
        
        
% *** Detecting strongest area -------------------------------------------
        
        clear min; % first high corralated cluster selection
        [w_max1_memb1,winner_id1]=min(std_value); % Finding which cluster is strong in this area
        if length(std_v(winner_id1).sel)>1 % if there is more than 1 candidate for any cluster
            [val sel_ID]=min(std_v(winner_id1).sel); % find again min one
            sc_memb1=save_detec(winner_id1).sel(sel_ID).data; % ex1tract time samples
        else
            sc_memb1=save_detec(winner_id1).sel.data; % Finding real time points position
        end
        
        temp=std_value(winner_id1);
        std_value(winner_id1)=10; % avoiding to select same one again
        

        % *** Saveing needed information ----------------------------------------
        
        % saving first component info
        
        SelStd_val(1,st,g)=w_max1_memb1;
        duration1=length(sc_memb1);
        start_point1=sc_memb1(1);
        end_point1=sc_memb1(length(sc_memb1));
        memb_comp1=sc_memb1;
        power_amp1=mean(x1(sc_memb1,:),1); % power of amplitude
        innerSim1 = corr(x1(sc_memb1,:)'); % correlation in specific area
        innerCorr(1).sti(st).data=innerSim1;
        
 
% %         % Plot the component1 --------------------------------------------
% %         
% %         figure('Renderer', 'painters', 'Position', [10 10 900 400])
% % 
% % %         set(gcf,'outerposition',get(0,'screensize'));
% %         subplot(1,2,1);
% %         topoplot(power_amp1,Chanloc)
% %         title(['Topography for ', compSet{com}, ', ', stimSet{st},', map ', int2str(winner_id1)]);
% %         set(gca,'fontsize',14);
% %         colorbar;
% %         cax1is([-1 1]);
% %         
% %         subplot(1,2,2)
% %         imagesc(innerSim1);
% %         title(['Correlation for ', compSet{com}, ', ', stimSet{st},', map ', int2str(winner_id1)]);
% %         x1label('Sample #');
% %         ylabel('Sample #');
% %         set(gca,'fontsize',14);
% %         colorbar;
% %         cax1is([-1 1]);
% % % %         
        
        
        
        % saving second component info
        l=0;
        for i=1:K
            emp_detect=isempty(save_detec(i).sel);
            if emp_detect~=1
                l=l+1;
            end
        end
        
        if l>1
            
            clear min; % second high corralated cluster selection
            [w_max1_memb2,winner_id2]=min(std_value); % Finding which cluster is strong in this area
            if length(std_v(winner_id2).sel)>1 % if there is more than 1 candidate
                [val sel_ID]=min(std_v(winner_id2).sel); % find again min one
                sc_memb2=save_detec(winner_id2).sel(sel_ID).data; % ex1tract time samples
            else
                sc_memb2=save_detec(winner_id2).sel.data; % Finding real time points position
            end
            
            
            SelStd_val(2,st,g)=w_max1_memb2;
            duration2=length(sc_memb2);
            start_point2=sc_memb2(1);
            end_point2=sc_memb2(length(sc_memb2));
            memb_comp2=sc_memb2;
            power_amp2=mean(x1(sc_memb2,:),1); % power of amplitude
            innerSim2 = corr(x1(sc_memb2,:)'); % correlation in specific area
            innerCorr(2).sti(st).data=innerSim2;
            comp_pow(st,2,:,g)=power_amp2;
             winnID(st,2,g)=winner_id2;

        else
           
            duration2=0;
            start_point2=0;
            end_point2=0;
            memb_comp2=0;
            winner_id2=0;
            comp_pow(st,2,:,g)=zeros(1,chan);
            winnID(st,2,g)=winner_id2;

            
        end
        
    % plot -------------------------------------------------------------

% % 
% %         figure('Renderer', 'painters', 'Position', [10 10 900 400])
% %         
% % %             set(gcf,'outerposition',get(0,'screensize'));
% %             subplot(1,2,1);
% %             %figure
% %             topoplot(power_amp2,Chanloc) % it is true
% %         title(['Topography2 for ', compSet{com}, ', ', stimSet{st},', map ', int2str(winner_id2)]);
% %             set(gca,'fontsize',14);
% %             
% %             colorbar;
% %             cax1is([-1 1]);
% %             
% %             subplot(1,2,2)
% %             imagesc(innerSim2);
% %         title(['Correlation2 for ', compSet{com}, ', ', stimSet{st},', map ', int2str(winner_id2)]);
% %             x1label('Sample #');
% %             ylabel('Sample #');
% %             set(gca,'fontsize',14);
% %             
% %             colorbar;
% %             cax1is([-1 1]);
% %             
            
            
            % ----------------------------------------------------------
    
% saving cluster1-cluster2-start1-end1-duration1-start2-end2-duration2 ----

        sel_info(st,1:8,g)=[winner_id1, start_point1, end_point1, duration1,...
            winner_id2,start_point2, end_point2, duration2];  % save needed info.
        
        comp_pow(st,1,:,g)=power_amp1;
         winnID(st,1,g)=winner_id1;

        
    end
% end