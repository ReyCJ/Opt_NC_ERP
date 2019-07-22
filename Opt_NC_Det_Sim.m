
% By: Reza Mahini Reza Mahini 2018 updated Apr 2019



function [Opt_NC,Opt_TW,corr_opt_NC,std_Opt_NC,corr_all]=Opt_NC_Det_Sim(compInfo_CC,a,b,St,g,Comp,Thr)


stimSet={'Cond1','Cond2'};
% % groupSet={'RS','HC'};
comp={'N2','P3'};


% for g=1:G % number of groups


count=1;


figure(g+2*Comp-2)


for st=1:St % number of stimuli
    for k=a:b % test range of clustering
        %       std_stim(k-(a-1))=CompMethod_all.std(m).stGr(k).data(st,g); % getting std_value for each stimulus wwithin range of clustering
        std_stim(2,k-(a-1))=compInfo_CC.std(k).data(st); % getting std_value for each stimulus wwithin range of clustering
        corr_avg(2,k-(a-1))=mean(mean(compInfo_CC.Corr(k).data.gr(st).data,2));
        std_stim(1,k-(a-1))=k;
        corr_avg(1,k-(a-1))=k;
    end
    corr_all(st,:)=corr_avg(2,:);
 
    clear min;
    clear max;
    [~,id]=min(std_stim(2,:)); % selecting minimum std or high correlation
    [~,M_id]=max(corr_avg(2,:));
    [Mval,~]=max(std_stim(2,:));
    [mval,~]=min(std_stim(2,:));
    [Maval,NC]=max(corr_avg(2,:));
    
 
    
    %%  Selecting optimal number of clusters --------------------------
    % ************ PLEASE SET THESE PARAMETERS CAREFULLY *********************

corThr=0.02; % Correlation stability parameter
maxThr=Thr; % Optimal selecting parameter
      
    fi=0; % not detected
    for i=2:size(corr_avg,2)
        if (corr_avg(2,i)>=Thr) && (abs(corr_avg(2,i-1)-corr_avg(2,i))<corThr)
            OPNC(g).stim(st).data=[corr_avg(1,i),corr_avg(2,i)];
            fi=1;
            break
        end
    end
    
    corThr1=corThr;  % We need to be flaxable for finding second choice
    while fi==0 && corThr1<0.1 % maximum tolorence
        
        for i=2:size(corr_avg,2)
            if (corr_avg(2,i)>=maxThr) && (abs(corr_avg(2,i-1)-corr_avg(2,i))<corThr1)
                OPNC(g).stim(st).data=[corr_avg(1,i),corr_avg(2,i)];
                fi=1;
                break
            end
        end
        maxThr=maxThr-0.0005; % decreasing the maximum threshold
        corThr1=corThr1+0.00001;
        
    end
    
    %%%
    [mival,~]=min(corr_avg(2,:));
    
    id=OPNC(g).stim(st).data(1); % OPNC ID
    
    %% Extracting information about selected optimal number of cluster
        
    if id~=0
        mval= std_stim(2,id-(a-1));
        Opt_NC(st,g)=OPNC(g).stim(st).data(1);           %   id+(a-1);  % method x stimulus
        corr_opt_NC(st,g)=OPNC(g).stim(st).data(2);      % M_id+(a-1);
        Opt_TW(st,:,g)=compInfo_CC.sel_TW(id).data(st,:,g); % optimal TW basen on optimum NC
        %       Opt_TW_Corr.stim(st).data=compInfo_CC.sel_TW(M_id+(a-1)).data(st,:,g); % optimal TW basen on optimum NC
        std_Opt_NC(st,g)=mval;
    else
        % there is no optimal
        msg = 'There is no optimal NC, please check the code!';
        error(msg)
    end
    
  %% Plot the inner similarity and OpNC
  
 
    subplot(1,St,st,'align');
    plot(corr_avg(2,:),'-o','LineWidth',1,'MarkerSize',5);
    xticks(1:(b-a+1));
    xticklabels(a:b);
    hold on
    %         plot([5,5],[mival,Maval],'--k');
    plot([(OPNC(g).stim(st).data(1)-a)+1,(OPNC(g).stim(st).data(1)-a)+1],[0 ,1.05],'--k');
    plot([0,(b-a)+1],[Thr,Thr],'--r');
    % %         gtext('Optimal number of clusters')
    title(['Optimal NC for ', comp{Comp}, ' Stim= ',stimSet{st}]);
    xlabel('Number of cluster #');
    ylabel('Correlation');

    count=count+1;
    set(gca, 'fontsize',12);
%     set(gcf,'outerposition',get(0,'screensize'));
    ylim([0 1.05]);

    
end

end
