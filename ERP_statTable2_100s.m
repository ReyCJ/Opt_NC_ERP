% Prepare Data for SPSS analysis ------------------------------------------
% By Reza Mahini may 2017

function [SPSS_tab_avg]=ERP_statTable2_100s(selected_TW,inData,anzChan,Subj,Stim,G)


SPSS_tab_avg=[];
nCh=length(anzChan);

for g=1:G % groups
    
    for s=1:Subj % subjects
        SPSS_tab_avg(Subj*(g-1)+s,1)=g;
        for st=1:Stim % stimuli
            time_window_tab(st,:)=selected_TW(st,2:3,g);
            mean_pow_sel_chan(:,st,s,g)=squeeze(mean(inData(anzChan,time_window_tab(st,1):time_window_tab(st,2),st,s,g),2));
            SPSS_tab_avg(Subj*(g-1)+s,nCh*(st-1)+2:st*nCh+1)=[squeeze(mean_pow_sel_chan(:,st,s,g))'];
            % SPSS_tab_avg(Subj*(g-1)+s,st*nCh-(nCh-1):st*nCh)
        end
    end
end


% xlswrite([filepath 'tab4statAnals.xls'],SPSS_tab_avg(:,:)); % ****

end
 