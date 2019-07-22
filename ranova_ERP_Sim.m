%% Statistical power analysis for within factors (Stim x Chann) -----------

% By: Reza Mahini Apr 2019


function [ranova_tbl]=ranova_ERP_Sim(SPSS_tab_avg)

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