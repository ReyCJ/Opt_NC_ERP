% The Function will prepare primary model of data by maeging all stimuli
% together
% Code by; Reza Mahini Oct 2017

function [ERP_Subj,inDaGA]=Data_Preparing(Data,Sabj,St,Sa,G) % N=number of subjects
% S= Number of stimuli, SA= number of samples

for g=1:G
    for subj=1:Sabj
        for st=1:St
            ERP_Subj(st*Sa-(Sa-1):st*Sa,:,subj,g)=squeeze(Data(:,:,st,subj,g))';
            % samples x channels x subjects x groups
        end
    end
    inGA(:,:,g)=squeeze(mean(ERP_Subj(:,:,:,g),3));
    
    
    % samples x channels x groups
end
for g=1:G
    inDaGA(g*size(inGA,1)-size(inGA,1)+1:g*size(inGA,1),:)=squeeze(inGA(:,:,g));
end
end




