function [Opt_NC,Opt_TW,InnSim_opt_NC,InnSim_all]=Opt_NC_Det_SIM(compInfo_CC,a,b,St,g,Comp,Thr)


stimSet={'Cond1','Cond2'};
comp={'N2','P3'};


% for g=1:G % number of groups

count=1;
P=figure(g+2*Comp-2);
set(P,'Renderer', 'painters', 'Position', [50 100 850 450]);


for st=1:St % number of stimuli
   for k=a:b % test range of clustering
      InnSim(2,k-(a-1))=compInfo_CC.innSimm(k).data(st);
      InnSim(1,k-(a-1))=k;
   end
   InnSim_all(st,:)=InnSim(2,:);

   clear min;
   clear max;
   [~,M_id]=max(InnSim(2,:));
   [Maval,NC]=max(InnSim(2,:));


   %%  Selecting optimal number of clusters --------------------------
   % ************ PLEASE SET THESE PARAMETERS CAREFULLY *********************

   corThr=0.02; % Correlation stability parameter
   maxThr=Thr; % Optimal selecting parameter

   fi=0; % not detected
   for i=2:size(InnSim,2)
      if (InnSim(2,i)>=Thr) && (abs(InnSim(2,i-1)-InnSim(2,i))<corThr)
         OPNC(g).stim(st).data=[InnSim(1,i),InnSim(2,i)];
         fi=1;
         break
      end
   end

   corThr1=corThr;  % We need to be flaxable for finding second choice
   while fi==0 && corThr1<0.1 % maximum tolorence

      for i=2:size(InnSim,2)
         if (InnSim(2,i)>=maxThr) && (abs(InnSim(2,i-1)-InnSim(2,i))<corThr1)
            OPNC(g).stim(st).data=[InnSim(1,i),InnSim(2,i)];
            fi=1;
            break
         end
      end
      maxThr=maxThr-0.0005; % decreasing the maximum threshold
      corThr1=corThr1+0.00001;
   end

   %%%
   [mival,~]=min(InnSim(2,:));

   id=OPNC(g).stim(st).data(1); % OPNC ID

   %% Extracting information about selected optimal number of cluster

   %  ****** important information about opNC and TW

   if id~=0
      Opt_NC(st,g)=OPNC(g).stim(st).data(1);           %   id+(a-1);  % method x stimulus
      InnSim_opt_NC(st,g)=OPNC(g).stim(st).data(2);      % M_id+(a-1);
      Opt_TW(st,:,g)=compInfo_CC.sel_TW(id).data(st,:,g); % optimal TW basen on optimum NC

   else
      msg = 'There is no optimal NC, please check the code!';
      error(msg)
   end

   %% Plot the inner similarity and OpNC


   subplot(1,St,st,'align');
   plot(InnSim(2,:),'-o','LineWidth',1,'MarkerSize',5);
   xticks(1:(b-a+1));
   xticklabels(a:b);
   hold on
   plot([(OPNC(g).stim(st).data(1)-a)+1,(OPNC(g).stim(st).data(1)-a)+1],[0 ,1.05],'--k');
   plot([0,(b-a)+1],[Thr,Thr],'--r');
   title(['Optimal NC for ', comp{Comp}, ' Stim= ',stimSet{st}]);
   xlabel('Number of cluster #');
   ylabel('Inner-similarity');

   count=count+1;
   set(gca, 'fontsize',12);
   ylim([0 1.05]);


end

end
