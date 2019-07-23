%  using of SOM (self-orgnizating map)neural network aim clustering
%  one stimuli as demo. the motivation of this demo we have used matlab
%  SOM toolbox for clustering ERP data
function [S_ERP]=SOM_GAERP(x,k)

% x is grand average data e.g. 4800 (observation )x30 (feature)

x=x';
m=1; n=k;

%----------------------Clustering with SOM ----------------------


net = selforgmap([m n]); % 2-dimension layer of 6 neurons arranged in an 2x3
% hexagonal grid for this example.

net.trainParam.showWindow=0;

[net,tr] = train(net,x); % Now the network is ready to be optimized with *train*.
% % % nntraintool

y = net(x);
cl_idx = (vec2ind(y)); % cluster indexes by SOM
S_ERP=cl_idx';

% --------------------------- membership plot --------------------------
% % % figure;
% % % imagesc(cl_idx');
% % % title(['The Membership for all Stimulus Grand Averaged Data']);
% % % set(gca,'fontsize',13);

% % % close all;

% % % nntraintool('close');
end



