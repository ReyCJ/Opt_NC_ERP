%--------------------- Using FCM for ERP Data ---------------------------
% The ERP data has well structure for clustering time sammples

function [F_ERP]=FCM_GAERP(x,k)


% ------------------------------- clustering ------------------------------
options=[1.1 120 NaN 0];
[center, U, obj_fcn] = fcm(x, k,options); % FCM with Matlab function

for i=1:k
    maxU = max(U);
    membership(i).idx = find(U(i, :) == maxU);
    
end

% ------------------------------ Membership -----------------------

[p,q]=size(membership);
d=1;
for i=1:q % number of clusters
    a=membership(1,i).idx;
    [a1,a2]=size(a);
    for j=1:a2
        FCM_m(a(j))=i;
    end
end

F_ERP=FCM_m';
c_idx=FCM_m';

end



