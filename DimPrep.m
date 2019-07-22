% Dimension adjustment for ERP data ---------------------------------
% By: Reza Mahini Apr 2019

% The Dimention of data must change to 
% Channel x Sample x Stim x Subject x Group 
 
function inData=DimPrep(Data,Chan,Sa,St,Subj,G)

O1= find(size(Data)==Chan);% finding position in Data for each dimention
O2= find(size(Data)==Sa);
O3= find(size(Data)==St);
O4= find(size(Data)==Subj);

if length(size(Data))==5
    O5= find(size(Data)==G);
    inData= permute(Data,[O1 O2 O3 O4 O5]);
end

inData= permute(Data,[O1 O2 O3 O4]); % Chan x Sa x St x Subj

% Test pemuted data ------------------------

disp('Original data size')
size(Data)

disp('Converted data size');
size(inData)

end