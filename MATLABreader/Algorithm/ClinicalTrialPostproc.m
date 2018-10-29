matlist     =dir('*.bind');
crop=81:680;
nZ=size(crop);
datasize=[1024 500 500];
for i=1:1%numel(matlist)
    [~,matfile,~]     =fileparts(matlist(i).name);
    %[~,tfold,~]   =fileparts(matlist(i).name);
    %[~,tfold,~]   =fileparts(tfold);

    [BL,STR]=CppReadBin(matfile,crop, datasize);
    tfold=['\',matfile,'\'];
    NLproject;

%     [BL,STR]=CppReadBin(matfile,crop);
%     folder=[tfold,'BL\'];
%     ProjectZeach;
%     
%     [BL,STR]=CppReadBin(matfile,crop);
%     folder=[tfold,'STR\'];
%     ProjectZeachSTR;
end

