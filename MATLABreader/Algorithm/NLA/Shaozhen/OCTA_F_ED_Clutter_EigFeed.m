%%% Siavash Yousefi

%% p_tis is the power of the input signal r_c
%% p_bld is the power of the filtered signal (clutter rejected signal) which
%% is showing the blood flow
%% Num_Ensembles is the number of A-scans per depth location (typically
%% 4) ensemble size
%% Lines_per_Frame is the number of A-scans in each B-image
%% Points_per_Aline is number of points in each A-scan, which is showing the
%% depth 
%% NumEV is the number of eigenvectors to be removed, typically by default
%% should be 2

function [p_tis,p_bld] = OCTA_F_ED_Clutter_EigFeed(r_c, NumEV)

Num_Ensembles = size(r_c,3);
Lines_per_Frame = size(r_c,2);
Points_per_Aline = size(r_c,1); 
p_bld=zeros(Points_per_Aline,Lines_per_Frame); 
p_tis=p_bld;
for i=1:Lines_per_Frame

    data= double(reshape(r_c(:,i,1:Num_Ensembles), Points_per_Aline, Num_Ensembles));
    
    Rcx=data'*data;
    Rcx = Rcx/Points_per_Aline;
    
    [VRcx, D] = eig(Rcx);
    
    Pkx = zeros(Num_Ensembles, Num_Ensembles);
    
    for ev=1:1:NumEV
        Pkx = Pkx +  VRcx(:,(Num_Ensembles-ev+1))*(VRcx(:, (Num_Ensembles-ev+1))') ;
    end
    
    Pk = eye(Num_Ensembles, Num_Ensembles) - Pkx;

    
    filtered_data = data*Pk'; 

    p_bld(:,i) = sum(abs(filtered_data),2)./Num_Ensembles;
    p_tis(:,i) = sum(abs(data*Pkx'),2)./Num_Ensembles;

end




