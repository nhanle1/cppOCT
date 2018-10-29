function [proj_c,projInten,ColorMap] = ProjMaxDepth_v3(IMS,gwin,MaxAve,GainInt,cmap,RgDyn)
% Max+Depth projection, V10012016 Shaozhen Song, credits: Jingjiang Xu
%
% Input:    IMS: Cropped 3D volume. Projection will be performed on 1st dim
%           cmap:colormap
% Output:   proj:indexed RGB image, color-coded from depth

nX=size(IMS,2);
IMsort=sort(IMS,1,'descend');
projInten=squeeze(mean(IMsort(1:MaxAve,:,:),1));

% projInten=squeeze(max(IMS));

if gwin
    Gf = fspecial('gaussian',gwin,2);  % Gaussian filter to be used
    projInten=imfilter(projInten,Gf,'same');
end
%% Generate DepthMap
[~,I]=max(IMS);
projDepth = squeeze(I)/size(IMS,1);

projDepth=(projDepth-RgDyn(1))/(RgDyn(2)-RgDyn(1));
projDepth(projDepth<0)=0;
projDepth(projDepth>1)=1;

L = size(cmap,1);
Gs = round(interp1(linspace(0,1,L),1:L,projDepth));
ColorMap = reshape(cmap(Gs,:),[size(Gs) 3]);
%% adjust intensity map
% projInten=bsxfun(@minus,projInten,mean(projInten,2));
% projInten=projInten-imgaussfilt(projInten,[160 1]);
if GainInt~=0
    projInten=mat2gray(projInten);
    projInten=imadjust(projInten,GainInt,[0 1]);
end
%%

IntenMap = repmat(projInten, [1,1,3]);
proj_c = ColorMap.*IntenMap;
end