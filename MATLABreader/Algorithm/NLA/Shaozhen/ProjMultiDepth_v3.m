function ProjMultiDepth_v3(IMS,gwin,MaxAve,GainInt,BandLines,filenameSeg1,cut1,cut2)
% Max+Depth projection, V10012016 Shaozhen Song, credits: Jingjiang Xu
%
% Input:    IMS: Cropped 3D volume. Projection will be performed on 1st dim
%           cmap:colormap
% Output:   proj:indexed RGB image, color-coded from depth
% L_band=floor(size(IMS,1)/nBand);

nY=size(IMS,2);
nX=size(IMS,3);
nBand=numel(BandLines);
for iBand=1:nBand
    if iBand==1
        IMband=IMS(1:BandLines(iBand),:,:);
    else
        if BandLines(iBand)>size(IMS,1), BandLines(iBand)=size(IMS,1);end
        IMband=IMS(BandLines(iBand-1)+1:BandLines(iBand),:,:);
    end 
        IMsort=sort(IMband,1,'descend');
        projInten=squeeze(mean(IMsort(1:MaxAve,:,:),1));
        newmax=projInten(cut1:end-cut2,:,:);
        projInten=imresize(newmax,[size(projInten,1),size(projInten,2)]);

if gwin
    Gf = fspecial('gaussian',gwin,2);  % Gaussian filter to be used
    projInten=imfilter(projInten,Gf,'same');
end

if GainInt
    projInten=mat2gray(projInten);
    projInten=imadjust(projInten,prctile(projInten(:),GainInt),[0 1]);
end

projInten=rot90(projInten);
if iBand==1
    projpriv=projInten;
    projpriv(projpriv<0)=0;
    projIntenDiff=projInten;
else
    projIntenDiff=projInten-projpriv;
    projIntenDiff(projIntenDiff<0)=0;
    projpriv=projInten;
end

% projIntenDiff=imadjust(projIntenDiff);
% bandcolor=permute(cmap(iBand,:),[1,3,2]);
% projC=repmat(projInten,[1,1,3]).*repmat(bandcolor,[nX,nY,1]);

% imwrite(projC,['.\color\','proj_C_',filenameSeg1(1:end-4),'_Layer',num2str(iBand),'.tiff']);

imwrite(projInten,['.\grayscale\','proj_bw_',filenameSeg1(1:end-4),'_Layer',num2str(iBand),'.tiff']);
imwrite(projIntenDiff,['.\diff\','proj_d_',filenameSeg1(1:end-4),'_Layer',num2str(iBand),'.tiff']);
end
end