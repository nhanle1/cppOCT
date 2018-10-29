function [IMG] = SSOCT_F_PhCompV3(IMG,IMGStt,Thr,do_medianshift)
% Phase-compensation without K-clock jitter fix for SSOCT, V10012016 Shaozhen Song
%
% Input:   Bs: Raw data of a cluster of B-scan. Dimension: Z-X-Repeat
%          K: K-space
%          KES: Even-spaced K
%          IMcropRg: [top bottom] to crop Z range
%          Thr: Threshold for low signal
% Output:   IMG: complex data after phase compensation

nZ=size(IMG,1);
nX=size(IMG,2);

Zvector=(IMGStt:IMGStt+nZ-1);
DetSft=-10:2:8;

for iR=2:size(IMG,3)
    %%
    pd=angle(IMG(:,:,iR-1).*conj(IMG(:,:,iR)));
    pd(250:end,:)=NaN;
    pd(abs(IMG(:,:,iR))<Thr)=NaN;  % Threshold to remove low OCT signal
    %fix wrapping
    
    indWrap=nanmedian(pd.^2)>pi;
    pdWrap=pd(:,indWrap);
    pdWrap(pdWrap<0)=pdWrap(pdWrap<0)+2*pi;
    pd(:,indWrap)=pdWrap;
    %determine the order wrapping
    phdifDet=bsxfun(@plus,repmat(pd,[1,1,numel(DetSft)]),permute(DetSft,[3,1,2]))./repmat(Zvector',[1,nX,numel(DetSft)]);
    phdifDetM=nanmedian(phdifDet,1);
    phdifDetM(isnan(phdifDetM))=0;
    phdifDetV=nansum(abs(bsxfun(@minus,phdifDet,phdifDetM)),1);
    [~,I]=min(phdifDetV,[],3);
    phE=phdifDetM(sub2ind([1,nX,numel(DetSft)],ones(1,nX),1:nX,I));
    phComp=(phE'*Zvector)';
    DetDrift=2*DetSft(I)*pi/3;
    phComp=bsxfun(@plus,phComp,DetDrift);
    IMG(:,:,iR)=IMG(:,:,iR).*exp(1i*(phComp));
    %%
    %Remove bult motion
    if do_medianshift
        pd=angle(((IMG(:,:,iR-1)).*conj(IMG(:,:,iR))));
        pd(abs(IMG(:,:,iR))<Thr)=NaN;
        Phmed=nanmedian(pd);
        Phmed(isnan(Phmed))=0;
        IMG(:,:,iR)=IMG(:,:,iR).*repmat(exp(1i*Phmed),nZ,1);
    end
end
