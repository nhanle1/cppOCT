
% folder=['\BLslice\'];
%folder=['\ED\',num2str(i),'\'];
if ~exist(matfile,'dir'),mkdir(matfile),end

outputdir=[cd,'\',folder];
mkdir(outputdir);

do_NM=1;
do_proj=1;

do_showSeg=0;
[nZ,nX,nY] = size(STR);


%% Save Dicom
BLn=double(STR);
RgFlowNM_ed=prctile(double(BLn(1:47:end)),[1 100]);
BLn=uint16(65536*(BLn-RgFlowNM_ed(1))/(RgFlowNM_ed(2)-RgFlowNM_ed(1)));




%% Projection
if do_proj
    Info.nZshift=-5;
    Info.nZseg=195; %total
    
    Info.WinSize=70;
    %Info.WinSize=200;
    thrsh_ratio=0.90;
    medfilt=[5 5];
    gaussfilt=[2 10];
    StrRg=[0 65535];
    seg_top=zeros(nX,nY);
    pb = ProgressBar(nY);
    parfor iY=1:nY
        cimgc=imgaussfilt(medfilt2(STR(:,:,iY),...
            medfilt),gaussfilt);
        %cimgc=imgaussfilt(double(edge(STR(:,:,iY),'Roberts')),gaussfilt);
        cimgMax=max(cimgc);
        for iX=1:nX
            vol=cimgc(:,iX);
            loc = find(vol>cimgMax(iX)*thrsh_ratio,1,'first');
            while isempty(loc)
                cimgMax(iX)=cimgMax(iX)*thrsh_ratio;
                loc = find(vol>cimgMax(iX)*thrsh_ratio,1,'first');
            end
            seg_top(iX,iY)=loc(1);
        end
        %             seg_top(:,iY)=smooth(seg_top(:,iY),0.7,'rloess');
        pb.progress;
    end
    pb.stop;
    %     seg_top     =smoothdata(seg_top,1,'rlowess',nX*0.3);
    %     seg_top     =smoothdata(seg_top,2,'rlowess',nY*0.3);
    seg_f=medfilt2(seg_top,medfilt,'symmetric');
    seg_fs=seg_f+Info.nZshift;
    %%
    if do_showSeg
        seg_f2=seg_fs+Info.nZseg;
        for iY=1:5:nY
            cla
            cimg=BL(:,:,iY);
            %             imagesc(cimg,StrRg);
            imagesc(cimg);
            hold on,plot(seg_fs(:,iY),'--r'),
            plot(seg_f2(:,iY),'--r'),
            pause(0.05)
        end
    end
    %% Projection
    [O,P,Q]=ndgrid(1:Info.nZseg,1:nX,1:nY);
    O=bsxfun(@plus,O,permute(round(seg_fs),[3,1,2]));
    O(O>size(STR,1))=size(STR,1);
    O(O<1)=1;
    ind=sub2ind(size(STR),O,P,Q);
    clear O P Q
    % Projection
    MaxAve=3;
    gwin=0;
    GainInt=0;
    IMF=squeeze(BLn);
    IMS=double(IMF(ind))/65535;
    %     Rg=prctile(IMS(1:20:end),[80 99.5]);
    %     IMS(IMS<Rg(1))=Rg(1);IMS(IMS>Rg(2))=Rg(2);
    %     IMS=(IMS-Rg(1))/(Rg(2)-Rg(1));
    %         cmap=flip(isolum(256),1);
    %         [Pj_c,Pj_g,CM]=ProjMaxDepth_v3(IMS,gwin,MaxAve,[0.05 0.95],cmap,[0.02 0.75]);
    %         %Pj_g=imadjust(Pj_g,[0.3 1],[0 1]);
    %         Pj_ci= my_interp(Pj_c, [1 -5 20 20 -5 1]/32);
    %         Pj_gi= my_interp(Pj_g, [1 -5 20 20 -5 1]/32);
    %         figure;imshow(Pj_ci);
    %         figure;imshow(Pj_gi);
    
    for i=1:1:Info.nZseg-5
        temp    =(mat2gray(squeeze(mean(IMS(i:i+2,:,:),1))));
        Rg=prctile(temp(1:27:end),[7 100]);
        temp(temp<Rg(1))=Rg(1);temp(temp>Rg(2))=Rg(2);
        temp=(temp-Rg(1))/(Rg(2)-Rg(1));
        imwrite(imrotate(squeeze(temp),180),...
            [outputdir,'c',num2str((i)),'.png'],...
            'png');
    end
end

clearvars thres temp maxAve gwin GainInt