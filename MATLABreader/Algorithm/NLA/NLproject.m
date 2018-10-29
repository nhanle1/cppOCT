% folder='\burnpig1\';
%folder=['\ED\',num2str(i),'\'];
if ~exist(matfile,'dir'),mkdir(matfile),end
outputdir=[cd,'\',matfile,'\'];
do_NM=1;
do_proj=1;
saveDicom=1;
do_showSeg=1;


%% Save Dicom
% STR=uint16(65536*(STR-min(STR(:)))/(max(STR(:))-min(STR(:))));
% BL=uint16(65536*(BL-min(BL(:)))/(max(BL(:))-min(BL(:))));


[nZ,nX,nY] = size(STR);

for i=3:5
    BLn     =zeros(nZ,nX,nY);
    filename=num2str(i);
    if do_NM
        BLn=(imgaussfilt(double(BL),1.2)-imgaussfilt(double(STR),1.5)*i/10);
        %% modify norm
%         for iY=1:nY
%             BLn(:,:,iY)     =NLnorm(BL(:,:,iY),STR(:,:,iY),4,4,i/10);
%         end
        %% end modify norm
        
        RgFlowNM_ed=prctile(double(BLn(1:37:end)),[10 100]);
        BLn=uint16(65536*(BLn-RgFlowNM_ed(1))/(RgFlowNM_ed(2)-RgFlowNM_ed(1)));
    end
    if i==5||i==10
        fprintf('saving DICOM....');
        dicomwrite(reshape(BLn, [nZ,nX,1,nY]),[outputdir,filename,'_NM_Flow_ed','.dcm']);
        if i==5
            dicomwrite(reshape(STR, [nZ,nX,1,nY]),[outputdir,filename,'-StruC.dcm']);
            dicomwrite(reshape(BL, [nZ,nX,1,nY]),[outputdir,filename,'-Flow_edC.dcm']);
        end
        %     else
        %         fprintf('no DICOM saved');
        %         dicomwrite(reshape(BLDNMimg_ed, [a,b,1,c]),[outputdir,filename,'_NM_Flow_ed','.dcm']);
    end
    %% Projection
    if do_proj
        nZshift=16; %start from 14 for retina and 25 for gum
        nZseg=115; %total

        WinSize=70;
        %WinSize=200;
        thrsh_ratio=0.95;
        medfilt=[5 5];
        gaussfilt=[5 5];
        StrRg=[0 65535];
        seg_top=zeros(nX,nY);
        pb = ProgressBar(nY);
        parfor iY=1:nY
            cimgc=imgradient(imgaussfilt(STR(:,:,iY),...
                3,'FilterSize' ,[13 13]));
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
        %         seg_top     =smoothdata(seg_top,1,'rlowess',nX*0.3);
        %         seg_top     =smoothdata(seg_top,2,'rlowess',nY*0.3);
        %         for iX=1:nX
        %             seg_top(iX,:)=smooth(seg_top(iX,:),0.5,'rloess');
        %         end
        pb.stop;
        fprintf('\n')
        seg_f=medfilt2(seg_top,medfilt,'symmetric');
        seg_fs=seg_f+nZshift;
        %%
        if do_showSeg
            seg_f2=seg_fs+nZseg;
            for iY=1:10:nY
                cla
                cimg=BLn(:,:,iY);
                %             imagesc(cimg,StrRg);
                imagesc(cimg);
                hold on,plot(seg_fs(:,iY),'--r'),
                plot(seg_f2(:,iY),'--r'),
                pause(0.05)
            end
        end
        %% Projection
        [O,P,Q]=ndgrid(1:nZseg,1:nX,1:nY);
        O=bsxfun(@plus,O,permute(round(seg_fs),[3,1,2]));
        O(O>size(STR,1))=size(STR,1);
        O(O<1)=1;
        ind=sub2ind(size(STR),O,P,Q);
        
        % Projection
        MaxAve=3;
        gwin=0;
        GainInt=0;
        IMS=mat2gray(BLn(ind));
        parfor iZ=1:nZseg
            IMS(iZ,:,:)=(imadjust(squeeze(IMS(iZ,:,:)),[0.5-i*0.05 ...
                1-i*0.05],[]));
        end
%         IMS=double(IMS(ind))/65535;
%         Rg=prctile(IMS(1:20:end),[70 99.98]);
%         IMS(IMS<Rg(1))=Rg(1);IMS(IMS>Rg(2))=Rg(2);
%         IMS=(IMS-Rg(1))/(Rg(2)-Rg(1));
        cmap=flip(isolum(256),1);
        %         [Pj_c,Pj_g,CM]=ProjMaxDepth_v3(IMS,gwin,MaxAve,[0 1],cmap,[0.02 0.75]);
        [Pj_c,Pj_g,CM]=ProjMaxDepth_v3(IMS,gwin,MaxAve,[0 1],cmap,[0.01 0.90]);
        %Pj_g=imadjust(Pj_g,[0.3 1],[0 1]);
        Pj_ci= my_interp(Pj_c, [1 -5 20 20 -5 1]/32);
        Pj_gi= my_interp(Pj_g, [1 -5 20 20 -5 1]/32);
        %         figure;imshow(Pj_ci);
        %         figure;imshow(Pj_gi);
        imwrite(imrotate(Pj_ci,180),[outputdir,'t',filename,'.tiff']);
%         imwrite(imrotate(Pj_gi,180),[outputdir,filename,'-Proj_g.tiff']);
    end
end