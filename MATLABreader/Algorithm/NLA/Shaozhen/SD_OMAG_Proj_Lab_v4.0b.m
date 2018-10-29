clear;clc;close all;
if ~exist('.\output\','dir'),mkdir('.\output\'),end
outputdir=[cd,'\output\'];

do_reg=0; % toggle Inter-Frame registration (for motion removal)
useref=0; % set 1 to use "refdata" as reference signal.

saveDicom=1;
show_img=1;
do_NM=1;
do_proj=1;
do_showSeg=1;

IMcropRg=100:450;
nZcrop=numel(IMcropRg);

imshowrgZ=1:nZcrop;
distcomp.feature( 'LocalUseMpiexec', false ); % parallel processing

tmp=dir('.\*OMAG.oct');
for iFile=1:size(tmp,1), fprintf(['[%d]\t',tmp(iFile).name(1:end-4),'\n'],iFile);end

% SPL=1024;
SPL=2048;
% nX=730;
% nY=730;
% nR=4;
% coefs=[-1.8429E+0	1.0561E+0	2.9183E-5	-2.7669E-8];
K=1:SPL;
KES=textread('KES.txt');
% coefs=fliplr(coefs);KES=polyval(coefs,1:SPL);
refname='ref_data';
use_autoRg=1;
RgStr=[90 130];
RgFlow=[38 65];
%%
tic

for iFile=1:size(tmp,1)%Specify data # to process
    %%
    filename=tmp(iFile).name;
    
    fid=fopen(filename);
    fseek(fid,-8,'eof');  
    Protof = fread(fid,1,'uint64');

    fseek(fid,Protof,'bof');  
    XscanRange= fread(fid,1,'double');
    YscanRange=fread(fid,1,'double');
    Xoffset=fread(fid,1,'double');
    Yoffset=fread(fid,1,'double');
    nX=fread(fid,1,'uint32');
    ImageNum=fread(fid,1,'uint32');
    nR=fread(fid,1,'double');
    FrameRate=fread(fid,1,'double');
    AlineRate=fread(fid,1,'double');
    Xcoef=fread(fid,1,'double');
    Ycoef=fread(fid,1,'double');
    camPhase=fread(fid,1,'double');
    cardPhase=fread(fid,1,'double');
%     Coef=fread(fid,4,'double');
%     KES=fread(fid,SPL,'double');
    
    nY=floor(ImageNum/nR);
    
    
    %%

%     nX=char(regexp(filename,'[_]\w[^_]*[A][_]','match'));nX=str2num(nX(2:end-2));
%     nR=char(regexp(filename,'[_]\w[^_]*[R][_]','match'));nR=str2num(nR(2:end-2));
%     nY=char(regexp(filename,'[_]\w[^_]*[P][_]','match'));nY=str2num(nY(2:end-2));
    
    if useref==0%use the first 50k a-lines to calc ref
        fid=fopen(filename);
        REF= repmat(mean(fread(fid,[SPL min(20000,nY*nX)],'ushort'),2),1,nX);
        fclose(fid);
    else
        REF = repmat(importdata(refname)',1,nX);
    end
    STR=zeros(nZcrop,nX,1,nY);
    BLD_ed=STR;  
    %%
    if show_img
        fig_flow=figure('position',[460 400 900 460]);
    end
    
    pb = ProgressBar(nY);
    bob = uint64(0);
%     for iY=1:nY
    parfor iY=1:nY
    IMG=zeros(nZcrop,nX,nR)+1i;
        fid=fopen(filename);
        fseek(fid,bob+(SPL*nX)*nR*(iY-1)*2,'bof');
        for iR=1:nR     
            A = fread(fid,[SPL nX],'ushort');
            SpecC=interp1(K,A - REF,KES,'linenumar',0);  % Interpolation + background subtraction
            Dep=fft(SpecC,SPL,1);
            IMG(:,:,iR)= Dep(IMcropRg,:);  %complex value
        end 
        %%  Sub-pixel motion compensation
        if do_reg
            Nsub=10;Colshift=0;
            [IMG,Dx] = OCTA_F_SubPixReg(IMG,Nsub,Colshift);
            %fprintf('iY=%i\t',iY);fprintf('%0.2f,',Dx);fprintf('\n');
        end
        %%
        p_tis=mean(abs(IMG),3);
        [~,p_bld_ed] = OCTA_F_ED_Clutter_EigFeed(IMG, 2);  % complex and linear scale for ED

        Stru=20*log10(p_tis);
        Flow_ed=20*log10(p_bld_ed);

        %%
        if show_img
            figure(fig_flow)
            imagesc(Flow_ed,RgFlow); ylim([imshowrgZ(1) imshowrgZ(end)])
            colormap gray;
            drawnow
        end
        %%
        STR(:,:,iY)=Stru;
        BLD_ed(:,:,iY)=Flow_ed;
        fclose(fid);
        pb.progress;
    end
    pb.stop;
%%
        fprintf('calculating range....'); 
        if use_autoRg
            nP=numel(STR);
            RgStr  =  prctile(STR(1:37:nP),[2.5 99.99]);
            RgFlow_ed=prctile(BLD_ed(1:37:nP),[2.5 99.99]);
        end
   %% Save Dicom
        fprintf('saving DICOM....'); 
        
        STRimg=uint16(65536*(STR-RgStr(1))/(RgStr(2)-RgStr(1)));
        BLnZ_ed=uint16(65536*(BLD_ed-RgFlow_ed(1))/(RgFlow_ed(2)-RgFlow_ed(1)));
        
        if do_NM
            BLDNM_ed=double(BLnZ_ed)-double(STRimg)*0.35;
            RgFlowNM_ed=prctile(BLDNM_ed(1:36:end),[8 99.99]);
            
            BLDNMimg_ed=uint16(65536*(BLDNM_ed-RgFlowNM_ed(1))/(RgFlowNM_ed(2)-RgFlowNM_ed(1)));
            if saveDicom
                dicomwrite(BLDNMimg_ed,[outputdir,filename,'_NM_Flow_ed','.dcm']);
            end
        end
    if saveDicom
            dicomwrite(STRimg,[outputdir,filename,'-StruC.dcm']);
            dicomwrite(BLnZ_ed,[outputdir,filename,'-Flow_edC.dcm']);
    end
    fprintf('Done.\n'); 
    %% Projection
    if do_proj
        Info.nZshift=5;
        Info.nZseg=234;
        ST=double(squeeze(STRimg));
        ST(1:5,:,:)=0;
        nX=size(ST,2);nY=size(ST,3);
        Info.WinSize=70;
        Info.thrsh_ratio=0.85;
        Info.medfilt=[17 19];
        Info.gaussfilt=[17 19];
        
        StrRg=[0 65535];
        seg_top=zeros(nX,nY);
        pb = ProgressBar(nY);
        for iY=1:nY
            cimgc=imgaussfilt(medfilt2(ST(:,:,iY),[17,3]),Info.gaussfilt);
            cimgMax=max(cimgc);
            for iX=1:nX
                vol=cimgc(:,iX);
                loc = find(vol>cimgMax(iX)*Info.thrsh_ratio,1,'first');
                while isempty(loc)
                    cimgMax(iX)=cimgMax(iX)*Info.thrsh_ratio;
                    loc = find(vol>cimgMax(iX)*Info.thrsh_ratio,1,'first');
                end
                seg_top(iX,iY)=loc(1);
            end
            pb.progress;
        end
        pb.stop;
        fprintf('\n')
        seg_f=medfilt2(seg_top,Info.medfilt,'symmetric');
        seg_fs=seg_f+Info.nZshift;
        %%
           if do_showSeg
                seg_f2=seg_fs+Info.nZseg;
                for iY=1:5:nY
                    cla
                    cimg=ST(:,:,iY);
                    imagesc(cimg,StrRg);
                    hold on,plot(seg_fs(:,iY),'--r'),
                    plot(seg_f2(:,iY),'--r'),
                    pause(0.05)
                end
           end
        %% Projection
        [O,P,Q]=ndgrid(1:Info.nZseg,1:nX,1:nY);
        O=bsxfun(@plus,O,permute(round(seg_fs),[3,1,2]));
        O(O>size(ST,1))=size(ST,1);
        O(O<1)=1;
        ind=sub2ind(size(ST),O,P,Q);
        clear O P Q
        % Projection
        MaxAve=3;
        gwin=0;
        GainInt=0;
        IMF=squeeze(BLDNMimg_ed);
        IMS=double(IMF(ind))/65535;
        Rg=prctile(IMS(1:20:end),[20 99.9]);
        IMS(IMS<Rg(1))=Rg(1);IMS(IMS>Rg(2))=Rg(2);
        IMS=(IMS-Rg(1))/(Rg(2)-Rg(1));
        cmap=flip(isolum(256),1);
        [Pj_c,Pj_g,CM]=ProjMaxDepth_v3(IMS,gwin,MaxAve,[0.05 0.95],cmap,[0.05 0.75]);
        Pj_g=imadjust(Pj_g,[0.3 1],[0 1]);
        Pj_ci= my_interp(Pj_c, [1 -5 20 20 -5 1]/32);
        Pj_gi= my_interp(Pj_g, [1 -5 20 20 -5 1]/32);
        Pj_ci=rot90(Pj_c);
        Pj_gi=rot90(Pj_g);
        figure(2);imshow(Pj_ci);
        figure(3);imshow(Pj_gi);
        imwrite(Pj_ci,[outputdir,filename,'-Proj_c.tiff']);
        imwrite(Pj_gi,[outputdir,filename,'-Proj_g.tiff']);
    end
end
toc
