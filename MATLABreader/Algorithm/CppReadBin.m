function [BL,STR] = CppReadBin(filename,zcrop,datasize)
fileD = fopen([filename,'.bind']);
fileS=fopen([filename,'.bins']);
nZ=datasize(1);
nX=datasize(2);
nY=datasize(3);


D = fread(fileD, [nZ nX*nY], 'float', 0);
S=fread(fileS, [nZ nX*nY], 'float', 0);
BL=D(zcrop,:);
STR=S(zcrop,:);

BL=log1p(reshape(BL,[size(zcrop,2) nX nY]));
STR=log1p(reshape(STR,[size(zcrop,2) nX nY]));

STR=uint16(65536*(STR-min(STR(:)))/(max(STR(:))-min(STR(:))));
BL=uint16(65536*(BL-min(BL(:)))/(max(BL(:))-min(BL(:))));

end

