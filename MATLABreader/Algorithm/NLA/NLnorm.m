function [ BLnorm ] = NLnorm( BL,STR,gaussfilt,strelsize,Nratio)
switch nargin
    case 2
        gaussfilt   =5;
        strelsize   =5;
        Nratio  =1;
    case 3
        strelsize   =5;
        Nratio  =1;
    case 4
        Nratio  =1;
end
%%normalize BL image
bg      = imopen(imgaussfilt(STR,gaussfilt),strel('disk',strelsize));
BLnorm  =BL-Nratio*bg;

end

