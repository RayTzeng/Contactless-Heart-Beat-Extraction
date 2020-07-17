function [BVP,BVP_I, PR, M, Sx, W] = CrPPG(VideoFile, FS, StartTime, Duration,  PlotTF)
% CHROM_DEHAAN The Chrominance Method from: De Haan, G., & Jeanne, V. (2013). Robust pulse rate from chrominance-based rPPG. IEEE Transactions on Biomedical Engineering, 60(10), 2878-2886. DOI: 10.1109/TBME.2013.2266196
%
%   Inputs:
%       VideoFile               = Video file path.
%       FS                      = Video framerate (fps).
%       StartTime               = Timepoint at which to start process (default = 0 seconds).
%       Duration                = Duration of the time window to process (default = 60 seconds).
%       ECGFile                 = File path to corresponding ECG data file (.mat) containing: 1) The waveform - ECGData.data, 2) The ECG sampling rate - ECGData.fs, 3) The ECG peak locations (in samples) - ECGData.peaks.
%       PPGFile                 = File path to corresponding PPG data file (.mat) containing: 1) The waveform - PPGData.data, 2) The PPG sampling rate - PPGData.fs, 3) The PPG peak locations (in samples) - PPGData.peaks.
%       PlotTF                  = Boolean to turn plotting results on or off.
%
%   Outputs:
%       BVP                     = Processed Blood Volume Pulse (BVP).
%       PR                      = Estimated Pulse Rate (PR) from processed BVP timeseries using peak in periodogram.
%       HR_ECG                  = Gold standard Heart Rate (HR) measured from the ECG timeseries R-waves for the window.
%       PR_PPG                  = Pulse Rate (PR) measured from the PPG timeseries systolic onsets for the window.
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio (SNR) calculated based on the ECG HR frequency using a method adapted from the method by G. de Haan, TBME, 2013.
%
%   Requires - Signal Processing Toolbox

%% Parameters
SkinSegmentTF=true;

LPF = 0.7; %low cutoff frequency (Hz) - specified as 40 bpm (~0.667 Hz) in reference
HPF = 2.5; %high cutoff frequency (Hz) - specified as 240 bpm (~4.0 Hz) in reference

WinSec=1; %(was a 32 frame window with 20 fps camera)
PRres=1; %resolution of PR default (1 data/5 second)
WSZPR=30*FS; %size of window of PR
cnt=0;

disRange=100;
scaleRange=0.6;

%% Add Backup Functions
if(~license('test', 'image_toolbox')&&SkinSegmentTF)
    addpath([cd '\optional\rgb2ycbcr.m']);%GNU GPL rgb2ycbcr.m function
end

%% Plot Control
if(PlotTF)
    PlotPRPSD = false;
    PlotSNR = true;
else
    PlotPRPSD = false;
    PlotSNR = false;
end

%% Load Video:
VidObj = VideoReader(VideoFile);
VidObj.CurrentTime = StartTime;

FramesToRead=floor(Duration*VidObj.FrameRate); %video may be encoded at slightly different frame rate

%% Read Video and Spatially Average:
T = zeros(FramesToRead,1);%initialize time vector
RGB = zeros(FramesToRead,3);%initialize color signal
POS= zeros(FramesToRead,2); %initialize position signal
FN = 0;

%% Processing

timer1=tic;

while hasFrame(VidObj) && (VidObj.CurrentTime <= StartTime+Duration)
    FN = FN+1;
    T(FN) = VidObj.CurrentTime;
    VidFrame = readFrame(VidObj);
    %VidFrame = imresize(VidFrame, 0.5);
    
    
    %position for optional face detection/tracking - originally specified in reference as using Viola Jones, 2004
    
    [img,face,imgbbox,x,y] = cropfacebbox(VidFrame);
      
    
    VidROI=img;
    
    
    if FN==1
        prevVidROI=VidROI;
        prevImgbbox=imgbbox;
        prevX=x;
        prevY=y;
        
    else
        if face==1
            dx=x-prevX;
            dy=y-prevY;
            dis=sqrt(dx^2+dy^2);
            %fprintf('%d %d %d %d\n',prevbboxPolygon(1),bboxPolygon(1),prevbboxPolygon(2),bboxPolygon(2));
            fprintf('face=%d x=%d y=%d dis=%d\n',face,x,y,dis);
            if cnt<FS
                if ((dis>disRange)||(size(prevVidROI,2)*scaleRange>size(VidROI,2)||size(prevVidROI,1)*scaleRange>size(VidROI,1)))
                    imgbbox=prevImgbbox;
                    x=prevX;  y=prevY;
                    for i=1:size(imgbbox,1)
                        VidROI=imcrop(VidFrame,imgbbox(i,:));
                    end
                    cnt=cnt+1;
                else
                    if(cnt>0)
                        cnt=cnt-1;
                    end
                end    
             else
                 cnt=0;
             end
        else
            imgbbox=prevImgbbox;
            x=prevX; y=prevY;
            for i=1:size(imgbbox,1)
                VidROI=imcrop(VidFrame,imgbbox(i,:));
            end
            
        end
        
        prevVidROI=VidROI;
        prevImgbbox=imgbbox;
        prevX=x;
        prevY=y;
        
        POS(FN,1)=x; %x position
        POS(FN,2)=y; %y position
 
    end
    
    str = 's01';
    imwrite(VidROI,['frames\',int2str(FN), '.jpg']);
%         imwrite(img,['croppedfaces\', str,'\',int2str(k), '.jpg']);
%         if oldimgbbox(3)*0.6<=imgbbox(3)
%             oldimgbbox=imgbbox;
%             imgSize=size(img);
%             VidROI = imcrop(img,[imgSize(2)*0.2,imgSize(1)*0.1,imgSize(2)*0.6,imgSize(1)*0.8]);
%             imwrite(VidROI,['doublecroppedfaces\', str,'\',int2str(k), '.jpg']);
%         end
%     else
%         VidROI=preVidROI;
%   end   
%     
%     k=k+1;
    
    if(SkinSegmentTF)%skin segmentation - not specified in reference 
        YCBCR = rgb2ycbcr(VidROI);
        Yth = YCBCR(:,:,1)>80;
        CBth = (YCBCR(:,:,2)>77).*(YCBCR(:,:,2)<127);
        CRth = (YCBCR(:,:,3)>133).*(YCBCR(:,:,3)<173);
        ROISkin = VidROI.*repmat(uint8(Yth.*CBth.*CRth),[1,1,3]);
        RGB(FN,:) = squeeze(sum(sum(ROISkin,1),2)./sum(sum(logical(ROISkin),1),2));
        kRGB=find(isinf(RGB(FN,:))|isnan(RGB(FN,:)));
        if(~isempty(kRGB))
            fprintf('oh no!!!!!!!!!!!!!!!!!!!!!!!');
            VidROI = VidFrame;
            imgbbox=[0,0,size(VidFrame,1),size(VidFrame,2)];
            RGB(FN,:) = sum(sum(VidROI,2)) ./ (size(VidROI,1)*size(VidROI,2));
        end
    else
        RGB(FN,:) = sum(sum(VidROI,2)) ./ (size(VidROI,1)*size(VidROI,2));
    end
    
    fprintf('%d\n', floor(FN/FramesToRead*100));
end%endwhile video

if(~license('test', 'image_toolbox')&&SkinSegmentTF)%remove path if added
    rmpath([cd '\optional\']);
end

fprintf('RGB channel processing time:%d\n', toc(timer1));

%% CHROM:

timer1=tic;

NyquistF = 1/2*FS;
[B,A] = butter(3,[LPF/NyquistF HPF/NyquistF]);%Butterworth 3rd order filter - originally specified as an a FIR band-pass filter with cutoff frequencies 40-240 BPM

%Window parameters - overlap, add with 50% overlap
WinL = ceil(WinSec*FS);
if(mod(WinL,2))%force even window size for overlap, add of hanning windowed signals
    WinL=WinL+1;
end
NWin = floor((FN-WinL/2)/(WinL/2));
S = zeros(NWin,1);
WinS = 1;%Window Start Index
WinM = WinS+WinL/2;%Window Middle Index
WinE = WinS+WinL-1;%Window End Index

for i = 1:NWin
    TWin = T(WinS:WinE,:);
    
    RGBBase = mean(RGB(WinS:WinE,:));
    RGBNorm = bsxfun(@times,RGB(WinS:WinE,:),1./RGBBase)-1;
    
    % CHROM
    Xs = squeeze(3*RGBNorm(:,1)-2*RGBNorm(:,2));%3Rn-2Gn
    Ys = squeeze(1.5*RGBNorm(:,1)+RGBNorm(:,2)-1.5*RGBNorm(:,3));%1.5Rn+Gn-1.5Bn
    
    kx=find(isinf(Xs)|isnan(Xs));
    ky=find(isinf(Ys)|isnan(Ys));
    hx=find(isfinite(Xs));
    hy=find(isfinite(Ys));
    
    Xs(kx)=mean(Xs(hx));
    Ys(ky)=mean(Ys(hy));
    
    Xf = filtfilt(B,A,double(Xs));
    Yf = filtfilt(B,A,double(Ys));
    
    Alpha = std(Xf)./std(Yf);
    
    SWin = Xf - Alpha.*Yf;
    
    %fprintf("%d %d\n", WinL,size(SWin));
    
    SWin = hann(WinL).*SWin;
    %overlap, add Hanning windowed signals
    if(i==1)
        S = SWin;
        TX = TWin;
    else
        S(WinS:WinM-1) = S(WinS:WinM-1)+SWin(1:WinL/2);%1st half overlap
        S(WinM:WinE) = SWin(WinL/2+1:end);%2nd half
        TX(WinM:WinE) = TWin(WinL/2+1:end);
    end
    
    WinS = WinM;
    WinM = WinS+WinL/2;
    WinE = WinS+WinL-1;
end

BVP=S;
T=T(1:length(BVP));

fprintf('CrPPG processing time:%d\n', toc(timer1));

%% ICA
[M,Sx,W,BVP_I]=PulseRateICA(BVP,POS(1:length(BVP),1),POS(1:length(BVP),2));

%% Estimate Pulse Rate from periodogram
%PR = prpsd(BVP,FS,40,240,PlotPRPSD);

timer1=tic;

fprintf('%d %d %d %d %d',FramesToRead,FN,NWin,WinL,length(BVP_I));

PR=PulseRateCal(BVP_I,PRres,WSZPR,FS);

fprintf('Peroidogram processing time:%d\n', toc(timer1));

end%end function

function [M,S,W,BVP_I]=PulseRateICA(BVP,Px,Py) 
    M=zeros(length(BVP),3);

    M(:,1)=BVP;
    M(:,2)=Px;
    M(:,3)=Py;

    [W,Sx] = ica(M',3);   %JADE

    S=transpose(Sx);

    MinPx=zeros(1,3);
    for c=1:3
        cc = abs(corrcoef(S(:,c),M(:,2)))+abs(corrcoef(S(:,c),M(:,3)));
        MinPx(1,c)=cc(1,2);
    end

    BVP_I = S(:,argmin(MinPx,2));
end

function PR=PulseRateCal(BVP_I,PRres,WSZPR,FS)
    TPR=floor(length(BVP_I)/FS);
    NPR=floor((length(BVP_I)-WSZPR)/(FS*PRres));
    WSZPRres=(FS*PRres);
    PR = zeros(NPR);
    threshold=10;
 
    movAvg=0;
    mu=0.2;

    for i = 1 : NPR
        SPR=(i-1)*WSZPRres+1;
        EPR=SPR+WSZPR-1;
    
        %fprintf('\n%d %d %d',SPR,EPR,NPR);
    
        PR(i,1) = prpsd(BVP_I((SPR:EPR)),FS,40,240,false); 
    
        if i==1
            movAvg=PR(i,1);
            PR(i,2)=movAvg;
        else
            differ=PR(i,1)-PR(i-1,2);
            if abs(differ)>threshold
                beta=exp(-abs(differ-threshold));
            PR(i,2)=(1-beta)*movAvg+beta*PR(i,1);
            else
                PR(i,2)=PR(i,1);
            end
            movAvg=movAvg*(1-mu)+PR(i,2)*mu;
        end 
    end
end