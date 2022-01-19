%% FIGURES PLOTTING ETC RAII

clear
clc
close all 

[FileToRead, pathname] = uigetfile('*.mat');  % opens window to data file
cd(pathname);
PROCMAT = importdata(FileToRead);

%% SMOOTHING
% SET SMOOTH NUMBER - NOT WORKING
% ex_smooth=0; % DEFAULT SMOOTHING
% smooth_num=5; % MATLAB DEFAULT - Sets a span of how many to average
% 
% for i=1:size(PROCMAT,3)   % TRAVERSING ARRAY
%     if ex_smooth==0
%         SmoothedReduced(:,:,i)=smooth(PROCMAT(:,:,i),smooth_num);
%     elseif ex_smooth==1
%         SmoothedReduced(:,:,i)=smooth(PROCMAT(:,:,i),smooth_num);
%     end
% 
% end

% WORKING 
%%
%SmoothedReduced = smooth3(PROCMAT); 
SmoothedReduced = ans; 
%% VIDEO
vidtemp = rot90(SmoothedReduced,3);
% vidtemp(:,230,:) = vidtemp(:,228,:);
% vidtemp(:,231,:) = vidtemp(:,229,:);
% vidtemp(:,187:189,:) = vidtemp(:,190:192,:);

vid=1;
vid_data=1;
normalise=0;
    max_range=(150:300);
    norm_val=20;
    norm_limit=20;
% STARTING POINT
TDC_first_plot_frame = 520;
vid_length=300;
% SCALE
cmax=10;
cmin=2; 
pixel_hists=vidtemp;
if vid == 1
    if normalise==1
        aviobj = VideoWriter( '_vid_NORM.avi', 'Motion JPEG AVI');
    else
        aviobj = VideoWriter( '_vid.avi', 'Motion JPEG AVI');
    end
    % FRAME RATE
    aviobj.FrameRate = 5;
    open(aviobj);
end

% CREATE VIDEO
if vid_data==1
    figure
    for i=1:vid_length
        pause(0.05)
        imagesc(pixel_hists(:,:,TDC_first_plot_frame-i)); 
        colormap('jet');
        colorbar
        
        c_range=caxis;
        if c_range(2)<cmax
            c_range(2)=cmax;
        end
        caxis([cmin c_range(2)]);
        
         
     title(TDC_first_plot_frame-i)
        if vid == 1
            vid_plot = getframe;
            writeVideo(aviobj,vid_plot);
        end
    end
end
if vid == 1
    close(aviobj);
end
%% REBIN3
%16 Pixel Bin width returns a 32 pixel array on the detector axis 

% X = PROCMAT(:,1:96,:);
% X(is(X))=0;
% x_width = 16;
% y_width = 3;
% Y = convn(X,ones(x_width,y_width),'valid');
% BinnedMat = Y(1:x_width:end,1:y_width:end,:); %Z dimensions ? 32x32x3



% Select bin Width such that we create 32 pixels on the second axis, 128
% points/4 = 32x32 array in both axis 

%Binned = nansum(PROCMAT(1:XBinWidth,1:YBinWidth,1));
%% BEBIN Works 21/10/20
A=PROCMAT(:,:,:);
% 
% A(:,[164,156,219],:) = A(:,170:172,:);
% A(:,364,:) = A(:,365,:);
% A(:,243,:) = A(:,245,:);
% A(:,231,:) = A(:,232,:);
% A(:,230,:) = A(:,229,:);
% A(:,156,:) = A(:,157,:);
% A(:,219,:) = A(:,220,:);
% A(:,188,:) = A(:,189,:);
% A(:,10,:) = A(:,11,:);
% A(:,15,:) = A(:,16,:);
% A(:,333,:) = A(:,334,:);
% A(:,240,:) = A(:,241,:);
% A(:,85,:) = A(:,86,:);
% A(:,115,:) = A(:,116,:);
% A(:,40,:) = A(:,41,:);
% A(:,72,:) = A(:,73,:);
% A(:,146,:) = A(:,145,:);

A(isnan(A))=0;
binsz = 32;
[m,n,p] = size(A);
mr = m/binsz;
nr = n/binsz;
A= reshape(A,mr,binsz,nr,binsz,p);
A = sum(A,[1 3]);
A = reshape(A,binsz,binsz,p);
%BinnedMat2(:,16,:) = BinnedMat2(:,17,:);

BinnedMat2 = A;


% Trying to apply gaussian filter 
%B = imgaussfilt3(PROCMAT,1.76);

%% TCSPC TRACE VIEWER

% SELECT TCSPC TRACE OF INTEREST
%  Linetest = SmoothedReduced(86,180,:);
%  plot(squeeze(Linetest))
%  xlim([ 0 size(PROCMAT,3)])

%% MANUAL TIME-POINT VIEWER
% 
figure(3)
for s = 300:600
    imagesc(SmoothedReduced(:,:,s))
   % imagesc(PROCMAT(:,:,s))
    %caxis([1 3])
    colorbar
    title(['Time Bin ' num2str(s(:)')])
    waitforbuttonpress
end 

%% MULTI-TRACE VIEWER

% Multi(1,:) = smooth(SmoothedReduced(270,200,:),5);
% Multi(2,:) = smooth(SmoothedReduced(185,330,:),5);
% 
% MultiLine = Multi'; 
% 
% Dim = size (Multi);
% %Multi2 = permute (Multi', [3,2,1]);
% %Multi3 = reshape (Multi2(:),Dim(3),[]);
% 
% [row, col] = find(ismember(Multi, max(Multi(1,:,:))));
% [row2, col2] = find(ismember(Multi, max(Multi(2,:,:))));
% 
% figure(4)
% plot (MultiLine(:,:)); 
% title('TCSPC Traces')
% xlabel('Time Bins')
% ylabel('Counts')
% xline(col(1),'--b',{ num2str(col(:))},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
% xline(col2(1),'--r',{ num2str(col2(1))},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
% legend('Bull','Calculator')
% set(gca,'FontSize',15)
% grid on 
% 
% %xline(4.5,'-',{'Acceptable','Limit'});
% 
% PeakDelayPS = col-col2;

%% Peak Location

%[row, col] = find(ismember(Multi, max(Multi(1,:,:))));
%[row2, col2] = find(ismember(Multi, max(Multi(2,:,:))));

  
%% Non-time resolved sum

 %CountsMF = sum(pixel_data_cube(:,:,40:70), 'all');


% Sum(:,1:319) = sum(PROCMAT(:,1:319,240:259),3);
% Sum(:,381:512) = sum(PROCMAT(:,381:512,240:259),3);
% Sum(:,320:380) = sum(PROCMAT(:,320:380,200:240),3);


 Sum = sum(PROCMAT(:,:,10:300),3);

 


% 
% SumB = sum(BinnedMat2(:,:,10:300),3);
SumSection = sum(PROCMAT(:,1:300,250:258),3);

% CLEANING
Sum(:,[164,156,219],:) = Sum(:,170:172,:);
Sum(:,364,:) = Sum(:,365,:);
Sum(:,243,:) = Sum(:,245,:);
Sum(:,231,:) = Sum(:,232,:);
Sum(:,230,:) = Sum(:,229,:);
Sum(:,156,:) = Sum(:,157,:);
Sum(:,219,:) = Sum(:,220,:);
Sum(:,188,:) = Sum(:,189,:);
Sum(:,10,:) = Sum(:,11,:);
Sum(:,15,:) = Sum(:,16,:);
Sum(:,333,:) = Sum(:,334,:);
Sum(:,240,:) = Sum(:,241,:);
Sum(:,85,:) = Sum(:,86,:);
Sum(:,115,:) = Sum(:,116,:);
Sum(:,40,:) = Sum(:,41,:);
Sum(:,72,:) = Sum(:,73,:);
Sum(:,146,:) = Sum(:,145,:);

%Sum(1:48,1:98,:) = Sum(1:48,400:497,:);
%Sum(1:48,108:200,:) = Sum(1:48,408:500,:);
% 
%  SumB(:,6,:) = SumB(:,5,:);
%  SumB(:,1,:) = SumB(:,2,:);
%  SumB(:,3,:) = SumB(:,4,:);
%  SumB(1:14,12,:) = SumB(1:14,13,:);
%  SumB(20:32,12,:) = SumB(20:32,13,:);
%  SumB(1:14,16,:) = SumB(1:14,17,:);
%  SumB(20:32,16,:) = SumB(20:32,17,:);
 
Rotate = rot90(Sum,3);

%%

figure(5)
%totaltest(totaltest>1000) = NaN;

%subplot(2,1,1)
SumAxis = nansum(Sum(:,:), 'all');
imagesc(Sum)
caxis([2 100])
colorbar
title('RAII')
xlabel('Pixel Number')
ylabel('Slice Number')
set(gca,'FontSize',15)

% subplot(2,1,2)
 %%
 
 figure(6)
imagesc(SumB)
caxis([0 10000])
colorbar
title('RAII')
xlabel('Pixel Number')
ylabel('Slice Number')
set(gca,'FontSize',15)

%% Multiple Time Windows

% Sum1 = sum(PROCMAT(:,:,280:290),3);
% Sum2 = sum(PROCMAT(:,:,270:280),3);
% Sum3 = sum(PROCMAT(:,:,250:270),3);
% Sum4 = sum(PROCMAT(:,:,240:260),3);
% 
% 
% 
% % CLEANING
% Sum1(:,[164,156,219],:) = Sum1(:,170:172,:);
% Sum1(:,364,:) = Sum1(:,365,:);
% Sum1(:,243,:) = Sum1(:,245,:);
% Sum1(:,231,:) = Sum1(:,232,:);
% Sum1(:,230,:) = Sum1(:,229,:);
% Sum1(:,156,:) = Sum1(:,157,:);
% Sum1(:,219,:) = Sum1(:,220,:);
% Sum1(:,188,:) = Sum1(:,189,:);
% Sum1(:,10,:) = Sum1(:,11,:);
% Sum1(:,15,:) = Sum1(:,16,:);
% Sum1(:,333,:) = Sum1(:,334,:);
% Sum1(:,240,:) = Sum1(:,241,:);
% Sum1(:,85,:) = Sum1(:,86,:);
% Sum1(:,115,:) = Sum1(:,116,:);
% Sum1(:,40,:) = Sum1(:,41,:);
% Sum1(:,72,:) = Sum1(:,73,:);
% Sum1(:,[164,156,219],:) = Sum1(:,170:172,:);
% Sum1(:,364,:) = Sum1(:,365,:);
% Sum1(:,243,:) = Sum1(:,245,:);
% Sum1(:,231,:) = Sum1(:,232,:);
% Sum1(:,230,:) = Sum1(:,229,:);
% Sum1(:,156,:) = Sum1(:,157,:);
% Sum1(:,219,:) = Sum1(:,220,:);
% Sum1(:,188,:) = Sum1(:,189,:);
% Sum1(:,10,:) = Sum1(:,11,:);
% Sum1(:,15,:) = Sum1(:,16,:);
% Sum1(:,333,:) = Sum1(:,334,:);
% Sum1(:,240,:) = Sum1(:,241,:);
% Sum1(:,85,:) = Sum1(:,86,:);
% Sum1(:,115,:) = Sum1(:,116,:);
% Sum1(:,40,:) = Sum1(:,41,:);
% Sum1(:,72,:) = Sum1(:,73,:);
% 
% Sum2(:,[164,156,219],:) = Sum2(:,170:172,:);
% Sum2(:,364,:) = Sum2(:,365,:);
% Sum2(:,243,:) = Sum2(:,245,:);
% Sum2(:,231,:) = Sum2(:,232,:);
% Sum2(:,230,:) = Sum2(:,229,:);
% Sum2(:,156,:) = Sum2(:,157,:);
% Sum2(:,219,:) = Sum2(:,220,:);
% Sum2(:,188,:) = Sum2(:,189,:);
% Sum2(:,10,:) = Sum2(:,11,:);
% Sum2(:,15,:) = Sum2(:,16,:);
% Sum2(:,333,:) = Sum2(:,334,:);
% Sum2(:,240,:) = Sum2(:,241,:);
% Sum2(:,85,:) = Sum2(:,86,:);
% Sum2(:,115,:) = Sum2(:,116,:);
% Sum2(:,40,:) = Sum2(:,41,:);
% Sum2(:,72,:) = Sum2(:,73,:);
% 
% Sum3(:,[164,156,219],:) = Sum3(:,170:172,:);
% Sum3(:,364,:) = Sum3(:,365,:);
% Sum3(:,243,:) = Sum3(:,245,:);
% Sum3(:,231,:) = Sum3(:,232,:);
% Sum3(:,230,:) = Sum3(:,229,:);
% Sum3(:,156,:) = Sum3(:,157,:);
% Sum3(:,219,:) = Sum3(:,220,:);
% Sum3(:,188,:) = Sum3(:,189,:);
% Sum3(:,10,:) = Sum3(:,11,:);
% Sum3(:,15,:) = Sum3(:,16,:);
% Sum3(:,333,:) = Sum3(:,334,:);
% Sum3(:,240,:) = Sum3(:,241,:);
% Sum3(:,85,:) = Sum3(:,86,:);
% Sum3(:,115,:) = Sum3(:,116,:);
% Sum3(:,40,:) = Sum3(:,41,:);
% Sum3(:,72,:) = Sum3(:,73,:);
% 
% Sum4(:,[164,156,219],:) = Sum4(:,170:172,:);
% Sum4(:,364,:) = Sum4(:,365,:);
% Sum4(:,243,:) = Sum4(:,245,:);
% Sum4(:,231,:) = Sum4(:,232,:);
% Sum4(:,230,:) = Sum4(:,229,:);
% Sum4(:,156,:) = Sum4(:,157,:);
% Sum4(:,219,:) = Sum4(:,220,:);
% Sum4(:,188,:) = Sum4(:,189,:);
% Sum4(:,10,:) = Sum4(:,11,:);
% Sum4(:,15,:) = Sum4(:,16,:);
% Sum4(:,333,:) = Sum4(:,334,:);
% Sum4(:,240,:) = Sum4(:,241,:);
% Sum4(:,85,:) = Sum4(:,86,:);
% Sum4(:,115,:) = Sum4(:,116,:);
% Sum4(:,40,:) = Sum4(:,41,:);
% Sum4(:,72,:) = Sum4(:,73,:);
% 
% Rotate1 = rot90(Sum1,3);
% Rotate2 = rot90(Sum2,3);
% Rotate3 = rot90(Sum3,3);
% Rotate4 = rot90(Sum4,3);



%% Multi Sum Comp

% figure(7)
%  subplot(2,2,1)
% imagesc(Rotate1)
% caxis([50 300])
% colorbar
% title('Non Time-Resolved Sum')
% xlabel('Pixel Number')
% ylabel('Slice Number')
% set(gca,'FontSize',15)
% 
%  subplot(2,2,2)
% imagesc(Rotate2)
% caxis([5 300])
% colorbar
% title('Non Time-Resolved Sum')
% xlabel('Pixel Number')
% ylabel('Slice Number')
% set(gca,'FontSize',15)
% % 
% subplot(2,2,3)
% imagesc(Rotate3)
% caxis([5 200])
% colorbar
% title('Non Time-Resolved Sum')
% xlabel('Pixel Number')
% ylabel('Slice Number')
% set(gca,'FontSize',15)
% % 
% subplot(2,2,4)
% imagesc(Rotate4)
% caxis([5 100])
% colorbar
% title('Non Time-Resolved Sum')
% xlabel('Pixel Number')
% ylabel('Slice Number')
% set(gca,'FontSize',15)


%% Filtering and Additional Plotting

 %G = imgaussfilt(PROCMAT(:,:,250:257),1);

%% ROI

% figure(11)
% sgtitle('Region of Interest')
% set(gca,'FontSize',15)
% 
% subplot(1,2,1)
% imagesc(Sum(205:230,260:285))
% title('Cental Point')
% caxis([1 200])
% 
% subplot(1,2,2)
% imagesc(Sum(40:65,195:220))
% caxis([1 200])
% title('Edge Point')
% colorbar

%% Counting 
% 
% zCountsRA = sum(PROCMAT(8:17,10:19,:), 'all');
% zCountsTotal = nansum(PROCMAT(:,:,:), 'all');
% 
% zCountsBinnedTotal = sum(BinnedMat2(:,:,:), 'all');
% 
% zCBinnedPP = zCountsBinnedTotal ./ 1024;
% zCBinnedPPPS = zCBinnedPP ./ 5.12e-3;
% 
% zCRA = zCountsRA ./ 4400;

CountsTotalAll = nansum(PROCMAT(:,:,10:480),'all');

CountsWindow = nansum(PROCMAT(:,:,150:299),'all');
%CountsWindowPP = CountsWindow ./ 262144;

CountsTotalBinned = nansum(BinnedMat2(:,:,10:480), 'all');
CountsTotalBinnedPP = CountsTotalBinned ./ 1024;

CountsWindowBinned = nansum(BinnedMat2(:,:,150:299), 'all');
CountsWindowBinnedPP = CountsWindowBinned ./ 1024;

% ClippedRAII = PROCMAT(:,:,11:1000);
% 
% ClippedRAIIWindow = PROCMAT(:,:,150:299);

% SumPosRAII = sum(sum(ClippedRAII(ClippedRAII >0)));
% SumPosRAIIWind = sum(sum(ClippedRAIIWindow(ClippedRAIIWindow >0)));

%% Single Slice Processing

CountsTotalAll = nansum(PROCMAT(2,:,10:480),'all');

CountsWindow = nansum(PROCMAT(2,:,150:299),'all');


%% TCSPC TRACES
% 
 Linetest2 = BinnedMat2(14,18,:);
 figure(13)
 Line = squeeze(Linetest2);
 
 %LineConversion = Line ./ 0.031240;
  LineConversion = Line .* 2;
 
  
 
% plot(squeeze(smooth(Linetest2,10)),'Displayname','RAII')
 hold on
 % plot(squeeze(Line))
%  plot(squeeze(mfline),'Displayname','MegaFrame')
   plot(squeeze(smooth(LineConversion)),'Displayname','RAII')
 %  plot(CompLineRA)
 xlim([ 0 size(BinnedMat2,3)])
 ylim([-5 6000])
  title('TCSPC Traces')
xlabel('Time Bins')
ylabel('Counts')
 grid on
 legend('show')
 
 %% TRACES
 
 ProcessedLine = BinnedMat2(14,18,:);
% RawLine = allHists(300,300,:);

 ProLine = smooth(squeeze(ProcessedLine));
% RLine = squeeze(RawLine);
 
 %LineConversion = Line ./ 0.031240;
  LineConversion = Line .* 1;
 
  figure(14)
% plot(squeeze(smooth(Linetest2,10)),'Displayname','RAII')
 hold on
 % plot(squeeze(Line))
%  plot(squeeze(mfline),'Displayname','MegaFrame')
   plot(ProLine,'Displayname','RAII')
   %   plot(RLine,'Displayname','Raw')
 xlim([1 600])
 ylim([-10 500])
  title('TCSPC Traces')
xlabel('Time Bins')
ylabel('Counts')
 grid on
 legend('show')
%% PSF THRESHOLDING
% 
% SumPSF = sum(PROCMAT(:,:,240:260),3);
% 
% for i = 40:size(SumPSF,1)
%     
% for j = 1:200 
% 
% if SumPSF(i,j) >15 && SumPSF(i,j) <300
% 
%     SumPSF(i,j) = SumPSF(i,j) .* 4;
% 
% elseif SumPSF(i,j) == SumPSF(i,j)
%     
% end
% end
% end
% 
% for i = 50:size(SumPSF,1)
%     
% for j = 400:512
% 
% if SumPSF(i,j) >10 && SumPSF(i,j) <300
% 
%     SumPSF(i,j) = SumPSF(i,j) .* 3;
% 
% elseif SumPSF(i,j) == SumPSF(i,j)
%     
% end
% end
% end
% 
% SumPSF(SumPSF(:,:)<5) = 1;
% 
% 
% figure(19)
% imagesc(SumPSF)
% caxis([0 500])
% colorbar
% title('Non Time-Resolved Sum')
% xlabel('Pixel Number')
% ylabel('Scan Position')
% set(gca,'FontSize',15)

%%

% figure(20)
% %surf(Sum(125:145,390:410))
% %surf(Sum(125:145,400:420))
% caxis([0 100])
% xlim([1 21])
% ylim([1 21])
% colorbar
% title('Single Point (Edge)')
% xlabel('Pixels')
% ylabel('Pixels')
% set(gca,'FontSize',15)
% 
% sgtitle('Focus Variation with Lens Position')
% set(gca,'FontSize',15)
% 
% subplot(2,2,1)
% surf(SumBlur(205:230,260:285))
% title('Out of Focus')
% colorbar
% xlim([1 25])
% ylim([1 25])
% caxis([1 150])
% view(45,55)
% 
% subplot(2,2,2)
% surf(SumBlur(205:230,260:285))
% title('Out of Focus')
% colorbar
% xlim([1 25])
% ylim([1 25])
% caxis([1 150])
% view(2)
% 
% subplot(2,2,3)
% surf(Sum(205:230,260:285))
% caxis([1 150])
% xlim([1 25])
% ylim([1 25])
% title('In Focus')
% colorbar
% view(45,55)
% 
% subplot(2,2,4)
% surf(Sum(205:230,260:285))
% caxis([1 150])
% xlim([1 25])
% ylim([1 25])
% title('In Focus')
% colorbar
% view(2)

% 


%%
% %% Multi PSF
% 
% % figure(21)
% % subplot(2,2,1)
% % surf(Sum(250:262,214:228,:))
% % caxis([0 50])
% % colorbar
% % title('Most Central')
% % 
% % subplot(2,2,2)
% % surf(Sum(180:192,214:228,:))
% % caxis([0 50])
% % colorbar
% % title('Two Steps from Central')
% % 
% % subplot(2,2,3)
% % surf(Sum(110:122,214:228,:))
% % caxis([0 50])
% % colorbar
% % title('Four Steps from Central')
% % 
% % subplot(2,2,4)
% % surf(Sum(8:20,214:228,:))
% % caxis([0 50])
% % colorbar
% % title('Edge Spot')
% 
% %%
% 
% % figure(99)
% % 
% % plot(Sum(260,:));
% 
% %%
% 

 r = Sum;
 %r = Sum(:,350:512); 
 z1 = squeeze(sum(r,1));
 
 
 z = squeeze(sum(r,2));
  % z = squeeze(sum(r(1:320,:),2));
 
%  figure(19)
%  plot(z);
 
% 
% %%
% 

%  
%  %[pks1,locs1,widths1,proms1] = findpeaks(z1,'MinPeakDistance',20);
%  
%  [pks1,locs1,widths1,proms1] = findpeaks(z1,'MinPeakProminence',5300);
%  
%   [pks5,locs5,widths5,proms5] = findpeaks(z5,'MinPeakProminence',1800);
%   
%    [pks6,locs6,widths6,proms6] = findpeaks(z6,'MinPeakProminence',1800);
%    
%     
% % figure(15)
% % plot(z)
% % text(locs+.02,pks,num2str((1:numel(pks))'))
% 
% %% PSF Variation w/ FWHM
% 

%%
% 
z66 = z;
% 
% 
% for i = 350:395
%  z66(i) = z66(i) + i .*6;
%  
% end
% 
% for i = 396:404
%  z66(i) = z66(i) + i .*6.5;
%  
% end
% 
% for i = 405:512
%  z66(i) = z66(i) + i .*7;
%  
% end
% 
 z66 = smoothdata(z66,'movmean',2);
% 
%  [pks66,locs66,widths66,proms66] = findpeaks(z66,'MinPeakProminence',1400);
% 
%%
% [pks,locs,widths,proms] = findpeaks(z66,'MinPeakProminence',650);

% figure(21)
% 
%  findpeaks(z66,'MinPeakProminence',490,'Annotate','extents')
% % findpeaks(z(1:320,:),'MinPeakProminence',2500,'Annotate','extents')
% for l = 1:size(widths)
%     
% text(locs(l,1)+3,pks(l,1), num2str(widths(l,1)))
% 
% end 
% title('PSF FWHM across FOV')
% xlabel('Pixel Number')
% ylabel('Intensity')
% set(gca,'FontSize',20)
% xlim([20 460])
 %rectangle('Position',[230 4800 19 1000], 'FaceColor','w', 'EdgeColor', 'none')
 
 %% COMP
 
%  n0 = z0 ./ max(z0);
%  n10 = z10 ./ max(z10);
%  n20 = z20 ./ max(z20);
%  n30 = z30 ./ max(z30);
%  n40 = z40 ./ max(z40);
%  
% %  
%  figure(21)
% hold on
% 
%  %plot(n0,'LineWidth',1.5, 'DisplayName', 'SM1Z + 0')
% 
%  %plot(n10 + 0.1,'LineWidth',1.5, 'DisplayName', 'SM1Z + 10')
%  plot(n20 ,'LineWidth',1.5, 'DisplayName', 'SM1Z + 20')
%  plot(n30 + 0.05 ,'LineWidth',1.5, 'DisplayName', 'SM1Z + 30')
%  plot(n40 + 0.1,'LineWidth',1.5, 'DisplayName', 'SM1Z + 40')
%  
%  hold off
% 
% title('PSF Variation with SM1Z Position (70cm Working Distance)')
% xlabel('Pixel Number')
% ylabel('Intensity')
% set(gca,'FontSize',20)
% xlim([20 460])
% 
% legend('show')
%  %rectangle('Position',[230 4800 19 1000], 'FaceColor','w', 'EdgeColor', 'none')
 %%
% 
% 
% hold on 
% 
% %  findpeaks(z5,'MinPeakProminence',1800,'Annotate','extents')
% % % findpeaks(z(1:320,:),'MinPeakProminence',2500,'Annotate','extents')
% % for l = 1:size(w5)
% %     
% % text(locs5(l,1),pks5(l,1)+500, num2str(w5(l,1)))
% % 
% % end 
% % 
% %  findpeaks(z4,'MinPeakProminence',1390,'Annotate','extents')
% % % findpeaks(z(1:320,:),'MinPeakProminence',2500,'Annotate','extents')
% % for l = 1:size(w4)
% %     
% % text(locs(l,1),pks(l,1) + 250, num2str(w4(l,1)))
% % 
% % end 
% 
% %  findpeaks(z7,'MinPeakProminence',2500,'Annotate','extents')
% % % findpeaks(z(1:320,:),'MinPeakProminence',2500,'Annotate','extents')
% % for l = 1:size(w7)
% %     
% % text(locs(l,1)+3,pks(l,1), num2str(w7(l,1)))
% % 
% % end 
% 
%  findpeaks(z66,'MinPeakProminence',1400,'Annotate','extents')
% % findpeaks(z(1:320,:),'MinPeakProminence',2500,'Annotate','extents')
% for l = 1:size(widths66)
%     
% text(locs66(l,1),pks66(l,1)+500, num2str(widths66(l,1)),'FontSize', 15)
% 
% end 
% 
% hold off 
% xlim([0 460])
% rectangle('Position',[155 4100 40 960], 'FaceColor','w', 'EdgeColor', 'none')
% rectangle('Position',[320 6500 28 450], 'FaceColor','w', 'EdgeColor', 'none')
% 
% title('PSF FWHM across FOV')
% xlabel('Pixel Number')
% ylabel('Intensity')
% set(gca,'FontSize',20)
% 
% %%
% figure(21)
% 
% findpeaks(z1,'MinPeakProminence',1000,'Annotate','extents')
% 
% for l = 1:size(widths1')
%     
% text(locs1(1,l)+3,pks1(1,l), num2str(widths1(1,l)))
% 
% end 
% % 
% 
% %%
% 
% n6 = smoothdata(z6,'movmean',2) ./ max(z6);
% 
% n6(1:140) = n6(1:140) + 0.32;
% 
% n5 = smoothdata(z5,'movmean',2) ./ max(z5);
% 
% n5(1:140) = n5(1:140) + 0.32;
% 
% n4 = smoothdata(z4,'movmean',2) ./ max(z4);
% 
% n4(1:140) = n4(1:140) + 0.32;
% 
% %%
% 
% % n6 = z6 ./ max(z6);
% % 
% % n5 = z5 ./ max(z5);
% % 
% % n4 = z4 ./ max(z4);
% 
% %%
% 
% figure(23)
% 
% plot(n6 +0.032,'LineWidth',1.5, 'DisplayName', 'SM1Z Cental + 1mm - Scan 6')
% 
% hold on
% 
% %plot(z7,'LineWidth',2, 'DisplayName', 'Scan7')
% 
%  plot(n5 + 0.205,'LineWidth',1.5, 'DisplayName', 'Central SM1Z - Scan 5')
%  plot(n4 + 0.322,'LineWidth',1.5, 'DisplayName', 'SM1Z Cental - 1mm - Scan 4')
% hold off 
% 
% xlim([0 460])
% 
% title('PSF FWHM variation with Lens Position')
% xlabel('Pixel Number')
% ylabel('Normalised Intensity (with Offset)')
% set(gca,'FontSize',20)
% %annotation('textbox', [0.8, 0.7, 0.100001, 0.1], 'String', "Average FWHM")
% legend('show')
% 
% 
% 
% 
% %%
% 
% figure(22)
% 
% plot(z6,'LineWidth',2, 'DisplayName', 'Scan6')
% 
% hold on
% 
% %plot(z7,'LineWidth',2, 'DisplayName', 'Scan7')
% plot(z5,'LineWidth',2, 'DisplayName', 'Scan5')
% plot(z4,'LineWidth',2, 'DisplayName', 'Scan4')
% hold off 
% 
% legend('show')
% 
% 
% %%
% figure(50)
% grid off
%  findpeaks(smoothdata(z,'movmean',2),'MinPeakProminence',1500,'Annotate','extents')
% % findpeaks(z(1:320,:),'MinPeakProminence',2500,'Annotate','extents')
% for l = 1:size(widths)
%     
% text(locs(l,1),pks(l,1)+500, num2str(widths(l,1)),'FontSize', 15)
% 
% end 
% 
% hold off 
% 
% rectangle('Position',[255 1000 40 4840], 'FaceColor','w', 'EdgeColor', 'none')
% %rectangle('Position',[320 6500 28 450], 'FaceColor','w', 'EdgeColor', 'none')
% grid on
% title('FWHM of PSF Across FOV')
% xlabel('Pixel Number')
% ylabel('Intensity')
% set(gca,'FontSize',20)

%%


%Clipped = pixel_data_cube(:,:,11:1000);

%SumPostive = sum(sum(Clipped(Clipped >0)));

% ClippedRAII = PROCMAT(:,:,11:1000);
% 
% ClippedRAIIWindow = PROCMAT(:,:,150:299);
% 
% SumPosRAII = sum(sum(ClippedRAII(ClippedRAII >0)));
% SumPosRAIIWind = sum(sum(ClippedRAIIWindow(ClippedRAIIWindow >0)));

%SumPosRAIICOMPENSATION = SumPosRAII .*2;

%%

%SumMF = sum(pixel_data_cube(:,:,10:280),3);
%SumB = sum(BinnedMat2(:,:,10:280),3);
SumR = sum(PROCMAT(:,:,10:280),3);

figure(20)
sgtitle('Sponge Scene Imaging (200us Exposure)')
 set(gca,'FontSize',15)
subplot(2,2,1)
imagesc(SumM)
caxis([ 0 50000])
colorbar
title('MF')

subplot(2,2,2)
imagesc(SumB)
caxis([ 0 50000])
colorbar
title('RAII Binned')

subplot(2,2,3)
imagesc(SumR)
caxis([ 0 1000])
colorbar
title('RAII Full Res')


%%

LineMid = PROCMAT(:,256,:);
SumMID = sum(LineMid(:,:,10:280),'all');


%% 20-11-20 Processing

CompLineRA = squeeze(BinnedMat2(:,16,:));

SumCompLineRA = CompLineRA(:,20:299);

SumCompLineRA(SumCompLineRA < 0) = NaN;

CompLineRASum = nansum(SumCompLineRA(:,:),2);


%CompLineRASumWindow = nansum(SumCompLineRA(150:280,:));

%% Pixel Variation 20-11

%CountsLine = BinnedMat2(:,20,:);

%CountsLine2 = sum(CountsLine(:,:,20:300),3);

CountsLine512 = PROCMAT(288,:,:);

CountsLine512(CountsLine512 < 0) = NaN;

CountsLine5122 = nansum(CountsLine512(:,:,10:300),3);

CountsLine5122Window = nansum(CountsLine512(:,:,150:300),3);

CountsLine512Sorted = CountsLine5122';

SumCountsLine = sum(CountsLine512Sorted);


SumCountsLineWindow = sum(CountsLine5122Window);
% BinnedLine = BinnedMat2(16,:,:);
% 
% BinnedLine(BinnedLine < 0) = NaN;
% 
% BinnedLineCounts = nansum(BinnedLine(:,:,20:300),3);
% 
% BinnedLineCountsSorted = BinnedLineCounts';


%plot(CountsLine7)
%%

SumM = sum(pixel_data_cube(:,:,10:300),3);

figure(110)
imagesc(SumM)
caxis([ 0 2000])


%%

  figure(200)
  
  BrightMF = squeeze(pixel_data_cube(11,17,1:600));
  DimMF = squeeze(pixel_data_cube(25,17,1:600));
  
  BrightRA = squeeze(smooth(BinnedMat2(17,20,1:600)));
  DimRA = squeeze(smooth(BinnedMat2(27,19,1:600)));

%subplot(1,2,1)
hold on
plot(BrightRA, 'displayname', 'Bright RA')
plot(BrightMF, 'displayname', 'Bright MF')
 xlim([1 600])
 ylim([-10 500])
  title('Bright Pixel Comparison')
xlabel('Time Bins')
ylabel('Counts')
 grid on
 legend('show')
 set(gca,'FontSize',15)

% subplot(1,2,2)
% hold on
% plot(DimRA, 'displayname', 'Dim RA')
% plot(DimMF, 'displayname', 'Dim MF')
% title('Bright')
% title('RAII Binned')
% legend('show')
% ylim([ -10 300])
%   


  
  %%
  
%   POSBrightRA = BrightRA;
%   POSBrightRA(POSBrightRA < 0) = NaN;
%   
%     POSBrightMF = BrightMF;
%   POSBrightMF(POSBrightMF < 0) = NaN;
%   
%   %%
%   
%   
%   ySumRAPIX = nansum(POSBrightRA(10:300,:),'all');
%   ySumMFPIX = nansum(POSBrightMF(10:300,:),'all');
%   
%     ySumRAPIXWIND = nansum(POSBrightRA(150:300,:),'all');
%   ySumMFPIXWIND = nansum(POSBrightMF(100:250,:),'all');
%   
  
  %%
%   CountsLine512 = PROCMAT(288,:,:);
% CountsLine512(CountsLine512 < 0) = NaN;
% CountsLine5122 = nansum(CountsLine512(:,:,10:300),3);
% CountsLine5122Window = nansum(CountsLine512(:,:,150:300),3);
% CountsLine512Sorted = CountsLine5122';
% SumCountsLine = sum(CountsLine512Sorted);
% SumCountsLineWindow = sum(CountsLine5122Window);
% 
% 
% %%
% CountsLine32 = squeeze(BinnedMat2(18,:,:));
% CountsLine32(CountsLine32 < 0) = NaN;
% CountsLine322 = nansum(CountsLine32(:,10:300),2);
% CountsLine32Window = nansum(CountsLine32(:,150:300),2);
% CountsLine32Sorted = CountsLine322';
% SumCountsLine32 = sum(CountsLine32Sorted);
% SumCountsLineWindow32 = sum(CountsLine32Window);
% 
% figure(137)
% %plot(CountsLine2)
% hold on
% plot(CountsLine322)
% 
%   %%
%   
%   MFCountsLine = squeeze(pixel_data_cube(:,17,:));
%   MFCountsLine(MFCountsLine < 0) = NaN;
%   
% MFCountsLine2 = nansum(MFCountsLine(:,10:300),2);
% MFCountsLine2Window = nansum(MFCountsLine(:,150:300),2);
% 
% MFSumCountsLine = sum(MFCountsLine2);
% MFSumCountsLineWindow = sum(MFCountsLine2Window);
  
%% Single Slice

SingleSlice = PROCMAT(2,:,:);

SumSingleSlice = sum(SingleSlice(:,:,10:300),3);
SumSingleSliceWindow = sum(SingleSlice(:,:,150:300),3);

Slice = SumSingleSlice';
SliceWindow = SumSingleSliceWindow';

SumSlice = nansum(Slice);
SumSliceW = nansum(SliceWindow);

figure(300)
hold on 

plot(Slice)
plot(old5)

legend('show')


%%

figure(11)

for i = 20:299
    t = 300-i;
    imagesc(ans(:,:,t))
    %imagesc(allHists(:,:,t))
    title(t)
    waitforbuttonpress
end 
