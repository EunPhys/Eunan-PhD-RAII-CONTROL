% RAII PROCESSING SCRIPT

clear
clc
close all 

%% DATA READ IN + DARK SUB

% SCANS ARE SAVED AS .MAT FILES - LOAD IN AND GATHER RAW MATRIX
% ShowStack(MAT) to visualise and flick through slides! w/ Gaussian Filter

[FileToRead, pathname] = uigetfile('*.mat');  % opens window to data file
cd(pathname);
Scan = importdata(FileToRead);

% RAW MATRIX PIXELSxPOSITIONxTIMEBIN, ALTERNATE 1,2 TO CHANGE ORIENTATION
rawmatrix = permute(Scan,[1 2 3]);

%Bins = 1:600;

rawmatrix = rawmatrix(:,:,1:1000);


% DARK SUBTRACTION
DCR_RANGE = 280:300;
DCR2D = nansum(rawmatrix(:,:,DCR_RANGE), 3) ./size(DCR_RANGE,2);
DCR3D = repmat(DCR2D, 1, 1, 1000);

% DARK REDUCED MATRIX
MINUSDARK = rawmatrix - DCR3D;

%% PROCESSING THRESHOLD + FIG

threshold1=max(DCR2D(:,1), [], 'all')/1.1;

%  figure(1)
%  scatter(1:size(DCR2D,2),DCR2D(1,:),'Displayname','DCR');
%  hold on;
%  %plot(ones(1,size(DCR2D,2))*threshold1,'Displayname','Threshold 1')
%  %plot(ones(1,size(DCR2D,2))*threshold2,'Displayname','Threshold 2')
%  plot(ones(1,size(DCR2D,2))*threshold1,'linewidth',2,'Displayname','Threshold')
%  hold off
%  legend('show')
%  title('Mean DCR/Pix (Range 400:450)')
%  xlim([0 520])
%  set(gca,'FontSize',20)

%% SCREAMER REMOVAL
screamers=find((DCR2D(1,:))>threshold1);

no_scream=MINUSDARK;

% LETS REPLACE THE SCREAMER COLUMNS
for t=1:size(MINUSDARK,2)
    for i=screamers
        no_scream(:,i,:) = NaN(1);
    end
    
    fill = no_scream;
    % FILL IN THE BLANKS
 
    for i=screamers
        if ((487)>i) && (i>25)
            fill(:,i,:) = ((no_scream(:,i+1,:) + no_scream(:,i-1,:) ./2));
      
        elseif (i>1) && (i<=25)
            fill(:,i,:) = ((no_scream(:,i+1,:) + no_scream(:,i-1,:) + no_scream(:,i+25,:) ./3));
            
        elseif (i<512) && (i>=487)
            fill(:,i,:) = ((no_scream(:,i+1,:) + no_scream(:,i-1,:) + no_scream(:,i-25,:)./3));
        end       
    end
end 

% REDUCED ARRAY
PROCMAT = fill; 

%% CLIP

%PROCMAT = PROCMAT(1:size(PROCMAT,1),:,:);
%%

save('/Users/eunanmcshane/Desktop/Camera Comparison/RAII OCTOBER-NOV/RAII/13-11-20/Processed Matrices/Comp/Sponge/Single Slice TCSPC/ProcessedScan12.mat','PROCMAT');

%%

%  Linetest2 = rawmatrix(300,300,:);
%  figure(13)
%  Line = squeeze(smooth(Linetest2));
%  plot(Line)