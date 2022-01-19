%% Load
load('C:\temp\TCSPC\ScanResults1.mat')
%% SUM
check = sum(ans(:,:,150:350),3);
figure(1)
imagesc(check)

xlabel('Pixel Number')
ylabel('Slice Number')
colorbar
caxis([10 500])
title('Raw Scan')


%% VIEW

figure(3)
for i = 10:300
    t = 301-i;
    imagesc(ans(:,:,t))
  title(t)
    waitforbuttonpress
end

%%
Line = smooth(allHists(70,400,:),5);


figure(100)
plot(squeeze(Line))
%%
figure(2)
imagesc(check2)

xlabel('Pixel Number')
ylabel('Slice Number')
colorbar
caxis([100 1000])
title('Raw Scan')

%%


%%

figure(4)
subplot(2,2,1)
imagesc(check)

xlabel('Pixel Number')
ylabel('Slice Number')
colorbar
caxis([10 500])
title('Raw Scan')

subplot(2,2,2)
imagesc(check1)

xlabel('Pixel Number')
ylabel('Slice Number')
colorbar
caxis([10 500])
title('Raw Scan')

subplot(2,2,3)
imagesc(check2)

xlabel('Pixel Number')
ylabel('Slice Number')
colorbar
caxis([10 500])
title('Raw Scan')

subplot(2,2,4)
imagesc(check4)

xlabel('Pixel Number')
ylabel('Slice Number')
colorbar
caxis([10 500])
title('Raw Scan')