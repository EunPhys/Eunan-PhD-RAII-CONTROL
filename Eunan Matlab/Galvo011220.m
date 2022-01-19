%%

daq.reset
daq.HardwareInfo.getInstance('DisableReferenceClockSynchronization',true); % This is to overcome a different issue
d = daq.getDevices;
src = daq.createSession('ni');
ch1=addAnalogOutputChannel(src,'Dev1','ao0','Voltage');
ch2=addAnalogOutputChannel(src,'Dev1','ao1','Voltage');
src.Rate = 1000; %speed it execiutes the que
%outputSingleValue = 0.2;
outputSingleScan(src, [ 0 0]);

%%
resolution=200; 
r=0.1;
theta = 0:2*pi/resolution:2*pi;
for i=1:resolution
% 

points(i,1)=r*cos(theta(i));
points(i,2)=r*sin(theta(i));
outputSingleScan(src, [points(i,1) points(i,2)]);
 disp(i) 

 disp(points(i,:)); 
% disp([r theta(i)])
waitforbuttonpress
end

% for k=1:1
% 
% for i=1:length(points)
% disp(i) 
% outputSingleScan(src,points(i,:));
% pause(0.1)
% disp(points(i,:)); 
% disp([r theta(i)]); 
% pause
% end
% 
% end
% 
% pause(2)
% queueOutputData(src,[outputSignal1 outputSignal2]);
% src.startBackground;

% %%
%%
%infile = 'Sample1_220721.txt';
infile = 'Nov-5cm-15-5.txt';
%infile = 'april59points.txt';
points = load(infile, '-ascii');

for p = 20
%for p = 1:5:50 
%for p = 1
  % pause(1)
outputSingleScan(src, [points(p,1) points(p,2)]);
 %fprintf('%d ', p);
 disp(p)
 waitforbuttonpress
 %pause(1)
end 

%%
% Find the rough point =  outputSingleScan(src, [0 0]);

%These are X and Y Values from your rough point
GridX = -0.5136; 
GridY = -0.3984;

%Max Size in X,Y you want to scan across
GridSize = 0.1;
% Res - Needs to be incremental
Resolution = GridSize/20;

% Outer Y loop going from bottom left of grid to top right
for y = (GridY-GridSize/2):Resolution:(GridY+GridSize/2)
        fprintf('\n')
        fprintf('y = %d ',y)
        fprintf('\n')
    
    %Move horizontally for each Y value
    for x = (GridX-GridSize/2):Resolution:(GridX+GridSize/2)
   
         fprintf('x = %d ',x)
    
         outputSingleScan(src, [x y]);
         waitforbuttonpress
   
    end 

end
%% 26
o