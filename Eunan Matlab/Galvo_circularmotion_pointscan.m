daq.reset
daq.HardwareInfo.getInstance('DisableReferenceClockSynchronization',true);  % This is to overcome a different issue
d = daq.getDevices;
src = daq.createSession('ni');
ch1=addAnalogOutputChannel(src,'Dev2','ao0','Voltage');
ch2=addAnalogOutputChannel(src,'Dev2','ao1','Voltage');
src.Rate = 1000; %speed it execiutes the que
resolution=100; 
r=0.43;
theta = 0:2*pi/resolution:2*pi;

for i=1:resolution
    
    points(i,1)=r*cos(theta(i));
    points(i,2)=r*sin(theta(i));
end
 
for k=1:1
    
for i=1:length(points)
    disp(i) 
    outputSingleScan(src,points(i,:));
    %pause(0.1)
    disp(points(i,:));   
    disp([r theta(i)]); 
    pause
end

end

pause(2)
queueOutputData(src,[outputSignal1 outputSignal2]);
src.startBackground;

%scanning all points SAMPLE 1
% outputSingleScan(src, [0.24 -0.12]);
% pause(1)
% outputSingleScan(src, [0.08 -0.23]);
% pause(1)
% outputSingleScan(src, [-0.14 -0.25]);
% pause(1)
% outputSingleScan(src, [0.26 0.16]);
% pause(1)
% outputSingleScan(src, [0.22 0.3]);
% pause(1)
% outputSingleScan(src, [0.07 0.31]);
% pause(1)
% outputSingleScan(src, [0.26 0.26]);
% pause(1)
% outputSingleScan(src, [-0.39 0.2]);
% pause(1)
% outputSingleScan(src, [-0.14 0.33]);
% pause(1)
% outputSingleScan(src, [-0.26 0.31]);
% pause(1)
% outputSingleScan(src, [-0.43 0.14]);
% pause(1)
% outputSingleScan(src, [-0.41 0]);