function [histAll] = TCSPC_iterative( obj, exposure_time_us, VBD_mV, Nframes_it, Nsteps, stepSize, Nspads )
% Close all existing figures
close all

number_of_frames = Nframes_it;
VQUENCH_mV = 1200;
% Set number of bins to be used for histogram generation
number_of_bins          = 600;              % 20 MHz (50 ns) laser, 50 ps, plus a few extra bins just in case...
% IF YOU CHANGE PREVIOUS ALSO CHANGE LINE 317 (EUNAN) 
TDC_resolution          = 50;                       % in [ps]
BIN_width               = 1 * TDC_resolution;       % in [ps] (no biining applied, so 1 is used)

% histogram using all iterations
histAll = zeros(512,number_of_bins);

% % number of iteration
% Niters = 10;
% % number of frames per iteration
% Nframes_it = 10000;

% % number of min DCR SPADs to use per pixel
% Nspads = 8;

% Get bitfile
bitfile = obj.bitfile;
fprintf(1, '\nbitfile = %s\n', bitfile);

% reset chip
fprintf(1, '\nReset chip\n');
obj.ChipReset();

% Set SPAD type
if obj.REDSEL == 1
    spad_type       = 'RED SPAD';
    spad_prefix     = 'red_SPADs';
    spad_VBD        = 'VBDRED';
    BLUESEL         = 0;
    VOP_max_mV      = 28500;
else
    spad_type       = 'BLUE SPAD';
    spad_prefix     = 'blue_SPADs';
    spad_VBD        = 'VBDBLUE';
    BLUESEL         = 1;
    VOP_max_mV      = 17500;
end

% Control flags
% save_results_flag   = 0;
% save_figures_flag   = 1;            % only if plotting 
save_frames_flag    = 0;
plot_flag           = 1;

% if save_results_flag
%     path_prefix     = sprintf('C:/Kufcsak/18_04_13/');
%     % Create results folder if it does not exist already
%     if( ~exist( path_prefix, 'dir' ) )
%         mkdir ( path_prefix );
%     end
% end

% no scanning option for now, update later if needed...

SPADs_enabled               = Nspads;       % Indicates the number of min DCR SPADs enabled per pixel

clk_period                  = 30.0;    % System clk period in ns (= 33.33MHz)

frame_size_updated = number_of_frames;


% Initialise delay codes for SPC_A_GATING
CODE_P1 =   0;
CODE_P2 =   8;

% Initialise delay codes for SPC_B_GATING
CODE_P3 =   0;
CODE_P4 =   0;

% Initialise delay code for PSTOP
CODE_PSTOP = 0;                                     % 0-7 -> By-passing the delay generator

% Initialise delay code for PCAL
CODE_PCAL =  CODE_PSTOP + 100;

obj.SetReg( 'CODE_P1',  CODE_P1 );
obj.SetReg( 'ENB_P1', 1 );                          % 0/1 -> Enable/Disable ( Negative Logic! )

obj.SetReg( 'CODE_P2',  CODE_P2 );
obj.SetReg( 'ENB_P2', 1 );                          % 0/1 -> Enable/Disable ( Negative Logic! )

obj.SetReg( 'CODE_P3',  CODE_P3 );
obj.SetReg( 'ENB_P3', 1 );                          % 0/1 -> Enable/Disable ( Negative Logic! )

obj.SetReg( 'CODE_P4',  CODE_P4 );
obj.SetReg( 'ENB_P4', 1 );                          % 0/1 -> Enable/Disable ( Negative Logic! )

obj.SetReg( 'CODE_PSTOP',  CODE_PSTOP );
obj.SetReg( 'ENB_PSTOP', 0 );                       % 0/1 -> Enable/Disable ( Negative Logic! )

obj.SetReg( 'CODE_PCAL',  CODE_PCAL );
obj.SetReg( 'ENB_PCAL', 0 );                        % 0/1 -> Enable/Disable ( Negative Logic! )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set PIX_BINNING_EN 
% This parameter must be defined even when pixel binning is not used 
% since it is required by the getData function !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PIX_BINNING_EN     = 0;                      
obj.SetReg( 'PIX_BINNING_EN', PIX_BINNING_EN );                      
fprintf(1, '\n');

% Set number of Frames per iteration
obj.SetReg( 'NUMBER_OF_FRAMES', number_of_frames ); 
fprintf(1, '\n');
% Check NUMBER_OF_FRAMES
fprintf(1, '\n Check the current setting for the number of frames:');
f = wireoutdata( obj.okComms, obj.bank, 'NUMBER_OF_FRAMES_REG' );
fprintf(1, '\n * NUMBER_OF_FRAMES_REG: %d\n\n', f);

% Set SPAD voltage
VOP_mV = min( VBD_mV, VOP_max_mV );                 % Make sure specified SPAD voltage is not too high !!!
obj.SetVoltage( spad_VBD, VOP_mV/1000 );

obj.SetVoltage( 'VDDOSC', 3.3 );

% Set DelayGenVDD voltage
obj.SetVoltage( 'DelayGenVDD', 1.2 );

obj.SetVoltage( 'VQUENCH', VQUENCH_mV/1000 );

%TDCREADOUTEN = 1;                                   % Only the first TDC word to be readout
TDCREADOUTEN = 3;                                   % Both TDC words to be readout
obj.SetReg( 'TDCREADOUTEN', TDCREADOUTEN );        
% Set number of bins to be readout for TCSPC mode
if TDCREADOUTEN == 1
    bins_enabled = 1;                               % Reading only first 11-bits TCSPC data
else
    bins_enabled = 2;                               % Reading full 17-bits TCSPC data (including the Left/Right SPAD flag)
end
obj.SetReg( 'bins_enabled', bins_enabled);

if TDCREADOUTEN > 0
    % Disable Histogram mode
    obj.SetReg( 'HIST_EN', 0 );
    
    % Make sure histogram bins are not enabled in raw TCSPC mode
    obj.SetReg( 'HISTREADOUTEN', 0 );               % No histogram bins are enabled
    
    obj.SetReg( 'SPC_LEFT_RIGHT_CNT_EN', 0 );       % 0/1 -> Left/Right SPC Disabled/Enabled (It has no effect in TCSPC mode)
    
    % SPC_A & SPC_B counters
    obj.SetReg( 'SPC_A_EN', 0);                     % It has no effect in TCSPC mode
    obj.SetReg( 'SPC_B_EN', 0);                     % It has no effect in TCSPC mode
    
    obj.SetReg( 'AUTO_SEQ', 0 );                    % 0/1 -> raw TCSPC mode / Histogram mode (Note: Histogram mode does not work with AUTO_SEQ = 0)
    obj.SetReg( 'FAST_MODE', 1 );                   % This must be set to 1 !!! When to use FAST_MODE = 0 ? Investigate this option!
end

obj.SetReg( 'OVERFLOW_PROT_EN', 1 );                % If enabled, all TDC bits will be set to HIGH if an overflow occurs
obj.SetReg( 'RESET_EXTEND', 0 );                    % Extends RESET to the rising edge of the STOP signal

% Settings for TIME/SPAD gating
TIME_GATING_EN  = 0;                                % 0/1 -> Disable/Enable
SPAD_GATING_EN  = 0;                                % 0/1 -> Disable/Enable
obj.SetReg( 'TIME_GATING_EN', TIME_GATING_EN );
obj.SetReg( 'SPAD_GATING_EN', SPAD_GATING_EN );    

obj.SetReg( 'DELGEN_GRO_EN', 1 );                   % 0/1 -> Disable/Enable

OPTCLK_SEL = 1;                                    % 0/1 -> OPTCLK_FROM_FPGA/OPTCLK_TO_FPGA (i.e. 0 -> OPTCLK generated by FPGA, 1 -> OPTCLK generated by laser source)
obj.SetReg( 'OPTCLK_SEL', OPTCLK_SEL );             
% Enable for generating OPTCLK (e.i. LASER_SYNC) by FPGA
if OPTCLK_SEL  == 0
    obj.SetReg( 'clk_10M_enable', 1 );              % 0/1 -> Disable/Enable (required for OPTCLK_FROM_FPGA)
else
    obj.SetReg( 'clk_10M_enable', 0 );              % 0/1 -> Disable/Enable (required for OPTCLK_FROM_FPGA)
end

% This controls enabling the TESTPULSE signal from the firmware side only
% To use TESTPULSE, SRAM also needs to be programmed accordingly (for setting TESTSEL bit) !!!
obj.SetReg( 'TESTPULSE_EN', 1 );                    % 0/1 -> no TESTPULSE generated / a TESTPULSE generated by the firmware

% Re-program SI with the above control register settings
fprintf(1, '\n\nRe-programming SI\n\n');
obj.SIProg();

if not(obj.REDSEL)
    % AK - TMP
    % for blue SPADs
    % row 4
    % SRAM_WORD_A     = hex2dec('10100000');
    % row 5
    % SRAM_WORD_A     = hex2dec('08080000');
     % 4, 5 : 18180000
    % SRAM_WORD_A     = hex2dec('18180000');
    % % 4,5,6 1C1C0000
    % SRAM_WORD_A     = hex2dec('1C1C0000');
    % % 3,4,5,6 3C3C0000
    % SRAM_WORD_A     = hex2dec('3C3C0000');
    % % 3,4,5,6,7: 3E3E0000
    % SRAM_WORD_A     = hex2dec('3E3E0000');
    % % 2,3,4,5,6,7: 7E7E0000
    % SRAM_WORD_A     = hex2dec('7E7E0000');
    % % 2,3,4,5,6,7,8 : 7F7F0000
    % SRAM_WORD_A     = hex2dec('7F7F0000');
    % 1-8: FFFF0000
    SRAM_WORD_A     = hex2dec('FFFF0000');
else
    % for red SPADs
     % row 4
    % SRAM_WORD_A     = hex2dec('00001010');
    % row 5
    % SRAM_WORD_A     = hex2dec('00000808');
    % % 4, 5 : 18180000
    % SRAM_WORD_A     = hex2dec('00001818');
    % % 4,5,6 1C1C0000
    % SRAM_WORD_A     = hex2dec('00001C1C');
    % % 3,4,5,6 3C3C0000
    % SRAM_WORD_A     = hex2dec('00003C3C');
    % % 3,4,5,6,7: 3E3E0000
    % SRAM_WORD_A     = hex2dec('00003E3E');
    % % 2,3,4,5,6,7: 7E7E0000
    % SRAM_WORD_A     = hex2dec('00007E7E');
    % % 2,3,4,5,6,7,8 : 7F7F0000
    % SRAM_WORD_A     = hex2dec('00007F7F');
    % 1-8: FFFF0000
    SRAM_WORD_A     = hex2dec('0000FFFF');
end

% SRAM_WORD_A     = hex2dec('0000FFFF');

obj.SetReg( 'SRAM_WORD_A', SRAM_WORD_A );

% %obj.SetReg( 'SRAM_WORD_B', hex2dec('41') );        % Both TDC and TESTSEL enabled
% %obj.SetReg( 'SRAM_WORD_B', hex2dec('01') );        % TESTSEL enabled
obj.SetReg( 'SRAM_WORD_B', hex2dec('40') );        % TDC enabled
% %obj.SetReg( 'SRAM_WORD_B', hex2dec('7E') );        % TDC enabled and all calibration bits set to 1
% %obj.SetReg( 'SRAM_WORD_B', hex2dec('00') );        % Non enabled (Note: with this setting we get timestamp = 0 for all pixels)

for i=0:511
    obj.SRAMProg();
    %fprintf(1,'Wait for SRAM_WORD_COMPLETE flag for Pixel %d\n', i);
    trig = 0;
    % obj.SetReg( 'SRAM_WORD_B', hex2dec('40')+mod(i,32) );        % TDC enabled
    while trig == 0
        trig = trigoutdata( obj.okComms, obj.bank, 'SRAM_WORD_COMPLETE' );
    end
    %fprintf(1,'SRAM_WORD_COMPLETE flag received for Pixel %d at %s\n', i, char(datetime('now')) );
end
     
% programDCRmap(obj,'X6Y18',VBD_mV,SPADs_enabled)    

% Re-program SI with the above control register settings
fprintf(1, '\n\nRe-programming SI\n\n');
obj.SIProg();

% Initialise arrays
frame_array             = zeros(512,number_of_frames);
    
if plot_flag
    %  Create 2D Mesh
    % Create axes
%     axes1 = axes(figure);
%     
%     time_axis = BIN_width/1000:BIN_width/1000:number_of_bins*BIN_width/1000;
%     
%     xlim(axes1, [min(time_axis), max(time_axis)])
%     ylim(axes1, [1 512])
%     view(axes1, 3)                                          % view(3) sets the default three-dimensional view, az = –37.5, el = 30.
%     
%     grid(axes1, 'on')
%     hold(axes1, 'on')
    plot_hist = imagesc(BIN_width/1000:BIN_width/1000:number_of_bins*BIN_width/1000,1:512,zeros(512,number_of_bins));
%     plot_hist = surf(BIN_width/1000:BIN_width/1000:number_of_bins*BIN_width/1000,1:512,zeros(512,number_of_bins));
    
    % plot_hist = imagesc(zeros(512,number_of_bins));
    
    colorbar
    %caxis([-1 0.9])
    xlabel('Time (ns)')
    ylabel('Pixels')
%     zlabel('Counts')
%     % TODO save histograms?

%         fig2 = figure;

    % fig_FWHM = figure;
    % xlim([1,512]);
end

% Start exposure for number_of_lines 
obj.SetExpMode( exposure_time_us, clk_period );



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                        EUNAN SCAN 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % MULTI POINT SCAN? 1 = YES      
 
multipoint = 0;
infile = '91pointsNEW.txt';
points = load(infile, '-ascii');

if multipoint == 1 
    
    % FIBER GALVO SYSTEM
    
daq.reset
daq.HardwareInfo.getInstance('DisableReferenceClockSynchronization',true); % This is to overcome a different issue
d = daq.getDevices;
src = daq.createSession('ni');
ch1=addAnalogOutputChannel(src,'Dev1','ao0','Voltage');
ch2=addAnalogOutputChannel(src,'Dev1','ao1','Voltage');
src.Rate = 1000; %speed it execiutes the que

points1 = points; 

elseif multipoint == 0

    points = 1; % INDUVIDUAL POINT SCAN
end

% For Each Point we want a scan
for p = 1:size(points,1)

   
outputSingleScan(src, [points1(p,1) points1(p,2)]);
 fprintf('%d ', p);
 disp(p)
  
                    
pos = 0:stepSize:Nsteps*stepSize; % Actual Motor Positions
NumPIX = 512;
NumBINS = 600;
% 
 allHists = zeros(Nsteps, NumPIX, NumBINS); % initiliase large 3D array, Slices x 512 x bins
 m1 = motor;                                % Talk to the motor
 a = motor.listdevices;
 
 connect(m1,a{1});


% Take Scan Based on Input Parameters

for n_pos = 1 % How many positions
   

      currentPos = pos(n_pos); %     the actual position at this step
      fprintf('%d ', n_pos);
      moveto(m1,currentPos);    %     Go there
    fprintf(1, '%s', n_pos);
     
     
     
%     this starts the acquisiton of a TCSPC histogram (many many lines of timestamps)
    trigger(obj.okComms, obj.bank, 'EXPOSURE_START');
    
  


    % For each frame, bins_enabled*8 bytes are readout in scanning mode and bins_enabled*2*512 bytes are readout in non-scanning mode (see getData function)
%     this is the raw data of the many may lines of timestamps that will be
% used to calculate the timestamp, and the that built to a histogram
    frames = obj.getData( bins_enabled, frame_size_updated ); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decode sensor data (i.e. raw data to timestamp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    data_words = reshape(frames, 512*bins_enabled, number_of_frames);               % Reshape into one frame per column format
    
    for frame_no=1:number_of_frames
    
        % Check if only the first TDC word is requested for readout ( i.e. TDCOUT<10:0> )
        if bins_enabled == 1
            
            
            % Invert FINE bits (i.e. bin0<2:0>) to fix a design bug
            bin0            = bitxor(bin0, 7);
            
            % Making sure higher 5 MSBs are all set to zero (so using only first 11-bits)
            % This might be redundant !!!
            frame_data      = double (bitand( 2047, bin0 ) );
            
        end    

        % Check if both TDC words were requested for readout (i.e. First word = TDCOUT<10:0>, Second word = { LR_OUT, TDCOUT<15:11> } )
        if bins_enabled == 2
            
            %
            
           
            
            % Get frame data
            frame64_16  = reshape(data_words(:,frame_no), 64, 16);              % Each column holds a block of 64 words (8 columns for BIN0 + 8 columns for BIN1)
            frame128_8  = reshape(frame64_16, 128, 8);                          % Each column holds a block of 128 words (64 rows for BIN0 + 64 rows for BIN1)
            
            bin0_64_8   = frame128_8(1:64, 1:8);                                % Get top 64 rows for BIN0
            bin0        = reshape(bin0_64_8, [], 1);                            % Reshape bin0_64_8 into a single column format
            
            bin1_64_8   = frame128_8(65:128, 1:8);                              % Get bottom 64 rows for BIN1
            bin1        = reshape(bin1_64_8, [], 1);                            % Reshape bin1_64_8 into a single column format
            
            % Invert FINE bits (i.e. bin0<2:0>) to fix a design bug
            bin0 = bitxor(bin0, 7);
            
            % Separate Left/Right SPAD flag which is MSB of the 2nd TDC word (i.e. bin1<5>)
            %LR_SPAD = bitget(bin1, 6);
            
            % We need to set higher MSBs to zero so that read data bits can be selected with dec2bin function
            bin0 = bitand( 2047, bin0 );            % Making sure that bit 11 and all higher bits are set to zero (i.e., using only 11 LSBs)
            bin1 = bitand(   31, bin1 );            % Making sure that bit 5 and all higher bits are set to zero (i.e., using only 5 LSBs)
            
          

            % extend vars to 32 bit unsigned integers
            bin0_32 = uint32(bin0);
            bin1_32 = uint32(bin1);

            % shift and concatenate
            frame_data = bitshift(bin1_32,11) + bin0_32;

        end % for < if ( bins_enabled == 2 ) >

        % Re-order data words
        frame_data64_8              = reshape( frame_data, 64, 8 );
        frame_data_ordered          = reshape( frame_data64_8', 512, 1 );
        
        % Save current frame
        frame_array( :, frame_no ) = frame_data_ordered;

    end % for loop: frame_no
    
%     toc_decoding_frames = toc(tStart_decoding_frames);
%     
%     fprintf(1, '\nElapsed time for decoding %d frames: %s seconds\n\n', number_of_frames, toc_decoding_frames );
    
% this is sving the line of timestamps, probably you  don't need?
    if save_frames_flag
        save(sprintf('frame_array_%d',n_iter),frame_array);
        % i = i+1;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process data - this is where the timestamps are built into a histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % histogramming
    counts = zeros(512,number_of_bins);
    clear N;
    
    
% this is the actual loop to build the histogram
    for pix=1:512
        [N, ~] = histcounts( frame_array(pix,:), 'BinMethod', 'integers' );    
        % Exclude BIN0 (i.e. zero counts which might be too high!)
%         counts(pix, 1:numel(N)) = N(1:end);       
        counts(pix, 1:numel(N)-1) = N(2:end);
%           counts(pix, :) = N(1:number_of_bins);
    end
    
    % plot histograms - if you need to, set plot_flag = 1
    if plot_flag
        
%         if max(max(counts))   
%           zlim(axes1, [0 max(max(counts))])
%         end
        
        plot_hist.CData = counts;
% S        plot_hist.ZData = counts;
        drawnow limitrate;
%       

    end

   

%     currentHist is the histogram per position
    currentHist = counts(:,1:number_of_bins);
    
      % you may want to save to a large 3D array, eg.
   allHists(n_pos,:,:) = currentHist;
   
 
   

end
    % Save point scan & repeat for 
    save(['C:\temp\TCSPC\ScanResults' num2str(p) '.mat'],'allHists');
end 