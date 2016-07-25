function [Speech_on] = finalbatch_delay(wav_name)

%-------------------------------------------------------------------------------
% File Name   : epdbatch.m
% Description : This program detects the onset and offset of the speech for the 
%               input wave files.
%-------------------------------------------------------------------------------              
% Author      : Prateek Bansal
%-------------------------------------------------------------------------------
% Revised by Zenzi Griffin June 3, 2003 to skip initial times (250 ms) too early to
% be real response.
%-------------------------------------------------------------------------------              

% clear all;
% close all;

window_size = 10;    % Its in Milliseconds
end_silence = 100;   % Its in Milliseconds
back_off    = 1;   % In milliseconds 
%modified May 24

% Create an output file
%Finput = input(' Name of the output data file (.txt):','s');
%Fop = fopen(Finput,'a');

% Directory Structure

wavfiles =wav_name;
nosfiles = length(wavfiles);
if nosfiles == 0
   fprintf(1,' No Wave Files found \n');
   return;
end   
filescell = sort({wavfiles(1:nosfiles).name});

kcount = 1;
for icount = 1:nosfiles
      FILE = char(filescell(icount));
      ft1 = char(FILE);
            
      % Read the Wave file
      [input_sig Fs] = wavread2(FILE);
      if size(input_sig,2) ==2
        input_sig = input_sig(:,2);
      end;
      % Normalized the input signal
      input_sig = input_sig/max(input_sig);
 
      % High Pass filtered at 100Hz to remove dc content and AC hum
      W1 = 100 *(2/Fs) ; 

      % Low Pass Filtered at 4KHz
      W2 = (4*10^3)*(2/Fs);

      % Design the bandpass filter
      [B,A]= butter(10,[W1,W2]);

      % Filtered Signal Output
      filt_inpsig = input_sig;

      % Clear the variables
      clear W1 W2 input_sig;
      
      % Normalize the filtered signal
      filt_inpsig = filt_inpsig/max(filt_inpsig);
      
      % Length of filtered data
      signal_length = length(filt_inpsig);

      % Define the time axis
      time_t = 0 : 1 : (signal_length-1); 

      % WindowSize in Samples
      win_samples = round((Fs*window_size)/1000);
      
      % Define the Hamming Window
      ham_window = hamming(win_samples);

      j = 1;
      for i = 1: win_samples : signal_length - win_samples
           temp = filt_inpsig(i:i+win_samples-1);
           temp2 = temp.*ham_window;
           temp2 = abs(filter(B,A,temp2));
           energy(j) = sum(temp2);
           j = j+1;
      end
      
      % Compute of End of Silence in Samples
      end_silence_sample = round((end_silence*Fs)/1000);

      % Define the Silence Range
      silenceRange = 1:length(1:win_samples:end_silence_sample - win_samples);
      
      % Clear the variables
      clear ham_window win_samples signal_length temp temp2

      % Energy Thresholds
      % IMN (silence energy) is average energy for initial 100ms of signal.

%      IMX=max(energy);                       % peak Energy
%      IMN=mean(energy(silenceRange));        % Silence energy
%      I1=0.03*(IMX-IMN)+IMN;                 

%      I2=4*IMN;
%      ITL=min(I1,I2);                        % Lower Threshold
       ITL = mean(energy) * 0.2;
       ITU=5*ITL;                             % Upper Threshold

      % Determination of End Points

      N1=0;        % Start point initial estimate
      N2=0;        % End point initial estimate

      duration=length(energy);

      done=0;

      % Estimation of the start point based on energy considerations
      % to start estimate at the very start of the soundfile rather than 
      % 250 ms later, set change 25 to 1 in "for m=25:duration"

      for m=back_off:duration
         if and(energy(m)>=ITL,~done)
            for i=m+2:duration
               if energy(i)<ITL
                  break
               else   
                  if energy(i)>=ITU
                     if ~done
                        N1=i-(i==m);
                        done=1;
                     end
                     break
                  end
               end
            end
         end   
      end
      
      done=0;

      % Estimation of the end point based on energy considerations

      for m=duration:-1:1
         if and(energy(m)>=ITL,~done)
            for i=m-2:-1:1
               if energy(i)<ITL
                  break;
               else
                  if energy(i)>=ITU
                     if ~done
                        N2=i+(i==m);
                        done=1;
                     end
                     break
                   end
                end
             end   
          end
       end

       warpRatio=round(length(filt_inpsig)/length(energy));
       N1_w=N1*warpRatio;
       N2_w=N2*warpRatio;

       % Plotting the signal
       % The position of the onset & offset lines in the plot don't tend to
       % look accurate, so probably best to keep this option off until figured 
		% out.
       p = 0;

       if p
          subplot(3,1,1);
          plot(time_t,input_sig);
          hold on;
          line([N1_w,N1_w],[-1,1],'Color',[1 0 1]);
          line([N2_w,N2_w],[-1,1],'Color',[1 0 1]);
          subplot(3,1,2);
          plot(1:length(energy),energy);
          subplot(3,1,3);
          plot(1:length(zero_cross),zero_cross);
       end
       
       % Evaluate Sppech_on and Speech_off
       Speech_on = round((N1_w*1000)/Fs);
       Speech_off = round((N2_w*1000)/Fs);
       Speech_on  = round(Speech_on - 0.05*Speech_on);
       
       %fprintf(Fop,'   %s       %s         %d         %d    \n',ft1(1:2),ft1(3:length(ft1)-4),Speech_on,Speech_off);
       %clear Ft1 Ft2 Ft3 Ft4 Speech_on Speech_off N1 N2 N1_w N2_w  warpRatio energy ITU ITL
end   

%fclose(Fop);
