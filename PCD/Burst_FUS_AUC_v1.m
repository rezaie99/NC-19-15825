%% Main Function for analyzing Burst FUS PCD data 
% Author Mehmet Ozdas, zdasm@ethz.ch, ETH Zurich, Neurotechnology
% Laboratory, 2020
% 2.5MHz H-147 Transducer, Sonic Concepts
% Power amplifier: E&I 325LA for Transducer.50dB 
% Function generator: Agilent 33210A to triger Transducer
% Y107   PCD, Sonic Concepts
% 20dB, Broadband RF Amplifier 1MHz-1GHz, Ramsey Electronics
% Sampling Freq=125MS/s, 15 bits resolution 
% Picoscope: 5042 for data acquisation
% PicoScope: 3205B used for TTL to drive PRF of Picoscope: 5042

clear all;
%% Put directory which contains all mat.file here
%cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_5_Burst_25mV/UC_Carriers')
cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_6_Burst_FUS_100mV/UC_Carriers')
%cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_4_Burst_FUS_200mV/UC_Carriers')
%% isert variables
Fs=125000000; % 125MS/s
recording_window=20; % 20 ms
pulse_duration=2 % ms
t = 0:1/Fs:1-1/Fs;
window_length=Fs/(1000/recording_window)
pulse_size=Fs*(pulse_duration/1000) % number of samples for recording_window
%fprintf('Windows Length in ms :%i\n',pulse_size )
delay_after_TTL=77*125 % 77 us
psdx_size=pulse_size/2

%% AUC Inputs
%center_f=2.5*1000000
%fn=center_f*3 % 7.5 MHz
%fn_window=center_f/50 % in kHz
%% Load .mat files in an array and take FFT
matfiles = dir('*.mat');
N = length(matfiles);
filenames = cell(N,1);
data_array=zeros(N,window_length);
psdx_array=zeros(N,psdx_size);

%% AUC 
%freq array has 125000 freqs for 62.5MHz, half of main Sampling rate.
%freq(5000)=2.5 MHz
data_auc_array_uh= zeros(N,1);  % ultraharmonics 3.75, 6.25 MHz
data_auc_array_ih= zeros(N,3);  % integer harmonics, 7.5, 10, 12.5 MHz
data_auc_array_bbe= zeros(N,3); % Broadband emissions, Area between 2.5-5, 5-7.5, 7.5-10

freqs_uh=[6.25].*2000  % multuplied with 2000 to reach corresponding index in freq
window_uh=30*2 % 50 kHz before and after uh based on the example data
freqs_ih=[7.5, 10, 12.5].*2000
window_ih=40*2 % 40kHz before and after ih based on the example data
freqs_bbe=[2.5 5 7.5].*2000
window_bbe=250*2 % cut 250kHz around ih.

for i = 1:N
   %thisfig = figure();
   filenames{i} = matfiles(i).name;
   transit = load(filenames{i});
   data_array(i,:)=transpose(transit.A);
   data_tr=data_array(i,:);
   size_A=size(transit.A);
   start_pulse=round(2*size_A(1)/20)+delay_after_TTL; % pulse starts from here, based on picoscope TTL
   x=data_tr([start_pulse:start_pulse+pulse_size-1]);
   x=x.*hamming(length(x))'; % applying hamming window
   n = length(x);
   xdft = fft(x);
   xdft = xdft(1:n/2);
   psdx = (1/(Fs*n)) * abs(xdft).^2;
   freq = 1:Fs/length(x):Fs/2;
   freq=freq/1000000;
   data_auc=psdx; 
   
   count=1
   for j=freqs_uh  % calculate AUC for ultraharmonics
        data_auc_array_uh(i,count)=trapz(data_auc([j-window_uh+1:j+window_uh+1]));
        disp(count)
        count=count+1
   end 
   
   count=1
   for j=freqs_ih  % calculate AUC for integer harmonics
        data_auc_array_ih(i,count)=trapz(data_auc([j-window_ih+1:j+window_ih+1]));
        disp(count)
        count=count+1
   end 

      count=1
   for j=freqs_bbe  % calculate AUC for broadband emissions, exclude 250 kHz before and after integer harmonics
        data_auc_array_bbe(i,count)=trapz(data_auc([j+window_bbe+1:j+2.5*2000-window_bbe+1]));
        disp(count)
        count=count+1
   end 
   
   psdx_array(i,:)=psdx;  % store in array
end

%% Bubbles Mean
% 
% figure();
% x0=10;
% y0=10;
% width=1200;
% height=600
% set(gcf,'position',[x0,y0,width,height]);
% set(gca,'fontsize',20);
%%

ucc_uh_mean=sum(mean(data_auc_array_uh))
ucc_ih_mean=sum(mean(data_auc_array_ih))
ucc_bbe_mean=sum(mean(data_auc_array_bbe))

mean_psdx=mean(psdx_array);

%plot(freq([start_freq+1:stop_freq+1]),10*log10(mean_psdx([start_freq+1:stop_freq+1])),'red');
%plot(freq,10*log10(mean_psdx),'red');
%xlim([6.2 6.3])
%ylim([-170 -60])
%xlabel('Frequency MHz')
%ylabel('Power')
%set(gca,'fontsize',16)
%hold on

%%
%cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_5_Burst_25mV/Baseline')
cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_6_Burst_FUS_100mV/Baseline')
%cd('/Users/mehmetozdas/Desktop/Scripts/pcd_scripts/PCD_data/2020_04_09_NBBB179_H147_9/Site_4_Burst_FUS_200mV/Baseline')
matfiles = dir('*.mat');
N = length(matfiles);
filenames = cell(N,1);
data_array=zeros(N,window_length);
psdx_array=zeros(N,psdx_size);

%% AUC 
%freq array has 125000 freqs for 62.5MHz, half of main Sampling rate.
%freq(5000)=2.5 MHz
data_auc_array_uh= zeros(N,2);  % ultraharmonics 3.75, 6.25 MHz
data_auc_array_ih= zeros(N,3);  % integer harmonics, 7.5, 10, 12.5 MHz
data_auc_array_bbe= zeros(N,3); % Broadband emissions, Area between 2.5-5, 5-7.5, 7.5-10
for i = 1:N
   %thisfig = figure();
   filenames{i} = matfiles(i).name;
   transit = load(filenames{i});
   data_array(i,:)=transpose(transit.A);
   data_tr=data_array(i,:);
   size_A=size(transit.A);
   start_pulse=round(2*size_A(1)/20)+delay_after_TTL; % pulse starts from here, based on picoscope TTL
   x=data_tr([start_pulse:start_pulse+pulse_size-1]);
   x=x.*hamming(length(x))'; % applying hamming window
   n = length(x);
   xdft = fft(x);
   xdft = xdft(1:n/2);
   psdx = (1/(Fs*n)) * abs(xdft).^2;
   freq = 1:Fs/length(x):Fs/2;
   freq=freq/1000000;
   data_auc=psdx; 
   
   count=1
   for j=freqs_uh  % calculate AUC for ultraharmonics
        data_auc_array_uh(i,count)=trapz(data_auc([j-window_uh+1:j+window_uh+1]));
        disp(count)
        count=count+1
   end 
   
   count=1
   for j=freqs_ih  % calculate AUC for integer harmonics
        data_auc_array_ih(i,count)=trapz(data_auc([j-window_ih+1:j+window_ih+1]));
        disp(count)
        count=count+1
   end 

      count=1
   for j=freqs_bbe  % calculate AUC for broadband emissions, exclude 250 kHz before and after integer harmonics
        data_auc_array_bbe(i,count)=trapz(data_auc([j+window_bbe+1:j+2.5*2000-window_bbe+1]));
        dsisp(count)
        count=count+1
   end 
   
   psdx_array(i,:)=psdx;  % store in array
end

%% Bubbles Mean
% 
% figure();
% x0=10;
% y0=10;
% width=1200;
% height=600
% set(gcf,'position',[x0,y0,width,height]);
% set(gca,'fontsize',20);
%%
baseline_uh_mean=sum(mean(data_auc_array_uh))
baseline_ih_mean=sum(mean(data_auc_array_ih))
baseline_bbe_mean=sum(mean(data_auc_array_bbe))

ucc_uh_mean/baseline_uh_mean
ucc_ih_mean/baseline_ih_mean
ucc_bbe_mean/baseline_bbe_mean