%% Vocoding workshop:

clear all
clc
close all
%% 0. Load speech to be vocoded:

[x, fs] = audioread('Female.wav');

% stim = audioplayer(x, fs, 16);
% playblocking(stim);

% 
% plot((0:length(x)-1)/fs, x);
% xlabel('Time (s)')
% ylabel('Amplitude')

y = zeros(size(x)); %y is going to be the filtered signal

%% 1. Create the filterbank (Analysis filters):

%Define audio sampling frequency:

%fs = 44100; %Hz

%Define number of bands (channels; electrodes)

nchs = 8; %Try 4, 16, 22 => these correspond to the number of electrodes in the system

%Create a list of log-spaced frequencies defining the cutoff frequencies of 8 frequency bands between 150 and 7500 Hz.

cutoffs = logspace(log10(150), log10(7500), nchs+1); %logspace(log10(150), log10(7500), nchs+1); %linspace(150, 7500, nchs+1); %Here you can change the spacing: log/lin/ something else, 
%and you can also try changing the frequency range. You have nchs, thus
%it gives nchs+1 cutoff frequencies

%Define filter order. This indicates the sharpness of your filter. The
%larger the order, the sharper your filter. You can play around with that
%to simulate channel interaction between neighbouring electrodes. Start
%with 8.
filt.ord = 8;

%Define the filterbank
filterbank = struct();
for i=1:nchs
      filterbank(i).cutoffs = [cutoffs(i) cutoffs(i+1)];
      [b, a] = butter(filt.ord/4, filterbank(i).cutoffs*2./fs); %Butterworth bandpass filter with flat top. The function 'butter' doubles the order, thus you need to feed it filt.ord/2.
      %Also, because we are going to use filtfilt which doubles the order
      %again, so we need to divide by another 2.
      filterbank(i).b = b;
      filterbank(i).a = a;
      
      %Plot freq response of each filterbank to see what your filters look
      %like
      [h,f] = freqz(b,a,1024,fs); %The first filter looks funny because there is something weird going on in the freqz func
      %bec filterbank(1).a is very small so you reach the limits of precision
      %of the machine.
      figure(1)
      plot(f,20*log10(abs(h))); ylim([-30 0])
      hold on
      
      %Run the signal through the filterbank
      filterbank(i).x = filtfilt(b, a, x);% Run the signal through the filterbank
      
      %You can try listening to the sum of the output of all filters to
      %listen for which components have been dropped out. Try changing the
      %cutoffs to only low (e.g. 200-2000 Hz) or only high frequencies (e.g. 3000-7500) and see what happens.
      y = y + filterbank(i).x;
      
      
    
      filterbank(i).y = filterbank(i).x;
      
%        stim = audioplayer(filterbank(i).y, fs, 16);
%        playblocking(stim)
      
      %Plot the output signal of each filter
      figure(2)
      plot(i+filterbank(i).y);
      hold on
      
      %Listen to the output of each filter
      %soundsc(filterbank(i).y, fs) %REPLACE WITH AN AUDIOPLAYER obj
end

% stim = audioplayer(y, fs, 16);
% playblocking(stim)

% stim = audioplayer(x, fs, 16);
% playblocking(stim)
        
        


%You can try listening to the sum of the output of all filters to
%listen for which components have been dropped out. Try changing the
%cutoffs to only low (e.g. 200-2000 Hz) or only high frequencies (e.g. 3000-7500) and see what happens.
%soundsc(y, fs)

%% 2. Envelope extraction by half-wave rectification and low-pass filtering


%Define the low-pass filter (LPF) to extract the envelope:
[b_low, a_low] = butter(1, 250*2/fs, 'low'); %This is a fourth order LPF, with a cutoff of 250 Hz

for i=1:nchs
    
    %1. Half-wave rectify the signal
    filterbank(i).hwr = filterbank(i).x;
    filterbank(i).hwr(filterbank(i).hwr < 0) = 0; %Use max(x,0) instead of lines 86-87
    
    %2. Low-pass filter the half-wave-rectified signal to obtain the
    %envelope
    filterbank(i).env = filtfilt(b_low, a_low, filterbank(i).hwr);
    
    %Plot the envelope...
    figure(3)
    plot(filterbank(i).x + .5*i, '-k');
    hold on
    plot(filterbank(i).env*2 + .5*i, '-r');
    
    
end


%% 3. Use envelope to modulate electric pulse train

% y = zeros(length(x),1);
% for i = 1:nbands
%     
%     %1. Create the pulse train:
%     pulse = zeros(length(filterbank(i).env),1);
%     f0 = 500; %Typical channel stimulation rate in Cochlear devices
%     nT = round(fs/f0);
%     pulse(1:nT:end) = 1;
%     pulse(2:nT:end) = -1;
%     pulse = cosgate(pulse,fs,5e-3); %cosine ramp the pulse train to avoid ringing in your filtered signal
%     
% %     figure()
% %     plot(pulse)
% 
% %create square wave:
% %bla = square(2*(500)/pi*(1:length(x)));
% 
%     
%     
%     filterbank(i).mod = filterbank(i).env.*pulse;
%     
%     filterbank(i).mod = filtfilt(filterbank(i).b, filterbank(i).a, filterbank(i).mod);
%     
%     %Plot the modulated signal
%     figure(4)
%     plot(i+filterbank(i).mod)
%     title('Instantaneous pulse with frequency not matched')
%     hold on
%     
%     %Sum output of all channels and listen to it. This is the vocoded speech
%     y = y + filterbank(i).mod;
%     
% end
% 
% y = y*rms(x)./rms(y);
% soundsc(y,fs)
% figure(); plot(y)

%That was quite robotic, wasn't it? and not very clear either... In
%reality, pulses are not quite instantaneous like we modelled them. They
%have a period, so they look more square...
%And in simulation, turns out that you need to create square/sine waves whose frequencies
%correspond to the center freq of the channel...

y = zeros(length(x),1);
t = (0:length(x)-1)'/fs;

for i = 1:nchs
    
    
    %1. Create the pulse train:
    filterbank(i).center_f = sqrt(filterbank(i).cutoffs(1) .* filterbank(i).cutoffs(2));
    
    %Try square wave
    pulse = sin(2*pi*(filterbank(i).center_f)*t); 
    pulse(pulse > 0) = max(pulse);
    pulse(pulse < 0) = min(pulse);
    %Square wave is too abrupt; needs to be smoother because sharpness of
    %the pulse leads to ringing effects in modulation => plot 'y' and
    %see...
    
    %Try sine wave because transitions are smoother... => In NH we need a
    %carrier signal to carry the signal to the right place...
    %pulse = sin(2*pi*(filterbank(i).center_f)*t); 

    
    
    filterbank(i).mod = filterbank(i).env.*pulse;
    
    filterbank(i).mod = filtfilt(filterbank(i).b, filterbank(i).a, filterbank(i).mod);
    
    %Plot the modulated signal
    figure(5)
    plot(i+filterbank(i).mod)
    title('Pulse with frequency matched to center frequency of the band')
    hold on
    
    %Sum output of all channels and listen to it. This is the vocoded speech
    y = y + filterbank(i).mod;
    
end

y = y*rms(x)./rms(y);
 soundsc(y,fs)
% figure()
% plot(y)

%Not v clear... Try increasing number of channels to 22 instead of 8... Notice that as you
%increase the number of channels, the signal becomes clearer. Increase even
%more to say 40. What do you notice? The pitch fluctuations of the speaker
%become clearer...

%You could also try increasing the order of the filters to 12 instead of 8.
%Set nbands back to 8 and
%try again...

%Because cochlear implant users do not perceive harmonic structures, a
%noise carrier is often preferred to a sine carrier in vocoder
%simulations...

y = zeros(length(x),1);
t = (0:length(x)-1)'/fs;

for i = 1:nchs
    
    
    %1. Create the broadband noise signal:
    %noise carrier:
    nz = rand(size(x))*2-1; %To make it go from -1 to 1 to make it look more like a waveform. This does not create peaks, thus no spectral properties
    %(harmonic structure) when u just have a broadband noise.
    
    filterbank(i).mod = filterbank(i).env.*nz;
    
    filterbank(i).mod = filtfilt(filterbank(i).b, filterbank(i).a, filterbank(i).mod);
    
    %Broadband noise sounds terrible, let's look at narrowband noise...
    %nz_narrow = filtfilt(filterbank(i).b, filterbank(i).a, nz);
    %filterbank(i).mod = filterbank(i).env.*nz_narrow;
    
    %This is much better right?

    
    
    
    
    %Plot the modulated signal
    figure(6)
    plot(i+filterbank(i).mod)
    title('Modulated with noise carrier')
    hold on
    
    %Sum output of all channels and listen to it. This is the vocoded speech
    y = y + filterbank(i).mod;
    
end

y = y*rms(x)./rms(y);
soundsc(y,fs)
 
  



%% 4. Changing Synthesis filters
%%In a CI, you need to carry the temp env to the specific place on the
%cochlea, thus you need a carrier, which is the pulse train.
%Location of most apical electrode at the best case is around 500 Hz; MEDEL
%has a longer array that reaches about 300 Hz; but you still need to
%preserve the information transmitted below 1 kHz, so you need to preserve
%them and transpose them to higher freqs (mapping 200 Hz from the freq map
%to 500-700 Hz along the cochlea) => The frequency mismatch case.

%Let's model that...

%Frequencies along the basilar membrane are defined by the Greenwood
%formula...

 bm_mm = linspace(frq2mm(500), frq2mm(16000), nchs+1); % Create linearly spaced locations in mm along basilar membrane from location of most apical to most basal electrode
 bm_freq = mm2frq(bm_mm); %Convert locations in mm back to frequency according to Greenwood formula.


y = zeros(length(x),1);
t = (0:length(x)-1)'/fs;

n_ord_synth = 8; %Define order of synthesis filters. It can be different from that of analysis filters

for i = 1:nchs
    
    
    %1. Create the broadband noise signal:
    %noise carrier:
    nz = rand(size(x))*2-1; %To make it go from -1 to 1 to make it look more like a waveform. This does not create peaks, thus no spectral properties
    %(harmonic structure) when u just have a broadband noise.
    
    
    
    
    %In typical CI simulations, you modulate first with broadband noise,
    %AND THEN you appy the synthesis filters:
    filterbank(i).mod = filterbank(i).env.*nz;
    
    filterbank(i).synth = [bm_freq(i) bm_freq(i+1)];
    [filterbank(i).synth_b, filterbank(i).synth_a] = butter(n_ord_synth/4, filterbank(i).synth.*2./fs);
    filterbank(i).mod = filtfilt(filterbank(i).synth_b, filterbank(i).synth_a, filterbank(i).mod);

    
    
    
    
    %Plot the modulated signal
    figure(7)
    plot(i+filterbank(i).mod)
    title('Simulating mismatch')
    hold on
    
    %Sum output of all channels and listen to it. This is the vocoded speech
    y = y + filterbank(i).mod;
    
end

y = y*rms(x)./rms(y);
soundsc(y,fs)

%Plot the spectrograms of vocoded and original signal and notice the
%terrible frequency and time (spectrotemporal) resolutions of the vocoded
%one relative to the original one.
figure();
n = 2^11;%1024; %window size; fft is optimized to powers fo 2. The smaller n is, the smoother time is and the coarser your frequency and vice versa, so optimize.
w = hann(n);
[S,F,T] = spectrogram(y,w,round(n*.9),n,fs,'yaxis');

logS = 20*log10(abs(S)); %transform ampl to dB
m = max(logS(:)); %To find the peak value in the matrix
surf(T, F, logS);
shading flat
xlim([0 T(end)])
ylim([0 8000])
view(0,90)
set(gca, 'CLim', [m-70, m])
colorbar()

title('Vocoded signal')
xlabel('Time (s)', 'FontSize', 25)
ylabel('Frequency (Hz)', 'FontSize', 25)
zlabel('Amplitude (dB)', 'FontSize', 25)


figure();
n = 2^11;%1024; %window size; fft is optimized to powers fo 2. The smaller n is, the smoother time is and the coarser your frequency and vice versa, so optimize.
w = hann(n);
[S,F,T] = spectrogram(x,w,round(n*.9),n,fs,'yaxis');

logS = 20*log10(abs(S)); %transform ampl to dB
m = max(logS(:)); %To find the peak value in the matrix
surf(T, F, logS);
shading flat
xlim([0 T(end)])
ylim([0 8000])
view(0,90)
set(gca, 'CLim', [m-70, m])
colorbar()

title('Original signal')
xlabel('Time (s)', 'FontSize', 25)
ylabel('Frequency (Hz)', 'FontSize', 25)
zlabel('Amplitude (dB)', 'FontSize', 25)

