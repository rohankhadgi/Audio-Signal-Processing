% ECE 3304 - Discrete Time Signal Analysis
% Dr. Tanja Karp
% Spring 2019 Take Home Examination Portion
% Submitted by: Rohan Khadgi
% R#: R11482763
% Date: 05/10/2019

clear;
clc;

[x, Fs] = audioread('noisy_signal.wav'); % reading the audio signal
T = 1/Fs; % time period (reciprocal of sampling frequency)
sound(x,Fs);
N = length(x);

%% Question 1:

% The 4-second audio signal has a sinusoid which significantly increases
% from low to high frequency throughout the duration of time and low 
% frequency white noise which remains constant.

%% Question 2:

% We know, Normalized Frequency(Omega) is given by the signal k
% multiplied by twice the physical frequency (Sampling Frequency Fs in Hz)
% and divided by the length of the given signal. Mathematically,
% Normalized Frequency:
% (Omega) = (Range of frequency X 2*Sampling Frequency)/Length of signal

freq=(1:N).*Fs./N;
x_fft = fft(x); % fast fourier transform of the signal

plot(x_fft)

% Graphing the FFT of x over normalized frequency
figure(1)
subplot(121)
plot(2*freq./Fs,abs(x_fft))
grid on
title('FFT of X')
xlabel('Normalized Frequency (\times\pi rad/sample)') 
ylabel('X[k]')  
xlim([0 1]) 

% Source Used:
% http://www.informit.com/articles/article.aspx?p=1967020&seqNum=4

%% Question 3:

% The range of frequencies that are part of the signal range from 
% 0 to 6000Hz and the frequency band in which the noise is present is 
% between 900Hz to 2100Hz. A DTF was performed to convert given signal in
% time domain to frequency domain. The FFT of the signal was graphed 
% against Normalized Frequency and it was observed in the graph that 
% the disturbance in noise was seen in between 900 to 2100Hz which 
% corresponded to the noise in the audio signal.

x_psd = abs(x_fft).^2; % power spectral density of the fft(x)

% Graphing Power Spectral Density over Normalized Frequency
figure(1)
subplot(122)
stem(2*freq./Fs,x_psd)
grid on
title('Power Spectral Density of X')
xlabel('Normalized Frequency (\times\pi rad/sample)') 
ylabel('Power Spectral Density abs(X[k])^2')
xlim([0 1])

% Source used:
% Dr. Tanja Karp - fft_example.m, fft_finite_x.m
% Video: https://www.youtube.com/watch?v=VFt3UVw7VrE&t=738s
%% Question 4:

% Using Bandstop filter with Chebyshev Type 1 design method

% Since the range of frequency in which the noise existed was in between 
% the range of frequency in which the sinusoidal signal was present 
% (overlapped where both frequency matched), it was necessary to attenuate
% only the range of frequency in which the noise was present. For this
% reason, I used a Bandstop IIR filter with Chebyshev Type 1 design method.
% The filter order was specified to be 6 as given in the question and
% frequency specifications were set as follows:
% Fs = 16000Hz (obtained when audioread was performed)
% Fpass1=850Hz (Frequency below which the signal is passed)
% Fpass2=2100Hz (Frequency above which the signal is passed)
% This way, the signal that lied in between 850Hz to 2100Hz was attenuated
% which included both the noise and a portion of the sinusoid between that
% range. Chebyshev Type 1 design method was used because this gave a ripple
% towards the edge of the stopband allowing some portion of the signal to
% be heard at the edges of the stopband making it less prominent that a
% small portion of sinusoid was removed along with the noise.

% The filter was found to be stable. The filter information mentioned in 
% filterDesign states the filter to be stable but it can be also verified
% by looking at the pole-zero graph as all poles lie inside the radius of
% the circle stating that the values are within the range, hence being
% stable.

Fil_Order = 6;
Fpass1 = 850;  % First Passband Frequency
Fpass2 = 2100;  % Second Passband Frequency
Apass  = 1;     % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h1  = fdesign.bandstop('N,Fp1,Fp2,Ap', Fil_Order, Fpass1, Fpass2, Apass, Fs);
Hd1 = design(h1, 'cheby1');

x_fil_chy = filter(Hd1,x); % Applying the filter on x in time domain
x_filchy_fft = fft(x_fil_chy); % performing FFT of the filtered signal

% Graphing the FFT of IIR filtered signal over normalized frequency
figure(2)
stem(2*freq./Fs,abs(x_filchy_fft))
grid on
title('FFT of IIR Filtered X using Chebyshev Method')
xlabel('Normalized Frequency (\times\pi rad/sample)') 
ylabel('Filtered X[k]')
xlim([0 1]) 

% Graphing the magnitude and frequency response of the IIR filter
freqz(Hd1)
% Graphing pole-zero plot of the IIR filter
zplane(Hd1)

% Source used:
% http://ijarcet.org/wp-content/uploads/IJARCET-VOL-4-ISSUE-6-2703-2709.pdf

%% Question 5:

% Using an equiripple design method, an FIR filter to attenuate the noise
% of the signal was designed specifying the order to be minimum. This
% resulted in the minimum order resulting to be 342 which is significantly
% higher than the filter of order 6 that we used for IIR design. The sound
% output of the filtered signal also retained less of the original signal
% comapared to the IIR design and also passed low amplitude of the noise
% unlike the IIR filter which completely removed the noise from the signal.
% Comparing the graphs for the magnitude and phase response of both the
% filters, I would prefer Chebyshev IIR filter design over Equiripple
% filter design both in terms of efficiency (less order of filter) and
% quality of output.

Fpass1 = 900;              % First Passband Frequency
Fstop1 = 1000;             % First Stopband Frequency
Fstop2 = 2100;             % Second Stopband Frequency
Fpass2 = 2200;             % Second Passband Frequency
Dpass1 = 0.057501127785;  % First Passband Ripple
Dstop  = 0.001;           % Stopband Attenuation
Dpass2 = 0.057501127785;  % Second Passband Ripple
dens   = 20;               % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[Fil_Order, Fo, Ao, W] = firpmord([Fpass1 Fstop1 Fstop2 Fpass2]/(Fs/2), [1 0 ...
                          1], [Dpass1 Dstop Dpass2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(Fil_Order, Fo, Ao, W, {dens});
Hd2 = dfilt.dffir(b);

x_fil_equ = filter(Hd2,x); % Applying the filter on x in time domain
x_filequ_fft = fft(x_fil_equ); % performing FFT of the filtered signal


% Graphing the FFT of FIR filtered signal over normalized frequency
figure(5)
stem(2*freq./Fs,abs(x_filequ_fft))
grid on
title('FFT of FIR Filtered X uxinb Equiripple Method')
xlabel('Normalized Frequency (\times\pi rad/sample)') 
ylabel('Filtered X[k]')
xlim([0 1])

% Graphing the magnitude and frequency response of the FIR filter
freqz(Hd2)
% Graphing pole-zero plot of the FIR filter
zplane(Hd2)

% Sources used:
% Dr. Tanja Karp : FIR_example.m, FIR_example_audio.m

%% Question 6:

% When the filter structure is converted to a single IIR filter of order 6,
% the magnitude response remained unchanged for both the filters but the
% phase response and pole-zero plot changed compared to that of the 2
% cascaded filter design. This change is occured because large order single
% stage filters are problematic because calculating the numerator and
% denominator from the poles and zeros is not well defined problem and is
% very unstable. Also, Group delay is always defined as the negative first 
% derivative of the phase of the transfer function vs. frequency.
% This is the reason why the phase response is shifted to negative y-axis.

Fil_Order = 6;
Fpass1 = 850;   % First Passband Frequency
Fpass2 = 2100;  % Second Passband Frequency
Apass  = 1;     % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.bandstop('N,Fp1,Fp2,Ap', Fil_Order, Fpass1, Fpass2, Apass, Fs);
Hd6 = design(h, 'cheby1');

% Get the transfer function values.
[b, a] = tf(Hd6);

% Convert to a singleton filter.
Hd6 = dfilt.df2(b, a);

% Applying the Single Section IIR filter on x in time domain
x_fil_chy1 = filter(Hd6,x); 
x_filchy1_fft = fft(x_fil_chy1); % performing FFT of the filtered signal

% Graphing the FFT of Single section IIR filtered signal over normalized frequency
figure(8)
stem(2*freq./Fs,abs(x_filchy1_fft))
grid on
title('FFT of Single section IIR filtered signal Using Chebyshev Method')
xlabel('Normalized Frequency (\times\pi rad/sample)') 
ylabel('Filtered X[k]')
xlim([0 1])

% Graphing the magnitude and frequency response of Single section IIR filter
freqz(Hd6)
% Graphing pole-zero plot of Single section IIR filter
zplane(Hd6)

% Source used:
% https://dsp.stackexchange.com/questions/22497/tradeoff-between-filter-order-and-sections

%% Question 7:

% Creating an audio file using filtered signal by IIR Chebyshev Design
% Method with Second-Order Sections
audiowrite('NoNoiseIIRCHE.wav',x_fil_chy,Fs); 

sound(x_fil_chy,Fs); % audio signal without the noise and some portion of
% sinusoid also lost using Second-Order Sections Chebyshev IIR  filter

pause(4)
% Creating an audio file using filtered signal by FIR Equiripple Design
% Method with Second-Order Sections
audiowrite('NoNoiseFIREQU.wav',x_fil_equ,Fs);

sound(x_fil_equ,Fs); % audio signal without the noise and some portion of
% sinusoid also lost using Equiripple FIR filter
pause(4)

% Creating an audio file using filtered signal by FIR Equiripple Design
% Method with Second-Order Sections
audiowrite('NoNoiseIIR1Section.wav',x_fil_chy1,Fs);

sound(x_fil_chy1,Fs); % audio signal without the noise and some portion of
% sinusoid also lost using single-section Chebyshev IIR filter
pause(4)    

% The signal was unrecoverable between 0.6 seconds to 1.4 seconds. This is
% because the noise and the increasing sinusoid co-incide with the same
% range of frequency between this time. Hence, removing the noise of this
% frequency also resulted in removal of the part of the increasing sinusoid
% of the same frequency.
% Yes, the noise for the remaining time of the signal is successfully
% removed.

%% Question 8:

% Creating a new signal between 2.5 to 3.5 seconds of the noisy_signal.wav 
% away from freq. range that co-incides with the noise

x1s = 1:(Fs+1); % array that holds 16000 values

% copying values from 2.5 to 3.5 seconds to the new array
for i=1:16000;
    x1s(i)=x(40000+i);
end

Fil_Order = 25;
Fc = 2300;  % Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h4  = fdesign.lowpass('N,F3dB', Fil_Order, Fc, Fs);
Hd4 = design(h4, 'butter');

% Applying the Single Section IIR filter on x in time domain
x1s_fil_butw = filter(Hd4,x1s); 

% Graphing the magnitude and frequency response of IIR Butterworth filter
freqz(Hd4)
% Graphing pole-zero plot of the IIR Butterworth filter
zplane(Hd4)

sound(x1s_fil_butw,Fs); % audio signal with only the noisy portion taken in
% between 2.5 to 3.5 seconds of the original signal

% Creating an audio file using filtered signal by FIR Equiripple Design
% Method with Second-Order Sections
audiowrite('OnlyNoiseIIRButw.wav',x1s_fil_butw,Fs);


% In order to remove the sinusoidal signal from the audio I used a Low-Pass
% butterworth IIR filter of order N = 25 with a cut-off frequency of 2300Hz
% because the noise lied between 950Hz to 2100Hz which is below the cut-off
% frequency ensuring that the noise lied in the passband of the filter.
% The frequency of the sinusoidal signal between 2.5 to 3.5 seconds was way
% above 2300Hz which also made sure that that the sinusoidal signal did not
% pass through the filter. This resulted in just the noisy portion being
% present when heard as sound output.
% It can also be concluded that the filter used is stable looking at the
% pole-zero plot. We can see that all poles lie within the radius of the
% circle stating that all values are bounded and absolutely summable.
% Hence, the filter is stable.