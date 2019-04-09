% Audio Signal CWK Ajil %
T = 1; % Period (seconds) window size 
fs = 10000; % Sample Frequency (Hz)
df= 1/T; % Frequency Resolution 
dt= 1/fs; % Time resolution
t=(0:dt:T-dt); % Time vector
Amp=1;
s_rate = 200; % Swept rate (Hz/s) how condensed
f0=0; % Initial value pitch control positive start of figure 2 starting frequency
f=f0+s_rate*t; % Array of values
sweptsin=Amp.*sin(2*pi*f.*t);% Multiply arr ays to to avoid using loops

% Plot for sweptsin wave
figure(1)
plot(t,sweptsin)
title('Swept Sine Signal')
xlabel('Time (s)')
ylabel('Amplitude')

% Change sweptsin: amplitude
figure(2)
subplot(2,1,1);
Amp=5;
sweptsin1=Amp.*sin(2*pi*f.*t);% Multiply arr ays to to avoid using loops
plot(t,sweptsin1)
title('Swept Sine Signal with changed amplitude')
xlabel('Time (s)')
ylabel('Amplitude')

% Change sweptsin:log
subplot(2,1,2);
sweptsin2=log(sweptsin);% Multiply arr ays to to avoid using loops
plot(t,sweptsin2)
title('Swept Sine Signal using log')
xlabel('Time (s)')
ylabel('Amplitude')

%Ammend values
T = 5; % Period (seconds) window size 
fs = 10000; % Sample Frequency (Hz)
df= 1/T; % Frequency Resolution 
dt= 1/fs; % Time resolution
t=(0:dt:T-dt); % Time vector
Amp=1;
s_rate = 200; % Swept rate (Hz/s) how condensed
f0=0; % Initial value pitch control positive start of figure 2 starting frequency
f=f0+s_rate*t; % Array of values
sweptsin3=Amp.*sin(2*pi*f.*t);% Multiply arr ays to to avoid using loops

% Change sweptsin: time vector
figure(3)
plot(t,sweptsin3)
title('Swept Sine Signal')
xlabel('Time (s)')
ylabel('Amplitude')

% Audio Signal CWK%
T = 1; % Period (seconds) window size 
fs = 10000; % Sample Frequency (Hz)
df= 1/T; % Frequency Resolution 
dt= 1/fs; % Time resolution
t=(0:dt:T-dt); % Time vector
Amp=1;
s_rate = 200; % Swept rate (Hz/s) how condensed
f0=0; % Initial value pitch control positive start of figure 2 starting frequency
f=f0+s_rate*t; % Array of values
sweptsin=Amp.*sin(2*pi*f.*t);% Multiply arr ays to to avoid using loops

% Plot for frequency: f(0)=0
figure(4)
plot(t,f)
title('Graph showing frequency against time')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

% Create .wav file and find bitdepth
sound(sweptsin,fs);
audiowrite('wave44100.wav',sweptsin,fs);  
[y,fs] = audioread('wave44100.wav'); % Reads the file created retrieving fs and data
info=audioinfo('wave44100.wav'); % Shows bit depth and sample frequency
figure(5)
T = 1; % Period (seconds) window size 
t=(0:dt:T-dt); % Time vector
plot(t,y) % Plot time signiture
xlabel('Time(s)')
ylabel('Amplitude');
title('Graph showing time signiture against time');

% Fourrier transform double sided
yyf=fft(sweptsin);
yf1=imag(yyf); % Finds imaginary numbers
yf2=real(yyf); % Finds real numbers
yf3=abs(yyf); % Finds magnitude 
N=length(sweptsin); 
freq_array =(0:N-1)*df; % Frequency array using frequency resolution 

% Plot imaginary
figure(6)
plot(freq_array,yf1);
xlabel('Frequency')
ylabel('Magnitude');
title('FFT using imaginary numbers');

% Plot real
figure(7)
plot(freq_array,yf2);
xlabel('Frequency')
ylabel('Magnitude');
title('FFT using real numbers');

% Plot magnitude
figure(8)
plot(freq_array,yf3);
xlabel('Frequency')
ylabel('Magnitude');
title('FFT using magnitude');

% Compare real and imaginary
figure(9)
plot(freq_array,yf1,'*g'); % Green colour plot
hold on; % Tells program to wait before plotting
plot(freq_array,yf2,'*b'); % Blue colour plot
xlabel('Frequency')
ylabel('Magnitude');
title('FFT comparing imaginary and real numbers');

% Plot fft with increasing window size for imaginary
size=length(sweptsin);
for i=1:size/5:size+1 % Splits the window size into 5
    yyf=fft(sweptsin,i); % Take fft of values upto window size
    yyfi=imag(yyf);
    Nf=numel(fft(sweptsin,i));
    i_freq_array =[0 : Nf-1]*df;
    figure(10)
    plot(i_freq_array,yyfi);
    xlabel('Frequency')
    ylabel('Magnitude');
    title('Increasing FFT window size using imaginary');
    pause(1);
end
 
 % Plot fft with increasing window size for chirp in time domain
 for T = 0:1:5 % Period (seconds)
    fs = 4100; % Sample Frequency (Hz)
    df= 1/T; % Frequency Resolution 
    dt= 1/fs; % Time resolution
    t=(0:dt:T-dt); % Time Vector
    s_rate = 200; % Swept rate (Hz/s) how condensed
    f0=0;
    f=f0+s_rate*t;
    sweptsin=sin(2*pi*f.*t);
    fft_r=fft(sweptsin); % Take fft for frequency domain 
    sweptsin_r=ifft(real(fft_r)); % Take ifft for time domian in real domain
    figure(13)
    plot(t,sweptsin_r);
    xlabel('Time (s)')
    ylabel('Magnitude');
    title('Increasing FFT window size using real time domain');
    pause(1); 
 end
%Press ctrl + enter to run program and it willdisplay all plots.


