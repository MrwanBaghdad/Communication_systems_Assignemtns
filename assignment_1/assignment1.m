close all;

%% Read audio signal 
[signal,Fs] = audioread('malcolm_x.wav');



%% sample sound to 8000 sample/sec
sampled=resample(signal,8000,Fs);
t=length(sampled);
t=linspace(0,15,t);
plot(t,sampled);
title('sampled signal')
%% plot spectrum
spectrum=fftshift(fft(sampled));
dF = Fs/length(sampled);
f = -Fs/2:dF:Fs/2-dF; 
figure;
plot(f,abs(spectrum));
title('Spectrum')
xlabel('Frequency')
%% Uniform quantization with 256 levels.
K=8;
y=sampled(:,1);
ymax=max(y);
ymin=min(y);  
q=(ymax-ymin)/(2^K);
quant= uencode(y, 8) ;
delta = ((max(y) - min(y)) / (2^(K) -1) )
quant = uencode(y,8);
quant = double(quant);
quantnoise=mean((y-quant).^2);
fprintf('\n mean square error for uniform quantization = %g',quantnoise)
figure;
plot(t,y)
title('Uniform quantization')
  %% A-law compandor
  A1=[];
  A1 = compand(y,10,ymax,'A/compressor');
  A1=A1';
  quantnoise(1)=mean((y-A1).^2);
  fprintf('\n mean square error for A-law compandor (A=10) = %g',quantnoise(1))
  A2=[];
  A2 = compand(y,87.6,ymax,'A/compressor');
  A2=A2';
  quantnoise(2)=mean((y-A2).^2);
  fprintf('\n mean square error for A-law compandor (A=87.6) = %g',quantnoise(2))
    A3=[];
  A3 = compand(y,1000,ymax,'A/compressor');
  A3=A3';
  quantnoise(3)=mean((y-A3).^2);
  fprintf('\n mean square error for A-law compandor (A=1000) = %g',quantnoise(3))
 A=[10 87.6 1000];
 figure;
 plot(A,quantnoise)
 title('values of A vs mean square error')
 xlabel('A values')
 ylabel('mean square error')
 %% U-law compandor
   U1=[];
  U1 = compand(y,10,ymax,'mu/compressor');
  quantnoise=[];
  quantnoise(1)=mean((y-U1).^2);
  fprintf('\n mean square error for U-law compandor (U=10) = %g',quantnoise(1))
  U2=[];
  U2 = compand(y,255,ymax,'mu/compressor');
  quantnoise(2)=mean((y-U2).^2);
  fprintf('\n mean square error for U-law compandor (U=255) = %g',quantnoise(2))
    U3=[];
  U3 = compand(y,1000,ymax,'mu/compressor');
  quantnoise(3)=mean((y-U3).^2);
  fprintf('\n mean square error for U-law compandor (U=1000) = %g',quantnoise(3))
 U=[10 87.6 1000];
 figure;
 plot(U,quantnoise)
 title('values of U vs mean square error')
 xlabel('U values')
 ylabel('mean square error')

 %% Quantization of A2 and adding noise 

pcm = uencode(A2, 8);
SNR = 10:2:30;
ber = []    
i = 1

 for i= 1:length(SNR)
received = (qamdemod(awgn(qammod(pcm, 256, 0,'bin'), SNR(i),'measured'), 256, 0,'bin'))
[n, b] = biterr(de2bi(pcm), de2bi(received));
ber(i) = b;
%% reconstruction 
% received signal is quantized multiply by delta
reconstruction = udecode(uint16(received), 8)
reconstruction_expanded = compand(reconstruction, 87.6, max(y), 'a/expander')
sound(reconstruction_expanded)
file_name = strcat(num2str(i),'.wav')
audioWrite(file_name, reconstruction_expanded, Fs)
 end
plot(SNR, ber)