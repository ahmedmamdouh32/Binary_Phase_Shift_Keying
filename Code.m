clear
clc
close all

Rb = 100;
Fc = 10 * Rb; %carrier frequency
Tb = 1/Rb;
Fs = 10 * Fc; %sampling frequency
Ts = 1/Fs; %sample time
nb = 1e6; %number of bits
t = 0:Ts:Tb; %total time 
plot_sample = 500; %how many samples we will draw 

phi = cos(2*pi*Fc*t);
data = randperm(nb);        %random order of numbers [1:1e6]
data_binary = mod(data,2);  %converting numbers to Zeros and Ones

data_binary_BaseBand = 2*data_binary -1;  %1--->1 , 0--->-1

data_binary_plot = complex(data_binary_BaseBand);
plot(real(data_binary_plot),imag(data_binary_plot),'ro','LineWidth',3);
hold on; % to save figure
grid on %to add grids to image
xlabel('Real');
ylabel('Imaginary');
title('Constellation Diagram');

[Row,Column] = size(data_binary_BaseBand);  %take care!!!!!!

data_binary_BaseBand_rep = repelem(data_binary_BaseBand,Row,size(phi,2)); %to give rect shape for data

phi_rep = repmat(phi,Row,Column); %repeate phi to number of data bits

modulated_signal = data_binary_BaseBand_rep.* phi_rep; %creating the modulated signal

t_total = (0:Ts:plot_sample*Tb);

figure; %creating a figure

%Drawing the carrier:
subplot(3, 1, 1);
plot(t_total(1:plot_sample), phi_rep(1:plot_sample), 'b', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Carrier Signal');

%Drawing the Data Signal:
subplot(3, 1, 2);
plot(t_total(1:plot_sample), data_binary_BaseBand_rep(1:plot_sample), 'k', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Baseband Signal');

%Drawing the Modulated Signal:
subplot(3, 1, 3);
plot(t_total(1:plot_sample), modulated_signal(1:plot_sample), 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Modulated Signal');


transmitted_signal = modulated_signal(1:plot_sample);  %transmitting a small sample of the modulated signal

samples = length(transmitted_signal); % Number of samples (5 samples)

f_axis = (-Fs/2:Fs/samples:(Fs/2-Fs/samples)); % Frequency response axis 

transmitted_signal_fft = fftshift(abs(fft(transmitted_signal)) / samples);


figure;  %creating a igure
plot(f_axis, transmitted_signal_fft, 'r', 'LineWidth', 1);
grid on;
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Spectrum');
xlim([-4*Fc, 4*Fc]); %frequency plotting range



Uniform_Noise = randn(1,size(modulated_signal,2)); %generating a uniform noise with same size of modulated signal

SNR_dB = 20; %SNR in dB
SNR = 10^(SNR_dB/10);
Signal_Power = var(modulated_signal); %variance of the modulated signal = E

Noise_Power = Signal_Power/SNR;

Noise_Equation = sqrt(Noise_Power/2)*Uniform_Noise;

Signal_Received = modulated_signal + Noise_Equation;

%calculating noise in theoretical & practical:

SNR_Values = 0:1:20;
P_Error_thr = zeros(1,length(SNR_Values)); %probability of error theorotical
P_Error_prc = ones(1,length(SNR_Values)); %probability of error practical
SNR_Real = 10.^(SNR_Values / 10); % Element-wise division and power operation to get SNR values
for i = 1:length(SNR_Values)

      SNR = 10^(SNR_Values(i)/10); % Converting SNR from dB to linear

      Signal_Power = var(modulated_signal);
      Noise_Power = Signal_Power/SNR;
      Noise_Equation = sqrt(Noise_Power/2)*Uniform_Noise;
      Signal_Received = data_binary_BaseBand_rep + Noise_Equation;
      Decision_Maker = sign(real(Signal_Received)); % 
      error = sum(Decision_Maker ~= data_binary_BaseBand_rep); % Count errors
      P_Error_prc(i) = error / nb; % Simulated BER  

end

for i = 1 : length(SNR_Values)
      SNR = 10^(SNR_Values(i)/10); % Converting SNR from dB to linear
      P_Error_thr(i) = (1/2)*erfc((sqrt(SNR)));
end

%drawing the difference between practical and theoretical
figure;
semilogy(SNR_Values, P_Error_thr, 'g-o');
xlabel('SNR dB');
ylabel('Bit Error Rate');
hold on;
semilogy(SNR_Values, P_Error_prc, 'r-o');
grid on;
legend('Theoretical','practical');
