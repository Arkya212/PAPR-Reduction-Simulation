clear all
clc
close
M = 2;                          
no_of_data_points = 256;     
block_size = 128;              %Number of Subcarriers 
cp_len = ceil(0.1*block_size);  %   length of cyclic prefix 
no_of_ifft_points = block_size; 
%transmitter
rng(42);  
% Generate a random bit sequence
data_source = randi([0, 1], 1, no_of_data_points);
qam_modulated_data = qammod(data_source, M, 'PlotConstellation',true);
num_cols=length(qam_modulated_data)/block_size;
data_matrix = reshape(qam_modulated_data, block_size, num_cols); %Arranges the qam_modulated_data in size (8X16)
%   Create empty matrix to put the IFFT'd data
cp_start = block_size-cp_len;
cp_end = block_size;
%   Operate column wise & do CP
for i=1:num_cols
 ifft_data_matrix(:,i) = ifft((data_matrix(:,i)),no_of_ifft_points);
 %   Compute and append Cyclic Prefix
 for j=1:cp_len
    actual_cp(j,i) = ifft_data_matrix(j+cp_start,i);
 end
 %   Add the CP to the existing block to create the actual OFDM block
 ifft_data(:,i) = vertcat(actual_cp(:,i),ifft_data_matrix(:,i));
end
%    Convert to serial stream for transmission
[rows_ifft_data cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;
%   Actual OFDM signal to be transmitted
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);
%    clipping as a PAPR reduction method    
avg=0.6;                            %clipping levels at 0.6,.7,.8,.9
clipped=ofdm_signal;
for i=1:length(clipped)
	if (clipped(i) > avg)
		clipped(i) = avg;
 end
 if clipped(i) < -avg
		clipped(i) = -avg;
 end
end
%   HPA
%To show the effect of the PA simply we will add random complex noise
%when the power exceeds the avg. value, otherwise it adds nothing.
% 1. Generate random complex noise
noise = randn(1,len_ofdm_data) +  sqrt(-1)*randn(1,len_ofdm_data);
%   channel:
%   Create a complex AWGN channel       
channel = randn(1,block_size) + sqrt(-1)*randn(1,block_size);
%receiver   
order = 4;        % Filter order
cutoff_frequency = 0.8;  % Cutoff frequency normalized to Nyquist frequency
[b, a] = butter(order, cutoff_frequency, 'low');
%   1.  Pass the OFDM signal through the channel
after_channel = filter(b, a, ofdm_signal);
%   2.   Add Noise
awgn_noise = awgn(zeros(1,length(after_channel)),0);
%   3.  Add noise to the signal...
recvd_signal = awgn_noise+after_channel;
%   4.  Convert Data back to "parallel" form to perform FFT
recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);
%   5.  Remove CP
recvd_signal_matrix(1:cp_len,:)=[];
%   6.  Perform FFT
for i=1:cols_ifft_data
 %   FFT
 fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_ifft_points);
end
%   7.  Convert to serial stream
recvd_serial_data = reshape(fft_data_matrix, 1,(block_size*num_cols));
% Assuming recvd_serial_data is your received QAM-modulated signal
% 1. Demodulate the received QAM data
qam_demodulated_data = qamdemod(recvd_serial_data, M, 'Gray', 'UnitAveragePower', true);
%   receiver of clipped signal
%   1.  Pass the OFDM signal through the channel
after_channel = filter(channel, 1, clipped);
%   2.   Add Noise
awgn_noise = awgn(zeros(1,length(after_channel)),0);
%   3.  Add noise to the signal...
recvd_signal = awgn_noise+after_channel;
%   4.  Convert Data back to "parallel" form to perform FFT
recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);
%   5.  Remove CP
recvd_signal_matrix(1:cp_len,:)=[];
%   6.  Perform FFT
for i=1:cols_ifft_data
 %   FFT
 fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_ifft_points);
end
%   7.  Convert to serial stream
recvd_serial_data = reshape(fft_data_matrix, 1,(block_size*num_cols));
qam_demodulated_data = qamdemod(recvd_serial_data, M);
% % Calculate PAPR for QAM-demodulated data
power_peak = abs(recvd_serial_data).^2;
power_mean = mean(abs(recvd_serial_data).^2);
PAPR = 10 * log10(power_peak./power_mean);
fprintf('PAPR after QAM Input: %.2f dB\n', PAPR);
%CCDF
[f_cdf, x_cdf] = ecdf(PAPR);
f_ccdf = 1-f_cdf;
% Plot CCDF of PAPR
figure(12);
semilogy(x_cdf, f_ccdf, 'LineWidth', 2);
grid on;
title('CCDF of PAPR for QAM-demodulated Signal');
xlabel('PAPR (dB)');
ylabel('Probability');

