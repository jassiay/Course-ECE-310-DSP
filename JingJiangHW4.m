% ECE 310 HW 4
% Jing Jiang

% 2
% a)
f = 10 * 10e3;
f_s = 50 * 10e3;
N = 256;
bin_sp = f_s/N;

% b)
k1 = f/bin_sp;
k2 = (f_s - f)/bin_sp;

% c)
bin_sp = bin_sp/10e3;
w = abs((2 * pi * (10 - bin_sp * 51)) / 50);
diric_0 = diric(0, 250);
diric_w = diric(w, 250);
strad_loss = 20 * log10(abs(diric_w./ diric_0))

% d)
hamming = hamming(250);
w_dc = sum(hamming);
w_d = sum(transpose(hamming) .* exp(-j * w * (1:250)));
ham_loss = 20 * log10(abs(w_d / w_dc))

% e)
N_n = 512;
bin_sp_n = f_s/N_n / 10e3;

w_n = abs((2 * pi * (10 - bin_sp_n * 51)) / 50);
diric_w_n = diric(w_n, 250);
strad_loss_n = 20 * log10(abs(diric_w_n./diric_0))

w_dc_n = sum(hamming);
w_d_n = sum(transpose(hamming) .* exp(-j * w * (1:250)));
ham_loss_n = 20 * log10(abs(w_d_n/ w_dc_n))

% 6
fs = 100e6;
t = (0:999)*1/fs;
f = 20e6;
sine = 2*sin(2*pi*f*t);
awgn = sqrt(0.2)*randn(size(t));
sig = sine + awgn;
w = chebwin(1000,30);
win_output = sig.*w';
k = -1024/2:(1024/2-1);
f2 = k * fs/1024;

fast_f = fft(win_output, 1024);
plot(f2, 20*log10(fftshift(fast_f)));
xlim([-fs/2 fs/2]);
title('FFT of Sine wave');
xlabel('Freq(Hz)');
ylabel('Magnitude(dB)');

% 7
% a)
w0 = zeros(1,1000);
for n = 0:249
    w0(4*n+3) = 1;
end

x_hat = (w0.*sig);

x_hat_fastf = fft(x_hat, 1024);
x_hat2 = 20*log10(fftshift(x_hat_fastf));

figure;
plot(f2, x_hat2);
xlim([-fs/2 fs/2]);
title('FFT of x_hat');
xlabel('Freq(Hz)');
ylabel('Magnitude(dB)');

% Peak frequencies appears at -45, -30, -20, -5, 5, 20, 30, 45 (MHz)

% b)
w0_fft = fft(w0, 1000);
figure;
stem(0:999, abs(w0_fft));
title('FFT of w0');
xlabel('k');
ylabel('Magnitude');

% w0 convolves with the original sine wave.