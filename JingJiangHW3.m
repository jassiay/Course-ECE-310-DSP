% ECE 310 HW 3
% Jing Jiang

% Analog filter
wp = [100e3 130e3]*2*pi;
ws = [90e3 140e3]*2*pi;

[nb_a,wn_a] = buttord(wp,ws,2,30,'s');
[zb_a,pb_a,kb_a] = butter(nb_a,wn_a,'bandpass','s');
[nc1_a,wp1_a] = cheb1ord(wp,ws,2,30,'s');
[zc1_a,pc1_a,kc1_a] = cheby1(nc1_a,2,wp1_a,'bandpass','s');
[nc2_a,wp2_a] = cheb2ord(wp,ws,2,30,'s');
[zc2_a,pc2_a,kc2_a] = cheby2(nc2_a,30,wp2_a,'bandpass','s');
[ne_a,we_a] = ellipord(wp,ws,2,30,'s');
[ze_a,pe_a,ke_a] = ellip(ne_a,2,30,we_a,'bandpass','s');

figure('Name','Analog')
subplot(2,2,1);
zplane(zb_a,pb_a)
grid
title('Butterworth Filter')
subplot(2,2,2);
zplane(zc1_a,pc1_a)
grid
title('Cheby 1 Filter')
subplot(2,2,3);
zplane(zc2_a,pc2_a)
grid
title('Cheby 2 Filter')
subplot(2,2,4);
zplane(ze_a,pe_a)
grid
title('Elliptic Filter')
hold on;

% Digital via bilinear transform
fs = 400e3;
wp_d = [100e3 130e3]*2/fs;
ws_d = [90e3 140e3]*2/fs;
[nb_d,wn_d] = buttord(wp_d,ws_d,2,30);
[zb_d,pb_d,kb_d] = butter(nb_d,wn_d,'bandpass');
[nc1_d,wp1_d] = cheb1ord(wp_d,ws_d,2,30);
[zc1_d,pc1_d,kc1_d] = cheby1(nc1_d,2,wp1_d,'bandpass');
[nc2_d,wp2_d] = cheb2ord(wp_d,ws_d,2,30);
[zc2_d,pc2_d,kc2_d] = cheby2(nc2_d,30,wp2_d,'bandpass');
[ne_d,we_d] = ellipord(wp_d,ws_d,2,30);
[ze_d,pe_d,ke_d] = ellip(ne_d,2,30,we_d,'bandpass');

figure('Name','Digital')
subplot(2,2,1);
zplane(zb_d,pb_d)
grid
title('Butterworth Filter')
subplot(2,2,2);
zplane(zc1_d,pc1_d)
grid
title('Cheby 1 Filter')
subplot(2,2,3);
zplane(zc2_d,pc2_d)
grid
title('Cheby 2 Filter')
subplot(2,2,4);
zplane(ze_d,pe_d)
grid
title('Elliptic Filter')
hold on;

% Digital via impulse invariance
[bb_a, ab_a] = butter(nb_a,wn_a,'s');
[bc1_a, ac1_a] = cheby1(nc1_a,2,wp1_a,'s');
[bc2_a, ac2_a] = cheby2(nc2_a,30,wp2_a,'s');
[be_a,ae_a] = ellip(ne_a,2,30,we_a,'s');

[bb_d,ab_d] = impinvar(bb_a,ab_a,fs);
[bc1_d,ac1_d] = impinvar(bc1_a,ac1_a,fs);
[bc2_d,ac2_d] = impinvar(bc2_a,ac2_a,fs);
[be_d,ae_d] = impinvar(be_a,ae_a,fs);

figure('Name','Digital impulse inv.')
subplot(2,2,1);
zplane(bb_d,ab_d)
grid
title('Butterworth Filter')
subplot(2,2,2);
zplane(bc1_d,ac1_d)
grid
title('Cheby 1 Filter')
subplot(2,2,3);
zplane(bc2_d,ac2_d)
grid
title('Cheby 2 Filter')
subplot(2,2,4);
zplane(be_d,ae_d)
grid
title('Elliptic Filter')
hold on;

% Analog filter magnitude response
w = linspace(0,200e3*2*pi,2000);
hb_a = freqs(bb_a,ab_a,w);
figure('Name','Analog Magnitude Response');
subplot(2,2,1);
plot(w,20*log10(hb_a));
axis([0 200e3*2*pi -50 1]);
title('Butterworth');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

hc1_a = freqs(bc1_a,ac1_a,w);
subplot(2,2,2);
plot(w,20*log10(abs(hc1_a)));
axis([0 200e3*2*pi -50 1]);
title('Cheby 1');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

hc2_a = freqs(bc2_a,ac2_a,w);
subplot(2,2,3);
plot(w,20*log10(abs(hc2_a)));
axis([0 200e3*2*pi -50 1]);
title('Cheby 2');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

he_a = freqs(be_a,ae_a,w);
subplot(2,2,4);
plot(w,20*log10(abs(he_a)));
axis([0 200e3*2*pi -50 1]);
title('Elliptic');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

% Digital filter magnitude response
f = linspace(0, fs/2, 2000);
[bb_d, ab_d] = butter(nb_d,wn_d);
[bc1_d, ac1_d] = cheby1(nc1_d,2,wp1_d);
[bc2_d, ac2_d] = cheby2(nc2_d,30,wp2_d);
[be_d,ae_d] = ellip(ne_d,2,30,we_d);

hb_d = freqz(bb_d,ab_d,f,fs);
figure('Name','Digital Magnitude Response');
subplot(2,2,1);
plot(f,20*log10(hb_d));
axis([0 fs/2 -50 1]);
title('Butterworth');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

hc1_d = freqz(bc1_d,ac1_d,f,fs);
subplot(2,2,2);
plot(f,20*log10(hc1_d));
axis([0 fs/2 -50 1]);
title('Cheby 1');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

hc2_d = freqz(bc2_d,ac2_d,f,fs);
subplot(2,2,3);
plot(f,20*log10(hc2_d));
axis([0 fs/2 -50 1]);
title('Cheby 2');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

he_d = freqz(be_d,ae_d,f,fs);
subplot(2,2,4);
plot(f,20*log10(he_d));
axis([0 fs/2 -50 1]);
title('Elliptic ');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

% Digital via impulse invariance magnitude response
hb_di = freqz(bb_d,ab_d,f,fs);
figure('Name','Digital via Inv. Magnitude Response');
subplot(2,2,1);
plot(f,20*log10(hb_di));
axis([0 fs/2 -50 1]);
title('Butterworth');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

hc1_di = freqz(bc1_d,ac1_d,f,fs);
subplot(2,2,2);
plot(f,20*log10(hc1_di));
axis([0 fs/2 -50 1]);
title('Cheby 1');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

hc2_di = freqz(bc2_d,ac2_d,f,fs);
subplot(2,2,3);
plot(f,20*log10(hc2_di));
axis([0 fs/2 -50 1]);
title('Cheby 2');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');

he_di = freqz(be_d,ae_d,f,fs);
subplot(2,2,4);
plot(f,20*log10(he_di));
axis([0 fs/2 -50 1]);
title('Elliptic');
ylabel('Magnitude (dB)');
xlabel('Freq (Hz)');


% Orders

% Order of Analog Butterworth bandpass filter is 18
% Order of Analog Cheby I bandpass filter is 10
% Order of Analog Cheby II bandpass filter is 10
% Order of Analog Elliptic bandpass filter is 8

% Order of digital Butterworth bandpass filter is 16
% Order of digital Chebyshev I bandpass filter is 10
% Order of digital Chebyshev II bandpass filter is 10
% Order of digital Elliptic bandpass filter is 6

% Order of digital butterworth bandpass via impulse invariance 18
% Order of digital Cheby 1 bandpass via impulse invariance is 10
% Order of digital Cheby 2 bandpass via impulse invariance is 10
% Order of digital Elliptic bandpass via impulse invariance is 10


% time constants
% Analog
tb_a = 1/abs(max(real(pb_a)));
tc1_a = 1/abs(max(real(pc1_a)));
tc2_a = 1/abs(max(real(pc2_a)));
te_a = 1/abs(max(real(pe_a)));

% Digital
T = 1/fs;
tb_d = T/abs(log(max(abs(pb_d))));
tc1_d = T/abs(log(max(abs(pc1_d))));
tc2_d = T/abs(log(max(abs(pc2_d))));
te_d = T/abs(log(max(abs(pe_d))));

% Digital via impluse inv.
[zb2_d,pb2_d,kb2_d] = tf2zp(bb_d,ab_d);
[zc12_d,pc12_d,kc12_d] = tf2zp(bc1_d,ac1_d);
[zc22_d,pc22_d,kc22_d] = tf2zp(bc2_d,ac2_d);
[ze2_d,pe2_d,ke2_d] = tf2zp(be_d,ae_d);

tb_di = T/abs(log(max(abs(pb2_d))));
tc1_di = T/abs(log(max(abs(pc12_d))));
tc2_di = T/abs(log(max(abs(pc22_d))));
te_di = T/abs(log(max(abs(pe2_d))));
% Yes, the impulse invariance method preserves time constant better

% Bilinear Exploration
wp_d = [100e3 130e3]*2*pi/fs; 
ws_d = [90e3 140e3]*2*pi/fs;
o_p = tan(wp_d/2);
o_s = tan(ws_d/2);
o = sqrt(o_p(1)*o_p(2));
B = o_p(1) - o_p(2);
o_stop = min([(abs((o_s(1)^2-o^2)/(B*o))) (abs((o_s(2)^2-o^2)/(B*o)))]);