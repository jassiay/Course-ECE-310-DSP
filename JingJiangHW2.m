% ECE 310 HW 2
% Jing Jiang

% Problem 1
% a)

wp_hi = 130000 * 2 * pi;
wp_lo = 100000 * 2 * pi;
ws_hi = 140000 * 2 * pi;
ws_lo = 90000 * 2 * pi;

B = wp_hi - wp_lo;
w0 = sqrt(wp_hi * wp_lo);

w_proto1 = abs((ws_lo^2 - w0^2) / (B * ws_lo));
w_proto2 = abs((ws_hi^2 - w0^2) / (B * ws_hi));

% b)
rp = 2;
rs = 30;
nButt = 1/2 * log10((10^(rs/10) - 1) / (10^(rp/10) - 1)) / log10(w_proto2);
nCheb = acosh(sqrt((10^(rs/10) - 1) / (10^(rp/10) - 1)))/acosh(w_proto2);

% c)
wp = [wp_lo wp_hi];
ws = [ws_lo ws_hi];

[nButt, wnButt] = buttord(wp, ws, rp, rs,'s');
[zButt, pButt, kButt] = butter(nButt, wnButt, 's');
[nCheb1, wpCheb1] = cheb1ord(wp, ws, rp, rs,'s');
[zCheb1, pCheb1, kCheb1] = cheby1(nCheb1, rp, wpCheb1, 's');
[nCheb2, wpCheb2] = cheb2ord(wp, ws, rp, rs,'s');
[zCheb2, pCheb2, kCheb2] = cheby2(nCheb2, rs, wpCheb2, 's');
[nElp, wpElp] = ellipord(wp, ws, rp, rs,'s');
[zElp, pElp, kElp] = ellip(nElp, rp, rs, wpElp, 's');

% d)
figure('Name','1. (d)')
subplot(2,2,1);
zplane(zButt, pButt)
grid
title('Butterworth Bandpass Digital Filter')

subplot(2,2,2);
zplane(zCheb1, pCheb1)
grid
title('Cheby 1 Bandpass Digital Filtr')

subplot(2,2,3);
zplane(zCheb2, pCheb2)
grid
title('Cheby 2 Bandpass Digital Filter')

subplot(2,2,4);
zplane(zElp, pElp)
grid
title('Ellip Bandpass Digital Filter')

% e)
[nElp2,wpElp2] = ellipord(1, w_proto2, rp, rs,'s');
[zElp2,pElp2,kElp2] = ellip(nElp, rp, rs, wpElp, 's');

figure('Name','1. (e)')
zplane(zElp2,pElp2)
grid
title('Ellip Lowpass Digital Filter')

% f)
w = linspace(0, 200e3*2*pi, 2000);

% Butterworth
[b_bw,a_bw] = butter(nButt,wnButt, 's');
h_bw = freqs(b_bw,a_bw,w);
figure('Name','1. (f) Butterworth');
subplot(2,1,1);
plot(w,db(h_bw))
hold on
hline = refline(0,0);
hline.Color = 'r';
hline = refline(0,-2);
hline.Color = 'r';
hline = refline(0, -30);
hline.Color = 'r';
hold off
axis([0 200e3*2*pi -100 20])
yticks(-100:10:20);
title('Butterworth Magnitude')
ylabel('Magnitude (dB)')
xlabel('W (rad/s)');
subplot(2,1,2);
plot(w,unwrap(angle(h_bw))*180/pi);
title('Butterworth Phase');
ylabel('Phase (deg)')
xlabel('W (rad/s)');

% Cheby 1
[b_c1,a_c1] = cheby1(nCheb1, rp, wpCheb1, 's');
h_c1 = freqs(b_c1,a_c1,w);
figure('Name','1. (f) Cheby 1');
subplot(2,1,1);
plot(w,db(h_c1))
hold on
hline = refline(0,0);
hline.Color = 'r';
hline = refline(0,-2);
hline.Color = 'r';
hline = refline(0, -30);
hline.Color = 'r';
hold off
axis([0 200e3*2*pi -100 20])
yticks(-100:10:20);
title('Cheby 1 Magnitude')
ylabel('Magnitude (dB)')
xlabel('W (rad/s)');
subplot(2,1,2);
plot(w,unwrap(angle(h_c1))*180/pi);
title('Cheby 1 Phase');
ylabel('Phase (deg)')
xlabel('W (rad/s)');

% Cheby 2
[b_c2,a_c2] = cheby2(nCheb2, rs, wpCheb2, 's');
h_c2 = freqs(b_c2,a_c2,w);
figure('Name','1. (f) Cheby 2');
subplot(2,1,1);
plot(w,db(h_c2))
hold on
hline = refline(0,0);
hline.Color = 'r';
hline = refline(0,-2);
hline.Color = 'r';
hline = refline(0, -30);
hline.Color = 'r';
hold off
axis([0 200e3*2*pi -100 20])
yticks(-100:10:20);
title('Cheby 2 Magnitude')
ylabel('Magnitude (dB)')
xlabel('W (rad/s)');
subplot(2,1,2);
plot(w,unwrap(angle(h_c2))*180/pi);
title('Cheby 2 Phase');
ylabel('Phase (deg)')
xlabel('W (rad/s)');

% Ellip
[b_el,a_el] = ellip(nElp, rp, rs, wpElp, 's');
h_el = freqs(b_el,a_el,w);
figure('Name','1. (f) Elliptic');
subplot(2,1,1);
plot(w,db(h_el))
hold on
hline = refline(0,0);
hline.Color = 'r';
hline = refline(0,-2);
hline.Color = 'r';
hline = refline(0, -30);
hline.Color = 'r';
hold off
axis([0 200e3*2*pi -100 20])

yticks(-100:10:20);
title('Ellip Magnitude')
ylabel('Magnitude (dB)')
xlabel('W (rad/s)');
subplot(2,1,2);
plot(w,unwrap(angle(h_el))*180/pi);
title('Ellip Phase');
ylabel('Phase (deg)')
xlabel('W (rad/s)');

% g)
figure('Name','1. (g) Elliptic');
plot(w,db(h_bw));
hold on;
plot(w,db(h_c1));
hold on;
plot(w,db(h_c2));
hold on;
plot(w,db(h_el));
axis([0 200e4 -100 20]);
yticks(-100:10:20);
legend('butter','cheby1','cheby2','ellip');

% h)
ws_edge = [ws_lo, ws_hi];
att_bw = db(freqs(b_bw,a_bw,ws_edge));
att_c1 = db(freqs(b_c1, a_c1, ws_edge));
att_c2 = db(freqs(b_c2,a_c2,ws_edge));
att_ell = db(freqs(b_el,a_el,ws_edge));

filter = {'Butterworth';'Cheby 1';'Cheby 2';'Elliptic'};
attenuation = [att_bw;att_c1;att_c2;att_ell];
table = table(filter,attenuation)

