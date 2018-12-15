% Jing Jiang
% HW8

% 4
% a)
order = 8;
rp = 1.5;
rs = 30;
[eb, ea] = ellip(4, rp,rs,[0.3 0.6]);
[he,we] = freqz(eb,ea, 800);

figure;
plot(we, 20*log10(abs(he)));
title('Elliptic Filter');
xlabel('Freq');
ylabel('Magnitude');

% b)
[z,p,k] = tf2zpk(eb,ea);
[sos_up, g_up] = zp2sos(z,p,k,'up', 'inf');
[sos_down, g_down] = zp2sos(z,p,k,'down', 'inf');

% c)
% up
[z1,p1,k1] = sos2zp(sos_up(1,:));
[z2,p2,k2] = sos2zp(sos_up(2,:));
[z3,p3,k3] = sos2zp(sos_up(3,:));
[z4,p4,k4] = sos2zp(sos_up(4,:));
up_zeros = [abs(z1);abs(z2);abs(z3);abs(z4)]; 
up_poles = [abs(p1);abs(p2);abs(p3);abs(p4)]; 
% magnitude of zeros are 1
% poles are indeed in increasing magnitude

% down
[z5,p5,k5] = sos2zp(sos_down(1,:));
[z6,p6,k6] = sos2zp(sos_down(2,:));
[z7,p7,k7] = sos2zp(sos_down(3,:));
[z8,p8,k8] = sos2zp(sos_down(4,:));
down_zeros = [abs(z5); abs(z6); abs(z7);abs(z8)]; 
down_poles = [abs(p5);abs(p6);abs(p7);abs(p8)]; 
% magnitude of zeros are 1
% poles are indeed in decreasing magnitude

% d)
% up
for i = 1:5
    if i == 1
        d1{1} = g_up/filt(sos_up(i,4:end), 1);
    elseif i ==5
        d1{5} = filt(sos_up(i-1,1:3), 1)*d1{i-1};
    else
        d1{i} = d1{i-1}*(filt(sos_up(i-1,1:3), 1)/filt(sos_up(i,4:end), 1));
    end
end

[b0, a0] = tfdata(d1{1});
[b1, a1] = tfdata(d1{2});
[b2, a2] = tfdata(d1{3});
[b3, a3] = tfdata(d1{4});

[h0, w] = freqz(b0{1}, a0{1}, 800);
[h1, w] = freqz(b1{1}, a1{1}, 800);
[h2, w] = freqz(b2{1}, a2{1}, 800);
[h3, w] = freqz(b3{1}, a3{1}, 800);

figure;
subplot(2,1,1)
hold on
plot(w/pi, 20*log10(abs(h0)));
plot(w/pi, 20*log10(abs(h1)));
plot(w/pi, 20*log10(abs(h2)));
plot(w/pi, 20*log10(abs(h3)));
hold off
title('Up SOS');
xlabel('Freq');
ylabel('Magnitude');

% down
for i = 1:5
    if i == 1
        d2{1} = g_down/filt(sos_down(i,4:end), 1);
    elseif i == 5
        d2{5} = filt(sos_down(4,1:3), 1)*d2{4};
    else
        d2{i} = d2{i-1}*(filt(sos_down(i-1,1:3), 1)/filt(sos_down(i,4:end), 1));
    end
end

[b5, a5] = tfdata(d2{1});
[b6, a6] = tfdata(d2{2});
[b7, a7] = tfdata(d2{3});
[b8, a8] = tfdata(d2{4});

[h5, w] = freqz(b5{1}, a5{1}, 800);
[h6, w] = freqz(b6{1}, a6{1}, 800);
[h7, w] = freqz(b7{1}, a7{1}, 800);
[h8, w] = freqz(b8{1}, a8{1}, 800);

subplot(2,1,2)
hold on
plot(w/pi, 20*log10(abs(h5)));
plot(w/pi, 20*log10(abs(h6)));
plot(w/pi, 20*log10(abs(h7)));
plot(w/pi, 20*log10(abs(h8)));
hold off
title('Down SOS');
xlabel('Freq');
ylabel('Magnitude');

% e)
scalar_4e = sos_up./flipud(sos_down)

% 5
% a)
bf = [0.1336, 0.0568, 0.0563, 0.1336];
af = [1, -1.5055, 1.2630, -0.3778];
b1 = [-0.4954, 1];
a1 = [1, -0.4954];
b2 = [0.7632, -1.0101, 1];
a2 = [1, -1.0101, 0.7632];

bf_q = fi(bf,1,5,3);
bf_q = bf_q.data;
af_q = fi(af,1,5,3);
af_q = af_q.data;

b1_q = fi(b1,1,5,3);
b1_q = b1_q.data;
a1_q = fi(a1,1,5,3);
a1_q = a1_q.data;
b2_q = fi(b2,1,5,3);
b2_q = b2_q.data;
a2_q = fi(a2,1,5,3);
a2_q = a2_q.data;

% b)
Hf_p = sum(bf) / sum(af);
Hf_n = sum(bf .* [1, -1, 1, -1]) / sum(af .* [1, -1, 1, -1]);
Hf_pq = sum(bf_q) / sum(af_q);
Hf_nq = sum(bf_q .* [1, -1, 1, -1]) / sum(af_q .* [1, -1, 1, -1]);
Ha_pq = .5*(sum(b1_q)/sum(a1_q) + sum(b2_q)/sum(a2_q));
Ha_nq = .5*(sum(b1_q.*[1, -1])/sum(a1_q.*[1, -1]) + 1);

Hf_z = 20*log10(abs(Hf_p));
Hf_pi = 20*log10(abs(Hf_n));
Hf_q_z = 20*log10(abs(Hf_pq));
Hf_q_pi = 20*log10(abs(Hf_nq));
Ha_q_z = 20*log10(abs(Ha_pq));
Ha_q_pi = 20*log10(abs(Ha_nq));

% Report the Errors
Hf_z_err = Hf_z - Hf_q_z
Hf_pi_err = Hf_pi - Hf_q_pi
Ha_z_err = Hf_z - Ha_q_z
Ha_pi_err = Hf_pi - Ha_q_pi

% c)
[Hf, w] = freqz(bf, af, 1e4);
[Hf_q, w] = freqz(bf_q, af_q, 1e4);
[H1_q, w] = freqz(b1_q, a1_q, 1e4);
[H2_q, w] = freqz(b2_q, a2_q, 1e4);
Hf_dB = 20*log10(abs(Hf));
Hf_q_dB = 20*log10(abs(Hf_q));
Ha_q_dB = 20*log10(abs(0.5*(H1_q + H2_q)));
figure;
hold on
plot(w, Hf_dB)
plot(w, Hf_q_dB)
plot(w, Ha_q_dB)
legend('Original', 'Quantized', 'Quantized Allpass sum')
xlabel('Freq')
xlim([0 pi])
ylabel('dB')
title('5c Comparison of magnitude responses')
ylim([-40 0])

% d)
% 1
max_dev_pb_q = max(abs(Hf_dB(1:3001) - Hf_q_dB(1:3001)))
max_dev_pb_a = max(abs(Hf_dB(1:3001) - Ha_q_dB(1:3001)))

% 2
ripple_pb = max(Hf_dB(1:3001)) - min(Hf_dB(1:3001))
ripple_pb_a = max(Ha_q_dB(1:3001)) - min(Ha_q_dB(1:3001))
ripple_sb = max(Hf_dB(4101:end))
ripple_sb_a = abs(max(Ha_q_dB(4101:end)))

% 3
max_gain_sb = max(Hf_dB(4101:end))
max_gain_sb_q = max(Hf_q_dB(4101:end))
max_gain_sb_a = max(Ha_q_dB(4101:end))

% e)
fa = .5*(tf(b1_q,a1_q) + tf(b2_q, a2_q));
ba = fa.Numerator{1,1};
aa = fa.Denominator{1,1};

[gd, w] = grpdelay(bf, af, 1e4);
[gd_q, w] = grpdelay(bf_q, af_q, 1e4);
[gd_a, w] = grpdelay(ba, aa, 1e4);

figure;
hold on
plot(w, gd)
plot(w, gd_q)
plot(w, gd_a)
legend('Original', 'Quantized', 'Quantized Allpass sum')
xlabel('Freq')
ylabel('Samples')
ylim([-20 20])
title('Group Delay Comparison')
% The group delay of the parallel allpass realization is better than the original

% f)
[z, p, k] = tf2zp(bf, af);
% report zeros and poles of the original filter
z
p
[zq, pq, kq] = tf2zp(bf_q, af_q);
% report zeros and poles of the quantized original filter
zq
pq
[za, pa, ka] = tf2zp(ba, aa);
% report zeros and poles of the quantized sum of allpass filter
za
pa
% Poles don't move outside the unit circle
% The zero move to the unit circle -1+0j



