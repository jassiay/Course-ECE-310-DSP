% ECE 310 HW 1 Problem 2
% Jing Jiang

% H(z)
b1 = [3, 4];
b2 = [1, 5];
a1 = [2, -1];
a2 = [1, 3];

Hz_b = conv(b1, b2);
Hz_a = conv(a1, conv(a2, a2));
[h1, w1] = phasez(Hz_b, Hz_a, 1000);

figure();
plot(w1, h1*180/pi);
hold on

% H_min(z)
b3 = 1;
a3 = [2, -1];

[h2, w2] = phasez(b3, a3, 1000);
plot(w2, h2*180/pi);

% A_p(z)
b4 = [1, 3/4];
b5 = [1, 1/5];
b6 = [1, 3];
a4 = [1, 4/3];
a5 = [1, 5];
a6 = [1, 1/3];
Ap_b = conv(conv(b6, b6), conv(b4, b5));
Ap_a = conv(conv(a6, a6), conv(a4, a5));
[h3, w3] = phasez(Ap_b, Ap_a, 1000);
plot(w3, h3*180/pi);

xlabel('Frequency');
ylabel('Phase (deg)');
legend('H(z)','Hmin(z)','Ap(z)');



