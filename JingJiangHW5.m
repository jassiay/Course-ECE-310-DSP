% ECE 310 HW 5
% Jing Jiang

% 2
% b)
w1= linspace(-pi,pi,100);
w2= linspace(-pi,pi,100);
h = 1/6 .* [1,4,1;4,-20,4;1,4,1];
H = freqz2(h,w1,w2);
max_hf_ima = max(abs(imag(H(:))));
max_hf_ima

H_real = real(H);
[ww1,ww2] = meshgrid(w1,w2);
a = 1/6 .* (-20 + 8*cos(ww1) + 8*cos(ww2) + 2*cos(ww1+ww2) + 2*cos(ww1-ww2));
diff = abs(H_real-a);
max_diff = max(diff(:));
max_diff

% c)
figure('Name','2.c contour')
contour(w1,w2,H_real)
hold on
plot([-1,1],[-1,-1],'k',[-1,1],[1,1],'k',[-1,-1],[-1,1],...
    'k',[1,1],[-1,1],'k');

title('Contour plot')
xlabel('w1')
ylabel('w2')

figure('Name','2.c surface')
surf(w1,w2,H_real)

title('Surface plot')
xlabel('w1')
ylabel('w2')

% it appears to be a bandpass filter

% d) 4
figure('Name','2.d4')
hold on
plot(w1,H(:,51));
plot(w1,-w1.^2);

title('Examination')
xlabel('w1')
ylabel('w2')

% e)
[H, f1, f2] = freqz2(h);
[ff1,ff2] = meshgrid(f1,f2);
lapacian2 = -(ff1.^2 + ff2.^2) * pi;

figure('Name','2.e');
surf(f1,f2,lapacian2);
hold on
surf(f1,f2,H);

xlabel('f1')
ylabel('f2')

% 3
% a)

% at the bottom of the script I defined the insertion function

% b)
H2 = insertion(h);
[H3,f1,f2] = freqz2(H2);

figure('Name','3.b')
contour(f1,f2,H3)
title('Contour plot')
xlabel('f1')
ylabel('f2')

figure('Name','3.b')
surf(f1,f2,H3);
title('Surface plot')
xlabel('f1')
ylabel('f2')

% This is a plausible approximation to the Lapacian operator.

% 4
% a)
h_y = 1/8 .* [-1,-2,-1;0,0,0;1,2,1];
[dgx,f1,f2] = freqz2(h_y');
figure('Name','4.a')
surf(f1,f2,abs(dgx))

title('dg/dx Surface plot')
xlabel('x')
ylabel('y')

[dgy,f1,f2] = freqz2(h_y);
figure('Name','4.a')
surf(f1,f2,abs(dgy))

title('dg/dy Surface plot')
xlabel('x')
ylabel('y')

% b)

% at the bottom of the script I defined the edgeDetection function

% c)
circuit = imread('circuit.tif');
figure('Name','circuit')

image(circuit);
image_d = im2double(circuit);
edg = edgeDetection(image_d,0.015, 1);
imtool(edg);
% Threshold value is 0.015

median(image_d(:));
edg2 = edgeDetection(image_d, median(image_d(:)), 1);
imtool(edg2);
% Median value too big

% d)
image_d = im2double(circuit);
edg_d = edgeDetection(image_d,0.001, 0);
imtool(edg_d);
% Threshold value is 0.001

median(image_d(:));
edg2_d = edgeDetection(image_d, median(image_d(:)), 0);
imtool(edg2_d);
% Median value too big

% 3 a) insertion functoin
function m = insertion(x)
    m = zeros(size(x,1)*2, size(x,2)*2);
    for row = 1 : size(x,1)
        for col = 1 : size(x,2)
            m(row*2, col*2) = x(row, col);
        end
    end
    
    m(1,:) = [];
    m(:,1) = [];
end

% 4 b) edgeDetection functoin
function m = edgeDetection(G, T, n)
    % n is the choice between L2 norm if n=1 or L1 norm if n=0
    h_y2 = 1/8 .* [-1,-2,-1;0,0,0;1,2,1];
    conv_x = conv2(G, h_y2', 'same');
    conv_y = conv2(G, h_y2, 'same');
    if n == 1
        result = (conv_x.^2 + conv_y.^2).^0.5;
    elseif n == 0
        result = conv_x.^2 + conv_y.^2;
    end
    m = zeros(size(G,1), size(G,2));
    for row = 1:size(G,1)
        for col = 1:size(G,2)
            if result(row,col) > T
                m(row,col) = 1;
            else
                m(row,col) = 0;
            end
        end
    end
end
