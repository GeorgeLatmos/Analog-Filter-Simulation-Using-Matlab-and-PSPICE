%Script to design Band-Pass Filter
%Giorgos Latmos
%AEM:8683
a3 = 8;
a4 = 3;
%Define Filter specifications
f0 = 900;
f1 = 650+25*a3;
f2 = (f0^2)/f1;
D = 2.2*(f0^2-f1^2)/f1;
f3 = (-D+sqrt(D^2+4*f0^2))/2;
f4 = f0^2/f3;
amin = 28+a4*5/9;
amax = 0.5+a3/36;
clear a3 a4 D

%Convert Hz to rad/sec
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

%Calculate center frequency
w0 = sqrt(w1*w2);
%Calculate Pass Band
bw = w2-w1;

%Normalized Frequencies for prototype Low Pass Filter
Omega_p = 1;
Omega_s = (w4-w3)/(w2-w1);

%Design prototype Low Pass Filter

%Calculate filter class
n = acosh(sqrt((10^(amin/10)-1)/(10^(amax/10)-1)))/acosh(Omega_s);
n = ceil(n);

%Calculate e parameter
e = sqrt(10^(amax/10)-1);

%Calculate Normalized 3 DB frequency
Omega_hp = cosh(1/n*acosh(1/e));

%Calculate a parameter
a = 1/n*asinh(1/e);

%Calculate Butterworth Angles and poles
y = zeros(n,1);
S = zeros(n,1);
for k=1:n
    phi = radtodeg((-(2*k-1)*pi)/(2*n));
    y(k) = 90-phi;
    S(k) = sinh(a)*cosd(y(k))+i*cosh(a)*sind(y(k));
end

%Transform poles, using Geffe Algorithm
Omega_final = zeros(n,1);
Q_final = zeros(n,1);
index = 1;
for k=1:2:n
    Sigma = abs(real(S(k)));
    Omega = abs(imag(S(k)));
    C = Sigma^2+Omega^2;
    qc = w0/bw;
    D = 2*Sigma/qc;
    E = 4+C/(qc^2);
    G = sqrt(E^2-4*D^2);
    Q = 1/D*sqrt(1/2*(E+G));
    K = Sigma*Q/qc;
    W = K+sqrt(K^2-1);
    Omega_final(index) = W*w0;
    Omega_final(index+1) = 1/W*w0;
    Q_final(index) = Q;
    Q_final(index+1) = Q;
    index = index+2;
end

clear Sigma Omega C qc D E G Q K W index

Strategy_Enhancement(Omega_final(1),Q_final(1),'Unit1');
Strategy_Enhancement(Omega_final(2),Q_final(2),'Unit2');
Strategy_Enhancement(Omega_final(3),Q_final(3),'Unit3');
Strategy_Enhancement(Omega_final(4),Q_final(4),'Unit4');

temp1 = series(Unit1.TF,Unit2.TF);
temp2 = series(Unit3.TF,Unit4.TF);
T = series(temp1,temp2);
clear temp1 temp2
clear k Omega_p phi

message = '....Band-Pass Chebyshev Filter....';
disp(message);
message = ['Total Number of Units: ', num2str(n)];
disp(message);
message = ['Unit 1: w0 = ',num2str(Unit1.w0),', Q = ', num2str(Unit1.Q)];
disp(message);
message = ['Unit 2: w0 = ',num2str(Unit2.w0),', Q = ', num2str(Unit2.Q)];
disp(message);
message = ['Unit 3: w0 = ',num2str(Unit3.w0),', Q = ', num2str(Unit3.Q)];
disp(message);
message = ['Unit 4: w0 = ',num2str(Unit4.w0),', Q = ', num2str(Unit4.Q)];
disp(message);
message = 'Circuit Parameters can be found in existing structs (Unit1, Unit2, Unit3, Unit4)';
disp(message);
clear message

s=j*w0;
a1 = 1.7782/norm(evalfr(T,s)); % Regulate Gain at 5 DB

T = a1*T;

ltiview(T);

clear a e k n Omega_final Omega_hp Omega_s Q_final s S y bw

w1_sim = w0-(w0-w1)/2;
w2_sim = w0+(w0+w1)/2;
w3_sim = 0.5*w3;
w4_sim = 2.4*w4;
w5_sim = 3.5*w4;

f1_sim = w1_sim/(2*pi);
f2_sim = w2_sim/(2*pi);
f3_sim = w3_sim/(2*pi);
f4_sim = w4_sim/(2*pi);
f5_sim = w5_sim/(2*pi);

[u1, t] = gensig('sin',1/f1_sim,15*10^-3,10^-4);
[u2, t] = gensig('sin',1/f2_sim,15*10^-3,10^-4);
u2 = u2*0.6;
[u3, t] = gensig('sin',1/f3_sim,15*10^-3,10^-4);
[u4, t] = gensig('sin',1/f4_sim,15*10^-3,10^-4);
u4 = u4*0.8;
[u5, t] = gensig('sin',1/f4_sim,15*10^-3,10^-4);
u5 = u5*0.4;

u=u1+u2+u3+u4+u5;

y = lsim(T,u,t);
y = fft(y);
y_new = zeros((length(y)+1)/2,1);
for k=1:(length(y)+1)/2
    y_new(k) = y(k);
end
plot(abs(y_new))
title('FFT of output signal');
xlabel('Frequency') % x-axis label
ylabel('Amplitude') % y-axis label








