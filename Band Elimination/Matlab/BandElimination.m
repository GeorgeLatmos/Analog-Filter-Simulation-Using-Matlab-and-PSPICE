%Script to design Band Elimination Filter
%Giorgos Latmos
%AEM:8683
a3 = 8;
a4 = 3;
%Define Filter specifications
f0 = 2400;
f1 = 1725+25*a4;
f2 = f0^2/f1;
D = 1/2.5*((f0^2-f1^2)/f1);
f3 = (-D+sqrt(D^2+4*f0^2))/2;
f4 = f0^2/f3;
amin = 26+a3*5/9;
amax = 0.5+a4/18;
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
Omega_s = (w2-w1)/(w4-w3);

%Design prototype Low Pass Filter

%Calculate filter class
n = acosh(sqrt((10^(amin/10)-1)/(10^(amax/10)-1)))/acosh(Omega_s);
n = ceil(n);

%Calculate e parameter
e = 1/(sqrt(10^(amin/10)-1));

%Calculate a parameter
a = 1/n*asinh(1/e);

%Calculate Normalized 3 DB frequency
Omega_hp = 1/(cosh(1/n*acosh(1/e)));

%Calculate Butterworth Angles and poles
y = zeros(n,1);
S = zeros(n,1);
for k=1:n
    phi = radtodeg((-(2*k-1)*pi)/(2*n));
    y(k) = 90-phi;
    S(k) = sinh(a)*cosd(y(k))+i*cosh(a)*sind(y(k));
end

%Calculate Inverse Chebyshev Poles
Omega_IC_norm = zeros(n,1);
for k=1:n
   Omega_norm = abs(S(k));
   Omega_IC_norm(k) = 1/Omega_norm*Omega_s;
end

%Calcuate TF Zeros
Omega_z = zeros(n/2,1);
index = 1;
for k=1:2:n
   Omega_z(index) = 1/cos(k*pi/8)*Omega_s;
   index = index+1;
end

%Invert Inverse Chebyshev Poles
Omega_IC_inv = zeros(n,1);
for k=1:n
   Omega_IC_inv(k) = 1/Omega_IC_norm(k);
end

%Invert Zeros
Omega_zeros_inv = zeros(n/2,1);
for k=1:n/2
   Omega_zeros_inv(k) = 1/Omega_z(k);
end

%Calculate HIGH PASS poles
S_HP = zeros(n,1);
for k=1:2:n
    Q = abs(abs(S(k))/(2*real(S(k))));
    real_part = -Omega_IC_inv(k)/(2*Q);
    imag_part = sqrt(Omega_IC_inv(k)^2-real_part^2);
    S_HP(k) = real_part+i*imag_part;
    S_HP(k+1) = real_part-i*imag_part;
end

%Transform HIGH PASS poles, using Geffe Algorithm
Omega_final = zeros(n,1);
Q_final = zeros(n,1);
index = 1;
for k=1:2:n
    Sigma = abs(real(S_HP(k)));
    Omega = abs(imag(S_HP(k)));
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

%Transform zeros
Omega_z_final = zeros(n,1);
index = 1;
for k=1:length(Omega_zeros_inv)
   K = 2+(Omega_zeros_inv(k)^2)/(qc^2);
   x = (K+sqrt(K^2-4))/2;
   Omega_z_final(index) = sqrt(x)*w0;
   Omega_z_final(index+1) = 1/sqrt(x)*w0;
   index = index+2;
end

message = ['Total Number of Units: ', num2str(n)];
disp(message);
for i=1:n
   message = ['Unit ', num2str(i),': w0 = ',num2str(Omega_final(i)/(2*pi)),', ','wz = ',num2str(Omega_z_final(i)/(2*pi)), ', Q = ', num2str(Q_final(i))];
   disp(message);
end

%Calculate Transfer Functions
T1 = tf([1 0 Omega_z_final(1)^2],[1 Omega_final(1)/Q_final(1) Omega_final(1)^2]);
T2 = tf([1 0 Omega_z_final(2)^2],[1 Omega_final(2)/Q_final(2) Omega_final(2)^2]);
T3 = tf([1 0 Omega_z_final(3)^2],[1 Omega_final(3)/Q_final(3) Omega_final(3)^2]);
T4 = tf([1 0 Omega_z_final(4)^2],[1 Omega_final(4)/Q_final(4) Omega_final(4)^2]);
temp1 = series(T1,T2);
temp2 = series(T3,T4);
T = series(temp1,temp2);
clear temp1 temp2
clear a bw e Omega_hp Omega_IC_inv Omega_IC_norm Omega_norm Omega_p Omega_s Omega_z Omega_zeros_inv qc S S_HP
clear x y i message phi real_part Sigma Omega C D E G Q K W index imag_part k Omega_inv_norm

%Filter Design

%Design 1st LPN filter --> BoctorLPN
BoctorLPN(Omega_final(2), Omega_z_final(2), Q_final(2),'Unit2');
%Design 2nd LPN filter --> BoctorLPN
BoctorLPN(Omega_final(4), Omega_z_final(4), Q_final(4),'Unit4');
%Design 1st HPN filter --> HPN (simple)
simpleHPN(Omega_final(1), Omega_z_final(1), Q_final(1),'Unit1');
%Design 2nd HPN filter --> HPN (simple)
simpleHPN(Omega_final(3), Omega_z_final(3), Q_final(3), 'Unit3');

clear Omega_final Omega_z_final Q_final

ltiview(T);

w1_sim = w0-(w0-w3)/2;
w2_sim = w0+(w0+w3)/3;
w3_sim = 0.4*w1;
w4_sim = 2.5*w2;
w5_sim = 3*w2;

f1_sim = w1_sim/(2*pi);
f2_sim = w2_sim/(2*pi);
f3_sim = w3_sim/(2*pi);
f4_sim = w4_sim/(2*pi);
f5_sim = w5_sim/(2*pi);

[u1, t] = gensig('sin',1/f1_sim,7*10^-3,10^-4);
u1 = 0.5*u1;
[u2, t] = gensig('sin',1/f2_sim,7*10^-3,10^-4);
u2 = 0.8*u2;
[u3, t] = gensig('sin',1/f3_sim,7*10^-3,10^-4);
u3 = 0.8*u3;
[u4, t] = gensig('sin',1/f4_sim,7*10^-3,10^-4);
u4 = 0.6*u4;
[u5, t] = gensig('sin',1/f5_sim,7*10^-3,10^-4);
u5 = 1.2*u5;

u=u1+u2+u3+u4+u5;
%lsim(T,u,t)

y = lsim(T,u,t);
y = fft(y);
y_new = zeros(ceil(length(y)/2),1);
for k=1:length(y)/2
    y_new(k) = y(k);
end
plot(abs(y_new))
title('FFT of output signal');
xlabel('Frequency') % x-axis label
ylabel('Amplitude') % y-axis label






















