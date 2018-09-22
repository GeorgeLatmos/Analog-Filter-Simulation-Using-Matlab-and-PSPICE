%Script to design High-Pass Filter
%Giorgos Latmos
%AEM:8683
m = 1;

%Define Filter specifications
fp = (4+m)*1000;
fs = fp/2.6;
amin = 24+3*6/9;
amax = 0.5+8/36;

%Convert Hz to rad/sec
ws = 2*pi*fs;
wp = 2*pi*fp;
Omega_s = wp/ws;

%Calculate filter class
n = log10((10^(amin/10)-1)/(10^(amax/10)-1))/(2*log10(Omega_s));
n = ceil(n);

%Calculate Normalized 3 DB frequency
Omega_hp = 1/((10^(amax/10)-1)^(1/10));

%High Pass Filter frequency
w0 = wp/Omega_hp;

%Normalized Butterworth poles for n=5 (Read from table)
P = zeros(n,1);
Q = zeros(n,1);
P(1) = -1;
Q(1) = 0.5;
P(2) = -0.809+j*0.5877;
Q(2) = 0.618;
P(3) = -0.809-j*0.5877;
Q(3) = 0.618;
P(4) = -0.309+j*0.951;
Q(4) = 1.618;
P(5) = -0.309-j*0.951;
Q(5) = 1.618;


%Calculate Unit-1 parameters
PassiveFilter(w0, 'Unit1');
%Calculate Unit-2 parameters
StrategyII(w0,Q(2),'Unit2');
%Calculate Unit-3 parameters
StrategyII(w0,Q(4),'Unit3');

temp = series(Unit1.TF,Unit2.TF);
T = series(temp,Unit3.TF);

ltiview(Unit1.TF, Unit2.TF, Unit3.TF, T)
%clear m n Omega_hp Omega_s P temp w0

w1_sim = 0.4*ws;
w2_sim = 0.9*ws;
w3_sim = 1.4*wp;
w4_sim = 2.4*wp;
w5_sim = 4.5*wp;

f1_sim = w1_sim/(2*pi);
f2_sim = w2_sim/(2*pi);
f3_sim = w3_sim/(2*pi);
f4_sim = w4_sim/(2*pi);
f5_sim = w5_sim/(2*pi);

clear w1_sim w2_sim w3_sim w4_sim w5_sim

%ltiview(T)

% %Simulation
[u1, t] = gensig('sin',1/f1_sim,6*10^-3,10^-4);
[u2, t] = gensig('sin',1/f2_sim,6*10^-3,10^-4);
u2 = u2*0.5;
[u3, t] = gensig('sin',1/f3_sim,6*10^-3,10^-4);
[u4, t] = gensig('sin',1/f4_sim,6*10^-3,10^-4);
u4 = u4*0.7;
[u5, t] = gensig('sin',1/f4_sim,6*10^-3,10^-4);
u5 = u5*0.5;

u=u1+u2+u3+u4+u5;

y=lsim(T,u,t);
y = fft(u);
y_new = zeros(ceil(length(y)/2),1);
for i=1:length(y)/2
   y_new(i) = y(i);
end
plot(abs(y_new));
title('FFT of input signal');
xlabel('Frequency') % x-axis label
ylabel('Amplitude') % y-axis label

clear i m n Omega_hp Omega_s P Q t temp






