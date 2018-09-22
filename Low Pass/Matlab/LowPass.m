%Script to design Low-Pass Filter
%Giorgos Latmos
%AEM:8683
m = 2;

%Define Filter specifications
fp = (3+m)*1000;
fs = 1.75*fp;
amin = 24+(8-5)*3/4;
amax = 0.75+(3-5)*0.25/4;

%Convert Hz to rad/sec
ws = 2*pi*fs;
wp = 2*pi*fp;
Omega_s = ws/wp;

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

%Calculate Omega and Q for each pole
Omega_final = zeros(n,1);
Q_final = zeros(n,1);
for k=1:n
   Q_final(k) = abs(abs(S(k))/(2*real(S(k))));
   Omega_final(k) = abs(S(k))*wp;
end


%1st order passive filter
PassiveFilter(Omega_final(3), abs(real(S(3))),'Unit1');
%Sallen-Key/Strategy II
StrategyII(Omega_final(1),Q_final(1),'Unit2')
%Sallen-Key/Strategy II
StrategyII(Omega_final(2),Q_final(2),'Unit3');

clear k m phi 

temp = series(Unit1.TF, Unit2.TF);
T = series(temp, Unit3.TF);
clear temp

%ltiview(Unit1.TF,Unit2.TF,Unit3.TF,T)

%Simulation
%Create square waves / T=0.5ms / ?=0.2ms
tau = 0.2*10^-3;
Tf = 0.2*10^-3;
Ts = 10^-5;
numOfWaves = 4;
[ u,t ] = generateSquareWave(tau, Tf, Ts, numOfWaves);

%Plot Results
y = lsim(T,u,t);
y = fft(u);
y_new = zeros(length(y)/2,1);
for k=1:length(y)/2
    y_new(k) = y(k);
end
plot(abs(y_new))
title('FFT of input signal');
xlabel('Frequency') % x-axis label
ylabel('Amplitude') % y-axis label

clear numOfWaves tau Tf Ts u t f a e k n Omega_hp Omega_s Q_final S wp ws y y_new Omega_final 






