clear all 
close all



%*********************************%
% ADAPTIVE FILTRATION APPLICATION %
%*********************************%



% STEPS: a) ... j)

% STEP:a)____________________________________________________
% load audio record
% record parameters: 44,1kHz, 16b, mono
[s,Fe] = audioread('signal_init.wav');

ns = length(s);

% build time array and plot
n = 0:1/(ns-1):1;

plot(n,s), grid									
	title('Semnalul vocal neprelucrat') 
	xlabel('Timp')
	ylabel('s(n)')

    
% STEP:b)____________________________________________________    
% normalize amplitude
% scale signal in range [-1,1]
scale = max(s(:)) - min(s(:));
s_01 = (s - min(s(:))) / scale;
s1 = 2*s_01 - 1;

figure
plot(n,s1), grid											
	title('Semnalul vocal scalat in [-1, 1]') 
	xlabel('Timp (s)')
	ylabel('s1(n)')
    
    
% STEP:c)____________________________________________________
% decimate signal

% fs_in = 44.1 kHz
% fs_out = 16 kHz
% decimate with rational factor
% R(resampling factor) = fs_out / fs_in = 16 / 44.1
% but must be type int => R = 160 / 441
interpolation = 160; %(interpolation factor)
decimation = 441; %(decimation factor)
% First do interpolation, second do decimation

s2 = resample(s1, interpolation, decimation);

% update signal params
Fe = 16000;
ns = length(s2);
n = 0:1/(ns-1):1;

figure
plot(n,s2), grid
	title('Semnalul s2(dupa resample)') 
	xlabel('Timp (s)')
	ylabel('s2(n)')
    
S2 = abs(fft(s2));
figure
plot((0:(ns-1))/ns*Fe, S2), grid                                
	title('Modul spetru discret semnalul S2') 
	xlabel('f_in / f_out')
	ylabel('|S2(k)|')
    

% STEP:d)____________________________________________________
% filter s2 in a FTJ, fc = 7kHz
% 7kHz is perfect for intelligible and naturalness
% (50Hz ... 7kHz is large band voice signal)

% create FTJ IIR
% 7KHz is double voice max fignal freq(2 * 3.4k = 6.8, round upper to 7k)
% freqz(b,a,no_points_fot_dft, fe)
[b, a] = butter(12,7000/(Fe/2));
[H, w] = freqz(b,a,[],Fe);
figure
plot(w,abs(H)),grid											
	title('Filtru digital de 7KHz') 
	xlabel('F(Hz)')
	ylabel('|H(z)|')
    
% filter signal s2
s = filter(b,a,s2);
figure
plot(n,s)											
	title('Semnalul s trecut prin filtrul FTJ') 
	xlabel('Time(s)')
	ylabel('s(n)')
    
S = abs(fft(s)); 
figure
plot((0:(ns-1))/ns*Fe, S), grid
 	title('Modul spectru discret, semnalul S') 
 	xlabel('Frecventa (Hz)')
 	ylabel('|S(k)|')
    
% save signal 's' on disk
audiowrite('C:\signal.wav', s, Fe);



% STEP:e)____________________________________________________
% create Distortion filter
% signal 's' is now input for Distortion filter
% Distortion filter order(as of any IIR filter) = min(N,M) - 1
% so in this case 16

% b and a are:
b = [0.1662, -0.0943, 0.2892, -0.1227, 0.2348, 0.0180, 0.0415, 0.1388, -0.0616, 0.1290, -0.0434, 0.0420, -0.0010, -0.0009, 0.0032, -0.0015, 0.0056];
a = [1.0000, -0.7548, 3.4400, -1.6385, 4.8436, -0.8156, 3.2813, 1.2582, 0.6571, 2.1922, -0.4792, 1.4546, -0.2905, 0.4693, -0.0208, 0.0614, 0.0120];
[Hd, wd] = freqz(b,a,[], Fe);
figure
plot(wd,abs(Hd))											
	title('Filtru de distorsionare') 
	xlabel('Frecventa(Hz)')
	ylabel('|Hd(z)|')

% pass signal 's' to the Distortion filter
x = filter(b,a,s);



% STEP:f)____________________________________________________
%inverstigate s and x in time and frequency
%s was already ploted above(at step d))

figure
plot(n,x)											
	title('Semnalul x (distorsionat)') 
	xlabel('Time(s)')
	ylabel('x(n)')

X = abs(fft(x)); 
figure
plot((0:(ns-1))/ns*Fe, X), grid
 	title('Modul spectru discret, semnalul X') 
 	xlabel('Frecventa (Hz)')
 	ylabel('|X(k)|')
    
% save signal 'x' on disk
audiowrite('C:\signal_distorted.wav', x, Fe);



% STEP:g)____________________________________________________
% implement adaptive filter using NLMS algorithm

%generate noise
s_aux = randn(size(n));
    
% pass s_aux to the distortion filter
x_aux = filter(b,a,s_aux);

% x_aux is distortion filter output
% then F.A. input in x_aux
% desired signal is s_aux
% and F.A. output is y_aux
x = x_aux;
d = s_aux;

% calculate energy algorithm
% choose N in range 10...30ms
N = 160; % window = 160 / 16k = 10ms or N = 480 for 30ms
M = fix(ns/N);	%no. of analysis windows during the signal s;

E = zeros(1,M);	%initialize
i = 1;

while i <= M
    for k = 1:N
        E(i) = E(i) + (x(k+(i-1)*N))^2;
    end
    i = i+1;
end
% end of energy calculus algorithm

% Plot pseudo-energy of input signal in F.A.
n1 = 0:1/(M-1):1;
figure
plot(n1,E)
    title('Energia semnalului de la intrarea F.A.') 
    xlabel('Timp (s)')
   	ylabel('Energie')
    
% assing parameters
mu = 0.1; sigma = sum(E); alpha = 0.01; L = 50;

% create F.A. with all above parameters and pass them to nlms.m function
b = zeros(1,L+1); px=0;
[y_aux,b,px] = nlms(x,d,b,mu,sigma,alpha,px);

SE = (d-y_aux).^2;	%square error
figure
plot([0:(length(n)-1)],SE), grid								
	title ('Eroarea patratica')
	xlabel('n')
	ylabel('SE')

    
    
% STEP:h)____________________________________________________
% investigate s_aux, x_aux and y_aux

% s_aux:
figure
plot(n,s_aux)											
	title('Zgomot gaussian standard(semnal cobai)') 
	xlabel('Time(s)')
	ylabel('s_aux(n)')

S_AUX = abs(fft(s_aux)); 
figure
plot((0:(ns-1))/ns*Fe, S_AUX), grid
 	title('Modul spectru discret, zgomot gaussian') 
 	xlabel('Frecventa (Hz)')
 	ylabel('|S_AUX(k)|')
    
    
% x_aux:
figure
plot(n,x_aux)											
	title('Semnalul x_aux(semnal cobai iesire filtru distorsionare)') 
	xlabel('Time(s)')
	ylabel('x_aux(n)')

X_AUX = abs(fft(x_aux)); 
figure
plot((0:(ns-1))/ns*Fe, X_AUX), grid
 	title('Modul spectru discret, semnalul X_AUX(zg gaussian)') 
 	xlabel('Frecventa (Hz)')
 	ylabel('|X_AUX(k)|')
    
    
% y_aux:
figure
plot(n,y_aux)											
	title('Semnalul y_aux: iesire F.A.') 
	xlabel('Time(s)')
	ylabel('y_aux(n)')

Y_AUX = abs(fft(y_aux)); 
figure
plot((0:(ns-1))/ns*Fe, Y_AUX), grid
 	title('Modul spectru discret, semnalul Y_AUX') 
 	xlabel('Frecventa(Hz)')
 	ylabel('|Y_AUX(k)|')

% compute product in frequency
[H,w] = freqz(b,1,[],Fe);
product = H .* Hd; % compute product as an scalar, by using mean value
product = sum(product) / 512; %(because H and Hd have 512 points)
% see product result in Workspace

% or compute product as a matrix
product2 = Hd.*H;
figure
plot(wd,abs(product2))
	title('Produs') 
	xlabel('Frecventa(Hz)')
	ylabel('|Hd(z) * H(z)|')



% STEP:i)____________________________________________________
% x -> |F.A.| -> y
    
% b coeffs are those returned from nlms
% a coffes(1...N) are 0 (Adaptive Filters are FIR)
% a0 = 1 because of "ecuatia cu diferente finite"
[H,w] = freqz(b,1,[],Fe);
figure
plot(wd,abs(Hd),':b','linewidth',2), grid				
hold on
	title ('Filtru adaptiv: modelare inversa')
	xlabel('Frecventa (Hz)')
	ylabel('Modulul raspunsului in frecventa')
plot(w,abs(1./H),'r','linewidth',1.5)
	legend('Sistem necunoscut','Filtru adaptiv')
	legend('Location','NorthEast')
hold off

[x,Fe] = audioread('signal_distorted.wav');
%sound(x, Fe);
%pause(8)
y_bossul_bossilor = filter(b,1,x);
%sound(y_bossul_bossilor, Fe)

audiowrite('C:\signal_recuperat.wav', x, Fe);



% STEP:j)____________________________________________________
% investigate s,x,y
% s and x were already ploted at some above steps
% but we plot it again here to be all grouped

% s: read it from disc, we want s before beeing distored
[s,Fe] = audioread('signal.wav');

figure
plot(n,s)											
	title('Semnalul util nedistorsionat') 
	xlabel('Time(s)')
	ylabel('s(n)')

S = abs(fft(s)); 
figure
plot((0:(ns-1))/ns*Fe, S), grid
 	title('Modul spectru discret semnalul util nedistorsionat') 
 	xlabel('Frecventa (Hz)')
 	ylabel('|S(k)|')
    
    
% x: read it from disc, we want x = s distored
[x,Fe] = audioread('signal_distorted.wav');
figure
plot(n,x)											
	title('Semnalul util distorsionat') 
	xlabel('Time(s)')
	ylabel('x(n)')

X = abs(fft(x)); 
figure
plot((0:(ns-1))/ns*Fe, X), grid
 	title('Modul spectru discret semnalul util distorsionat') 
 	xlabel('Frecventa (Hz)')
 	ylabel('|X(k)|')
    
    
% y:
figure
plot(n,y_bossul_bossilor)											
	title('Semnalul util recuperat prin F.A.') 
	xlabel('Time(s)')
	ylabel('y(n)')

Y_bossul_bossilor = abs(fft(y_bossul_bossilor)); 
figure
plot((0:(ns-1))/ns*Fe, Y_bossul_bossilor), grid
 	title('Modul spectru discret semnalul util recuperat prin F.A.') 
 	xlabel('Frecventa (Hz)')
 	ylabel('|Y(k)|')
    
