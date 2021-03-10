%DVB-T 2K Transmission
%The available bandwidth is 8 MHz
%2K is intended for mobile services
clearvars;
close all;

% -------------------------------------------------------------------------
%DVB-T Parameters
%--------------------------------------------------------------------------
Tu=224e-6;      % useful OFDM symbol period
T=Tu/2048;      % baseband elementary period
G=0;            % choice of 1/4, 1/8, 1/16, and 1/32
delta=G*Tu;     % guard band duration
Ts=delta+Tu;    % total OFDM symbol period
Kmax=2048;      % number of subcarriers
Kmin=0;
FS=4096;        % IFFT/FFT length
q=10;           % carrier period to elementary period ratio
fc=q*1/T;       % carrier frequency
Rs=4*fc;        % simulation period   % fc   2fc > 
t=0:1/Rs:Tu- 1/Rs;    % 226us   0 
tt=0:T/2:Tu - T/2;               
% -------------------------------------------------------------------------
%  FIR filters 
% -------------------------------------------------------------------------
p=1/Rs:1/Rs:T/2;                        
g=ones(length(p),1);

bt = 0.15;
span = 5;
sps = 4;

h1 = rcosdesign(bt, span, sps);
h2 = gaussdesign(bt, span, sps);
h1 = (h1(2:end)/sum(h1))*sum(g)';
h2 = (h2(2:end)/sum(h2))*sum(g)';

% Low pass filters consequently transmitter and reciever
[b,aa] = butter(12,1/20);
[B,AA] = butter(5,1/20);

% -------------------------------------------------------------------------
%  BER and SER score lists
% -------------------------------------------------------------------------
filter_BER = [];
filter_SER = [];
len_monte_carlo = 20;
BER_list = zeros(len_monte_carlo, 1);
SER_list = zeros(len_monte_carlo, 1);

% -------------------------------------------------------------------------
% Loop parameters
% -------------------------------------------------------------------------
fir_filter = g;
rng("shuffle")                 % Use default random number generato
snrrr = 0:2:20;
listt = [2 4 8 16 64];

% -------------------------------------------------------------------------
%                       Pilot index
% -------------------------------------------------------------------------
P = 256; % Number of Pilot carrier
pilot_step = Kmax/P;
pilot_value = 3+3j; 

pilot_carriers = 1:pilot_step:Kmax;
pilot_carriers = [pilot_carriers  Kmax];
data_carriers = 1:Kmax;
all_carriers = 1:Kmax;
data_carriers(pilot_carriers) = [];

% -------------------------------------------------------------------------

for kkkk = 1: length(listt)
    M = listt(kkkk);           % Size of signal constellation        
    k = log2(M);               % Number of bits per symbol
    n = (Kmax-P-1)*k ;         % Number of bits to process

for kkk = 1:length(snrrr)
snrr = snrrr(kkk)              % increase SNR gradualy

for kk = 1: len_monte_carlo
  
% Data generation 
dataIn = randi([0 1],n,1);                           % Generate vector of binary data
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary k-tuples, k = log2(M)
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers
data_freq = qammod(dataSymbolsIn,M);                 % Gray coding, phase offset = 1

% Pilot and Data Instertion
data = zeros(Kmax,1);
data(pilot_carriers) = pilot_value;
data(data_carriers) = data_freq;


% Zero padding for ifft
% -------------------------------------------------------------------------
A=length(data);   % 1706
info=zeros(FS,1); % 4096 ->  
info(1:(A/2)) = data(1:(A/2)).'; %Zero padding    % [data(:853)+ zeros + 
info((FS-((A/2)-1)):FS) =  data(((A/2)+1):A).';    %  data(854:end)  ] 

% Subcarriers generation (B) (OFDM)
%--------------------------------------------------------------------------
carriers=FS.*ifft(info,FS);

%       Add CP
% -------------------------------------------------------------------------
if G ~= 0  
    CP_len = length(carriers)*G;
    data_CP = carriers(fix(end-CP_len)+1: end);
    carriers = [data_CP ; carriers];
    t=0:1/Rs:Ts- 1/Rs;    % 226us   0 
end 

% ------------------------------------------------------------------------
% DAC
%--------------------------------------------------------------------------
% figure(1);
% subplot(211)
% plot(h1)
% subplot(212)
% plot(h2)
% --------------------------------------------------------------------------

L = length(carriers);
chips = [ carriers.';zeros((2*q)-1,L)];
p=1/Rs:1/Rs:T/2;
g=ones(length(p),1);
dummy=conv(fir_filter, chips(:));
u=[dummy; zeros(46,1)];
u = gain(carriers, u)*u;

uoft = filter(b,aa,u);

% Plots
% -------------------------------------------------------------------------
% figure(2);
% subplot(221);
% plot(t(80:480),real(u(80:480)));
% subplot(222);
% plot(t(80:480),imag(u(80:480)));
% 
% subplot(223);
% plot(t(80:480),real(uoft(80:480)));
% subplot(224);
% plot(t(80:480),imag(uoft(80:480)));

% Upconverter
% -------------------------------------------------------------------------
delay=64; %Reconstruction filter delay
s_tilde=(uoft(delay+(1:length(t))).').*exp(1i*2*pi*fc*t);
s=real(s_tilde);  

% Plots
%--------------------------------------------------------------------------
% figure(3);
% subplot(421);
% stem(tt(1:20),real(carriers(1:20)));
% title('Carrier signal real part 1-20 element')
% xlabel('time');
% subplot(422);
% stem(tt(1:20),imag(carriers(1:20)));
% title('Carrier signal imag part 1-20 element')
% xlabel('time');
% subplot(423);
% plot(t(1:400),real(u(1:400)));
% title('upsampled signal real part 1-400 element')
% xlabel('time');
% subplot(424);
% plot(t(1:400),imag(u(1:400)));
% xlabel('time');
% title('upsampled signal imag part 1-400 element')
% subplot(4,2,5);
% plot(t(1:400),real(uoft(1:400)));
% subplot(4,2,6);
% plot(t(1:400),imag(uoft(1:400)));
% title('low pass signal first 400 element')
% xlabel('time');
% subplot(4,2,[7 8]);
% plot(t(1:400),(s(1:400)));
% title('transmitted signal first 400 element')
% xlabel('time');

% ------------------------------------------------------------------------
% Frequency response plot from box filter
% ff=(Rs)*(1:(q*FS))/(q*FS);
% figure(4);
% subplot(311);
% plot(ff,abs(fftshift(fft(u,q*FS)))/FS);
% xlabel('f(MHz)')
% title('Frequency response with alising')
% subplot(312);
% plot(ff,abs(fftshift(fft(uoft,q*FS)))/FS);
% xlabel('f(MHz)')
% title('Frequency response after anti-alising filter')
% subplot(313);
% plot(ff,abs((fft(s,q*FS)))/FS);
% xlabel('f(MHz)')
% title('Up converted frequency response')

%--------------------------------------------------------------------------
%                                   Channel
%--------------------------------------------------------------------------
% z = [0.7 -0.5 0.4]; % channel example
z = [1, 0, 0.3+0.3j]; % channel example
s = filter(z,1,s);

s = awgn(s, snrr, 'measured');
%--------------------------------------------------------------------------
%                                RECEIVER 
%--------------------------------------------------------------------------

%   Downconversion
% -------------------------------------------------------------------------
r_tilde=exp(-1i*2*pi*fc*t).*s; %(F)

% -------------------------------------------------------------------------
% ff=(Rs)*(1:(q*FS))/(q*FS);

%--------------------------------------------------------------------------
%  Carrier suppression  low pass
%-------------------------------------------------------------------------
r_info = 2*filter(B, AA, r_tilde);

% Plots
% -------------------------------------------------------------------------
% figure(5);
% subplot(211);
% plot(t,real(r_info));
% axis([0 12e-7 -inf inf]);
% grid on;
% subplot(212);
% plot(t,imag(r_info));
% axis([0 12e-7 -inf inf]);
% grid on;
% -------------------------------------------------------------------------
% figure(6);
% f =(2/T)*(1:(FS))/(FS);
% subplot(211);
% plot(ff,abs(fft(r_tilde,q*FS))/FS);
% grid on;
% subplot(212);
% plot(ff,abs(fftshift(fft(r_info,q*FS)))/FS);
% grid on;


% -------------------------------------------------------------------------
% Sampling
% -------------------------------------------------------------------------

% Downsampling
r_data = decimate(r_info, 2*q);

% Removing cyclic prefix 
% -------------------------------------------------------------------------
if G ~= 0
    r_data = r_data( CP_len + 1 : end); 
end

% Plots
% -------------------------------------------------------------------------
% figure(7);
% subplot(321);
% plot((real(carriers(1:20))));
% grid on;
% title('rx real')
% subplot(322);
% plot((imag(carriers(1:20))));
% grid on;
% title('rx imag')
% subplot(323);
% plot((real(r_info(1:440))));
% title('tx upreal')
% grid on;
% subplot(324);
% plot((imag(r_info(1:400))));
% title('tx upimag')
% grid on;
% subplot(325);
% plot((real(r_data(1:20))));
% grid on;
% title('tx real')
% subplot(326);
% plot((imag(r_data(1:20))));
% grid on;
% title('tx imag')
ff=(Rs)*(1:(q*FS))/(q*FS);
% f =(2/T)*(1:(FS))/(FS);
% subplot(211);
% plot(ff,abs(fftshift(fft(r_tilde,q*FS)))/FS);
% grid on;
% subplot(212);
% plot(ff,abs(fftshift(fft(r_info,q*FS)))/FS);
% grid on;

% ------------------------------------------------------------------------
% FFT
%-------------------------------------------------------------------------
info_2N=(1/FS).*fft(r_data,FS); % (I)
info_h=[info_2N(1:A/2) info_2N((FS-((A/2)-1)):FS)];

%-------------------------------------------------------------------------
%               Channel Estimation
%-------------------------------------------------------------------------
pilots = info_h(pilot_carriers);
Hest_at_pilots = pilots / pilot_value ;

Hest_abs = interp1(pilot_carriers,abs(Hest_at_pilots),all_carriers,'PCHIP');
Hest_phase  = interp1(pilot_carriers,angle(Hest_at_pilots),all_carriers,'PCHIP');
Hest = Hest_abs .* exp(1j*Hest_phase);

info_h = info_h ./ Hest; % eliminate channel effect

% Plots
% -------------------------------------------------------------------------
% figure(8);
% f=(2/T)*(1:(FS))/(FS);
% subplot(221)
% plot(f,abs(fft(carriers,FS))/FS);
% subplot(222)
% stairs(real(data))
% subplot(223)
% plot(f,abs(fft(r_data,FS))/FS);
% subplot(224)
% stairs(real(info_h))

% selecting data carriers
data_freq_out = info_h(data_carriers);

%-------------------------------------------------------------------------
% Construction data 
%-------------------------------------------------------------------------
dataSymbolsOut = qamdemod(data_freq_out,M)';
% ------------------------------------------------------------------------
% Construcation binary data
dataOutMatrix = de2bi(dataSymbolsOut);                 % Convert to integers
dataOut = reshape(dataOutMatrix,[],1);   % Reshape data into binary k-tuples, k = log2(M)

% Plots
% ------------------------------------------------------------------------
% figure(9)
% subplot(211)
% plot(real(data(1:100)));
% subplot(212)
% plot(real(info_h(1:100)));

% ------------------------------------------------------------------------
% figure(10)
% subplot(211)
% stairs(dataIn(1:40),"b")
% axis([1 40 -0.5 1.5])
% title("DataIn First 40 bit")
% subplot(212)
% stairs(dataOut(1:40), "b")
% title("DataOut First 40 bit")
% axis([1 40 -0.5 1.5])


% Collection information
BER = biterr(dataIn, dataOut);
SER = symerr(dataSymbolsIn, dataSymbolsOut);

BER_list(kk) = BER;
SER_list(kk) = SER;
end
BER_mean(kkk) = mean(BER_list);
SER_mean(kkk) = mean(SER_list);
end
filter_BER = [filter_BER; BER_mean];
filter_SER = [filter_SER; SER_mean];
end
filter_SER = filter_SER/length(dataSymbolsOut);
filter_BER = filter_BER/length(dataOut);



% BER and SER plots
figure()
labels = ["2","4","8","16", "64"];
subplot(121)
plot(snrrr, filter_BER);

legend(labels);
xlabel("SNR (dB)")
ylabel("Bit Error Rate")
title("Window filter 2k Mod, awgn SNR versus SER rate")
set(gca, 'YScale', 'log')

subplot(122)
plot(snrrr, filter_SER);

legend(labels);
xlabel("SNR (dB)")
ylabel("Sym Error Rate")
title("Window filter 2k Mod, awgn SNR versus BER rate")
set(gca, 'YScale', 'log')

% Auto gain calculator
function g = gain(s1, s2)
g = min(real(s1))/min(real(s2));
g = abs(g);
end
