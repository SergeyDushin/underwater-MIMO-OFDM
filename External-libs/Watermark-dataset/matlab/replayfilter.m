function y = replayfilter(x, fs_x, h, fs_t, fs_tau, fc);

% Synopsis: y = replayfilter(x, fs_x, h, fs_t, fs_tau, fc);
%
% Input :  x      : input signal
%          fs_x   : sampling rate of x  (Hz)
%          h      : h(t,tau) -> time varying impulse response
%          fs_t   : sampling rate of h in time (t)
%          fs_tau : sampling rate of h in delay (tau)
%          fc     : center frequency of the channel (i.e., of the channel probe signal)
%        
% Output:  y  : input signal convolved with h(t,tau)
%
% Watermark version 1.0
% Forsvarets Forskningsinstitutt, 
% 03.11.2016


% bring signal to baseband
t=[0:length(x)-1]'/fs_x;
x=x.*exp(-2*pi*fc*1i*t);

% resample signal to sampling frequency of channel estimate
[N,D]=rat(fs_tau/fs_x);
x=resample(x,N,D);

y=0*x; % initiate output signal
h=flipud(h.'); % rotate and flip
[K dummy]=size(h);
x=[zeros(K-1,1); x]; % zero padding
Ly=length(y);

% Direct-replay channel simulation
for k=0:Ly-1
    f=k/K-floor(k/K); 
    n=floor(k/K)+1; % impulse response counter
    ir=(1-f)*h(:,n)+f*h(:,n+1); % linear interpolation
    y(k+1)=x(k+1:k+K).'*ir; % filtering
end

% resample waveform to sampling frequency of input signal
[N,D]=rat(fs_x/fs_tau);
y=resample(y,N,D);

% shift to passband
t=[0:length(y)-1]'/fs_x;
y=real(y.*exp(2*pi*fc*1i*t));

return
