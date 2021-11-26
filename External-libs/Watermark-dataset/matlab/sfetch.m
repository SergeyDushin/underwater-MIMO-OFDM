function [y, fs] = sfetch(signal, channel, packetNumber, SNR);

% Synopsis: [y, fs] = sfetch(signal, channel, packetNumber, SNR);
%
% Serial fetch: retrieve packets for temporal data.
%
% Input :  signal      : string variable denoting the signal (e.g. 'my_signal')
%          channel     : string variable denoting the acoustic channel (e.g. 'NOF1')
%          packetNumber: an integer value in the range [1, Nmax], where Nmax is the total number
%                        of filtered packets (bk.nPackets)
%          SNR         : numerical value -> Eb/N0 value in dB 
%                        Omit this argument to retrieve the packet without additive noise
%        
% Output:  y  : time series of the distorted signal 
%          fs : sampling frequency
%
%          NB. When an SNR is specified, the returned time series is ~10 seconds longer than 
%              the actual packet, which arrives at a random offset between ~ 4 and 6 s. 
%
% Watermark version 1.0
% Forsvarets Forskningsinstitutt, 
% 03.11.2016



%% Set path to Watermark base directory
P=mfilename('fullpath');
tmp = findstr(lower(P), 'matlab');
waterMarkPath=P(1:tmp(end)-1);

% load bookkeeping information
filename = fullfile(waterMarkPath, 'output', channel, signal, 'bookkeeping');
try
    load(filename);
catch
    disp([filename ' not found']);
    return
end

% check validity of packet number
if packetNumber <= 0 | packetNumber > bk.nPackets
    fprintf(['Invalid packet number. Use a value between ' sprintf('%d', 1) ' and ' sprintf('%d', bk.nPackets) '.\r']);
    y=[]; fs=[];
    return
end

% Find location of the desired packet
i_sounding = ceil(packetNumber/bk.nPacketsPerSounding);
i_packet = packetNumber - (i_sounding-1)*bk.nPacketsPerSounding;

fullfile(waterMarkPath, 'output', channel, signal, [channel '_' sprintf('%03d', i_sounding)]);
filename = fullfile(waterMarkPath, 'output', channel, signal, [channel '_' sprintf('%03d', i_sounding) '.wav']);
%fprintf(['Fetching packet ' sprintf('%d', packetNumber) '/' sprintf('%d', bk.nPackets) '\r']);
try
    [data, fs]=audioread(filename, [bk.packetIndex(i_packet,1), bk.packetIndex(i_packet,2)]);
catch
    disp(['Error reading data from ' filename]);
    return
end

%% A Doppler shift V0 has been removed from the data prior to the channel estimation.
% This has been done by resampling the raw acoustic data with a resampling factor 1-V0/c. 
% This Doppler shift will now re-instated. The used sign convention is such that a positive 
% velocity corresponds to a positive range rate, i.e. time dilation of the signal. 
% Small values may be entirely due to clock frequency offsets of transmitter and receiver
% hardware, which cause an apparent Doppler shift.  
% 
% NB. Resampling factors close to one may lead to memory problems, but the NCS1,NOF1,BCH1,
% KAU1, and KAU2 channels should be fine on contemporary work stations.


c = 1500;  % nominal sound speed in m/s
V0 = bk.V(i_sounding);  % Doppler velocity in m/s
resamplingFactor = 1/(1-V0/c);
[N1,D1] = rat(resamplingFactor);
try
    data = resample(data, N1, D1);
catch
    % Resampling factor is too close to 1. 
    disp('Warning: Doppler shift is not re-instated.');
end

if nargin == 4

    Eb = sum(data.^2)/bk.nBits; % energy per bit
 
    % Detection and synchronization are essential tasks for a communication receiver. 
    % In order to avoid a predictable signal start, the packet is zero padded before noise is added.
    % The packet has a random offset (starting somewhere between 4 and 6 s in the output waveform.
    rng('default');
    rng(packetNumber)
    i_start = round((4+2*rand)*fs);
    y = zeros(length(data) + round(10*fs), 1);
    y(i_start:i_start+length(data)-1) = data;
    
    % Create random white noise vector with corresponding power spectral density N0 = 2;
    awgn = randn(length(y), 1);
    N0 = 2; 
    
    % Scale noise and add to signal
    scalingfactor = 10*log10(Eb/N0) - SNR;
    awgn = awgn*10^(scalingfactor/20);
    y = y + awgn; % noisy signal at the specified SNR

else
    
    % Return packet as is
    y=data;

end
   

return

