function [y, fs] = pfetch(signal, channel, packetNumber);

% Synopsis: [y, fs] = pfetch(signal, channel, packetNumber);
%
% Parallel fetch: retrieve packets for array data.
%
% Input :  signal      : string variable denoting the signal (e.g. 'my_signal')
%          channel     : string variable denoting the acoustic channel (e.g. 'KAU1')
%          packetNumber: an integer value in the range [1, Nmax], where Nmax is the number
%                        of packets per channel file (bk.nPacketsPerSounding)
%        
% Output:  y  : multichannel time series of the distorted signal 
%          fs : sampling frequency
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
if packetNumber <= 0 | packetNumber > bk.nPacketsPerSounding
    fprintf(['Invalid packet number. Use a value between ' sprintf('%d', 1) ' and ' sprintf('%d', bk.nPacketsPerSounding) '.\r']);
    y=[]; fs=[];
    return
end

i_packet = packetNumber;
% Read array data
for i_sounding = 1:bk.nSoundings

    fullfile(waterMarkPath, 'output', channel, signal, [channel '_' sprintf('%03d', i_sounding)]);
    filename = fullfile(waterMarkPath, 'output', channel, signal, [channel '_' sprintf('%03d', i_sounding) '.wav']);
    %fprintf(['Fetching packet ' sprintf('%d', packetNumber) '/' sprintf('%d', bk.nPackets) '\r']);
    try
        [dummy, fs]=audioread(filename, [bk.packetIndex(i_packet,1), bk.packetIndex(i_packet,2)]);
        if i_sounding==1 % Dimensioning array
            y=zeros(length(dummy),bk.nSoundings);
        end
        data(:,i_sounding)=dummy;
    catch
        disp(['Error reading data from ' filename]);
        return
    end
end
    
%% A Doppler shift V0 has been removed from the data prior to the channel estimation.
% This has been done by resampling the raw acoustic data with a resampling factor 1-V0/c. 
% This Doppler shift will now re-instated. The used sign convention is such that a positive 
% velocity corresponds to a positive range rate, i.e. time dilation of the signal. 
% Small values may be entirely due to clock frequency offsets of transmitter and receiver
% hardware, which cause an apparent Doppler shift.  
% 
% NB. Resampling factors close to one may lead to memory problems, but the channel data 
% released with Watermark V1.0 should be fine fine on contemporary work stations.

c = 1500;  % nominal sound speed in m/s

% The Doppler velocity (m/s) should be the same for all channels for array data.
% Use the value given for the first hydrophone channel.
V0 = bk.V(1);  
resamplingFactor = 1/(1-V0/c);
[N1,D1] = rat(resamplingFactor);
try
    y = resample(data, N1, D1);
catch
    % Resampling factor is too close to 1. 
    disp('Warning: Doppler shift is not re-instated.');
end

return
