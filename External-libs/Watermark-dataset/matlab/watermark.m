function watermark(signal, channel, howmany)

% Synopsis: watermark(signal, channel, howmany);
%
% Input :  signal      : string variable denoting the signal (e.g. 'my_signal')
%          channel     : string variable denoting the channel (e.g. 'NOF1')
%          howmany     : use 'all' to process all channel files
%                      : use 'single' to process a single channel file
%                        
%
% This function takes an input signal, stacks multiple copies of it, and passes the
% stack through replay channels. Results are stored in Watermark's output directory. 
%
% Watermark version 1.0
% Forsvarets Forskningsinstitutt, 
% 03.11.2016


% Set path to Watermark base directory
P=mfilename('fullpath');
tmp = findstr(lower(P), 'matlab');
waterMarkPath=P(1:tmp(end)-1);

% load input signal
try
    load(fullfile(waterMarkPath, 'input', 'signals', signal));
catch
    disp(['Signal ' signal ' not found']);
    return
end

channel = upper(channel);

x=x(:); % ensure column vector
signalDuration=length(x)/fs_x;
effectiveBitRate=nBits/signalDuration;

fprintf('\r');
fprintf('Watermark V1.0\r');
fprintf('--------------\r');
fprintf(['Signal parameters for ' signal ':\r']);
fprintf(['message size = ' sprintf('%d', nBits) ' user bits\r']);
fprintf(['total duration = ' sprintf('%.3f', signalDuration) ' s\r']);
fprintf(['effective bit rate = ' sprintf('%.2f', effectiveBitRate) ' bit/s\r']);

fileList=dir(fullfile(waterMarkPath, 'input', 'channels', channel, 'mat', [channel '_*.mat']));
if isempty(fileList)
    fprintf(['Channel archive ' channel ' not found\r']);
    return
end

fprintf(['\rChannel parameters for ' channel ':\r']);

% Determine number of soundings to use
switch howmany
        case 'single'
            nSoundings=1;
        case 'all'
            nSoundings = length(fileList);
        otherwise
            disp('Invalid argument');
            return
end


% Load first file in the list and extract parameters
% All .mat files in this directory should be Watermark compatible channel soundings of the same type
load(fullfile(waterMarkPath, 'input', 'channels', channel, 'mat', [channel '_' sprintf('%03d',1)]));

dt = 1/fs_t; % time increment
dtau = 1/fs_tau; % delay increment
channelLength = (size(h,1)-1)*dt;
fprintf(['center frequency = ' sprintf('%.0f', fc) ' Hz\r']);
fprintf(['sounding duration = ' sprintf('%.1f', channelLength) ' s\r']);
fprintf(['number of soundings = ' sprintf('%d', nSoundings) '\r\r']);

% Determine maximum waveform length that 'fits into the channel' 
[s1 s2]=size(h);
maxSamples=floor((s1-1)*dt*fs_x);

Lx = length(x);

% The maximum delay spread in the replay channel is given by the delay coverage of h
L1 = ceil(s2*dtau*fs_x); 

% To be safe, add 2 ms of silence at either side of the signal. This
% accomodates possible sidelobes occurring in the filtering process.
L2 = ceil(2E-3*fs_x); 

% Zero padding of the signal
x_padded = [zeros(L2,1); x; zeros(L1+L2,1)];

% The number of packets that fits into the channel
nPacketsPerSounding = floor((maxSamples)/(Lx+L1+2*L2));

% Create a train of signal packets
signalTrain=repmat(x_padded, nPacketsPerSounding, 1);

% Compute the corresponding start and end samples of each packet
packetIndex=1+[0:nPacketsPerSounding-1]'*(Lx+L1+2*L2);
packetIndex=[packetIndex packetIndex+Lx+L1+2*L2-1];

nPackets = nPacketsPerSounding*nSoundings;
fprintf(['Number of packets per sounding = ' sprintf('%d', nPacketsPerSounding) '\r']);
fprintf(['Total number of simulated packets = ' sprintf('%d', nPackets) '\r\r']);

if nPackets==0
    maxSignalDuration=(maxSamples-L1-2*L2)/fs_x;
    fprintf(['The input signal is too long. It should be shorter than ~ ' sprintf('%3.1f', maxSignalDuration) ' s\r']);
    return
end

% Set path for the result files / create if necessary
outdir = fullfile(waterMarkPath, 'output', channel, signal);
if exist(outdir)~=7
      mkdir(outdir);
end

% Delete previous results, if any.
fprintf('Deleting previous results, if any ... ');
delete(fullfile(outdir, 'bookkeeping*.mat'));
delete(fullfile(outdir, [channel '_*.wav']));
fprintf('done.\r\r');


% Initiate vector with Doppler values
V = zeros(nSoundings,1);

% Apply the channel replay filter to all files in the channel archive
for n=1:nSoundings
    if n>=1
        load(fullfile(waterMarkPath, 'input', 'channels', channel, 'mat', [channel '_' sprintf('%03d',n)]));
    end
    
    V(n) = V0;
  
    fprintf(['Filtering ' signal ' x ' fileList(n).name(1:end-4) ' ... ' ]);
    y=replayfilter(signalTrain, fs_x, h, fs_t, fs_tau, fc);
    
    outfilename=fileList(n).name;
    outfilename(end-2:end)='wav';
    
    % Normalize to an RMS value of ~ 1 for the first file and store results while preserving 
    % relative signal levels.
    if n==1, normalizationFactor=sqrt(mean(y.^2)); end
    audiowrite(fullfile(outdir, outfilename), y/normalizationFactor, fs_x, 'BitsPerSample', 32);

    fprintf(['done. \r']);

    clear h fc V0
    
end

% Create and store bookkeeping struct
bk.nPackets = nPackets;
bk.nPacketsPerSounding = nPacketsPerSounding;
bk.packetIndex = packetIndex;
bk.nSoundings = nSoundings;
bk.nBits = nBits;
bk.effectiveBitRate=effectiveBitRate;
bk.V = V;
save(fullfile(outdir, 'bookkeeping.mat'), 'bk');

return
