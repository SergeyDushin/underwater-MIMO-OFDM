% =====================================================================    
% Hydroacoustic MIMO-OFDM system model
%
% Classes:
% dataGenerator - generator of data to be transmitted
% transmitter - class for transmitter side 
% receiver - class for receiver side
% channel - class for channel
% perfomanceAnalyzer - class for received data anlysis 
%
% Ojects:
% DG,TX, Ch, RX - object for link parts
% RA - object for perfomance analysis and vizualization
% =====================================================================    
% Sergey Dushin, nihsuds@gmail.com
% (c) 2021 by V.A. Trapeznikov Institute of control sciences of RAS
% www.ipu.ru
%
% =====================================================================
% MIT license
% =====================================================================
%version 0.1
%======================================================================
  
% clear
clear; close all;

%%
% ======================================================================
% Simulation parameters and its validation  
%=======================================================================
% Data generator and timings
DgSimulationTime=10;
DgDataSampleRate=100;
DgWaveForm='wgn';


% transmitter parameters
TxFrameLength=128;
TxCentralFrequency =35000;
TxQAMOrder=4;
TxNumberOfSubcarriers=256;
TxSubcarrierSpacing=20;
TxNumberOfAntenna=4;
TxInputDataSampleRate=DgDataSampleRate;
TxPreambulaLength=16;
TxCyclicPrefixLength=0.0001;
TxZeroGuardTimeLength=0;


% receiver  prarameters
RxFrameLength=128;
RxCentralFrequency =35000;
RxQAMOrder=4;
RxNumberOfSubcarriers=256;
RxNumberOfAntenna=4;
RxEqualizationMethod='OneTap';
RxInputDataSampleRate=DgDataSampleRate;
RxPreambulaLength=16;
RxCyclicPrefixLength=0.0001;
RxZeroGuardTimeLength=0;


% channel parameters
ChNumberOfInput=4;
ChNumberOfOutput=4;
ChType='Watermark_BCH1';


% Choose link configurations.
% Channel,TX,RX can be excluded from link for dedug 
% Data generator and analyzer are always presented
% Relevant configurations:
% DataGenerator-->analyzer
% DataGenerator-->TX-->analyzer
% DataGenerator-->TX-->Channel-->analyzer
% DataGenerator-->TX-->RX-->analyzer
% DataGenerator-->TX-->Channel-->RX-->analyzer

% Choose parts of link witch you can add for simulation:
%IncludeTX = 'Yes';
IncludeTX = 'No';
%IncludeChannel = 'Yes';
IncludeChannel = 'No';
%IncludeRX = 'Yes';
IncludeRX = 'No'; % TX, Ch, RX are off - for generator-->analyzer dedug

%Validate and print the selected link configuration: 
if strcmp(IncludeTX,'No')&&strcmp(IncludeChannel,'No')&&strcmp(IncludeRX,'No')
    disp('Link configuration: DataGenerator-->analyzer');
elseif strcmp(IncludeTX,'Yes')&&strcmp(IncludeChannel,'No')&&strcmp(IncludeRX,'No')
    disp('Link configuration:DataGenerator-->TX-->analyzer');
elseif strcmp(IncludeTX,'Yes')&&strcmp(IncludeChannel,'Yes')&&strcmp(IncludeRX,'No')
    disp('Link configuration:DataGenerator-->TX-->Channel-->analyzer');
elseif strcmp(IncludeTX,'Yes')&&strcmp(IncludeChannel,'No')&&strcmp(IncludeRX,'Yes')
    disp('Link configuration:DataGenerator-->TX-->RX-->analyzer');
elseif strcmp(IncludeTX,'Yes')&&strcmp(IncludeChannel,'Yes')&&strcmp(IncludeRX,'Yes')
    disp('Link configuration:DataGenerator-->TX-->Channel-->RX-->analyzer');
else
    error ('Invalid link configuration. You cannot switch on Ch and/or RX if TX is off')
end

% validate parameters 
% check that the parameters for DG, TX, CH, RX, Analyzer are compatible
% TBD


%%
% ======================================================================
% Create objects for simulation
% The program crate all objects in the link even some of they will be
% skiped during simulation
%=======================================================================
% data
DG=dataGenerator(...
    DgSimulationTime,...
    DgDataSampleRate,...
    DgWaveForm...
    ); 

% transmitter
TX=transmitter (...
    TxFrameLength,...
    TxCentralFrequency,...
    TxQAMOrder,...
    TxNumberOfSubcarriers,...
    TxSubcarrierSpacing,...
    TxNumberOfAntenna,...
    TxInputDataSampleRate,...
    TxPreambulaLength,...
    TxCyclicPrefixLength,...
    TxZeroGuardTimeLength...
    ); 

% channel
Ch=channel(...
    ChNumberOfInput,...
    ChNumberOfOutput,...
    ChType...
    );

% receiver
RX=receiver(...
    RxFrameLength,...
    RxCentralFrequency,...
    RxQAMOrder,...
    RxNumberOfSubcarriers,...
    RxNumberOfAntenna,...
    RxEqualizationMethod,...
    RxInputDataSampleRate,...
    RxPreambulaLength,...
    RxCyclicPrefixLength,...
    RxZeroGuardTimeLength...
    );

%%
% ======================================================================
% Run simulation   
%=======================================================================

% No need any additional validation. It is performed in "parameters" section

% Start the simulation

%create signal to be transmitted. It is created always
[InputData,GenFs]=DG.generateData();

%create TX passband signal if defined
if strcmp(IncludeTX, 'Yes')
    TxOut=TX.TransmitData(InputData);
end

%apply channel if defined
if strcmp(IncludeChannel, 'Yes')
    ChOut=CH.ApplyChannel;
end

%receive the signal if defined
if strcmp(IncludeRX, 'Yes')
    RxOut=RX.ReceiveTheSignal;
end

%%
%=======================================================================
% Create and run perfomance analyzer
% The first version of analyzer receive only Output signal. 
% Do it more flaxeble in next versions (for receiving of any set of the signal).  
%=======================================================================

% Analyzer Settings
if strcmp(IncludeRX, 'Yes')
    PaMode='ReceiverOut';
elseif strcmp(IncludeTX, 'Yes')
    PaMode='TransmitterOut';
else
    PaMode='DataGeneratorOut';   
end

% Assign Link Output for analyzer.
switch PaMode
    case 'DataGeneratorOut'
        SimOut=InputData;
        SampleRate=GenFs;
    case 'TransmitterOut'
        SimOut=TxOut;
    case 'Channel'
        SimOut=ChOut;
    case 'ReceiverOut'
        SimOut=RxOut;
    otherwise
        error('Simulation Output is invalid');
end

%Create Analyzer object
RA=perfomanceAnalyzer(...
    PaMode...
    );

% Estimate the perfomance. It is performed always (and depend on PaMode)
RA.AnalyzeData(InputData,SimOut,SampleRate,DgSimulationTime);
disp('Program done');





