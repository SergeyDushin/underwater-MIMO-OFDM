% =====================================================================    
% Hydroacoustic MIMO-OFDM implemetation
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
DgDataSampleRate=100000;
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



% Here is different link configurations:
% Channel,TX,RX can be excluded from link for dedug 
% Data generator and analyzer allways presented
% Relevant configuration:
% DataGenerator-->analyzer
% DataGenerator-->TX-->analyzer
% DataGenerator-->TX-->Channel-->analyzer
% DataGenerator-->TX-->RX-->analyzer
% DataGenerator-->TX-->Channel-->RX-->analyzer
%IncludeTX = 'Yes';
IncludeTX = 'No';
%IncludeChannel = 'Yes';
IncludeChannel = 'No';
%IncludeRX = 'Yes';
IncludeRX = 'No'; % TX, Ch, RX are off - for generator-->analyzer dedug


% Analysis and vizaulization parameters (depends on simulation modes)
% Channel Out = TX Out in terms of analyzer
if strcmp(IncludeRX, 'Yes')
    PaMode='ReceiverOut';
elseif strcmp(IncludeTX, 'Yes')
    PaMode='TransmitterOut';
else
    PaMode='DataGeneratorOut';   
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

RA=perfomanceAnalyzer(...
    PaMode...
    );


%%
% ======================================================================
% Run simulation   
%=======================================================================

% validate defined simulation mode to exclude wrong link configurations 
if strcmp(IncludeTX, 'No')
    if strcmp(IncludeRX, 'Yes')
        error('Invalid simulation mode: RX must be off, when TX is off')
    elseif strcmp(IncludeChannel, 'Yes')
        error('Invalid simulation mode: Channel must be off, when TX is off')
    end 
end

% Start the simulation

%create signal to be transmitted. It is created always
InputData=DG.generateData();
SimOut=InputData; % default simulation output if other parts are off

%create TX passband signal if defined
if strcmp(IncludeTX, 'Yes')
    %TX.TransmittData; maybe several methods for transmitted signal creation
    SimOut=TxOut;
end

%apply channel if defined
if strcmp(IncludeChannel, 'Yes')
    %CH.ApplyChannel; maybe several methods to apply channel
    SimOut=ChOut;
end

%receive the signal if defined
if strcmp(IncludeRX, 'Yes')
    %RX.ReceiveTheSignal; maybe several methods for receive the signal
    SimOut=RxOut;
end


% Estimate the perfomance. It is performed always (and depend on PaMode)
RA.AnalyzeData(InputData,SimOut);





