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
TxInputDataSampleRate=DataSampleRate;
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
RxInputDataSampleRate=DataSampleRate;
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
IncludeTX = 'Yes';
IncludeChannel = 'Yes';
IncludeRX = 'Yes';


% Analysis and vizaulization parameters (depends on simulation modes)
if strcmp(IncludeRX, 'Yes')
    PaMode='ReceiverOut';
elseif strcmp(IncludeTX, 'Yes')
    PaMode='TransmitterOut';
else
    PaMode='DataGeneratorOut';   
end
    
% validate parameters 
% check that the parameters for TX,CH, RX is compatible
% TBD


%%
% ======================================================================
% Create objects for simulation
% The program crate all objects in the link even part of this will be
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

% analyze link configuration
% TBD parsing of the configuration using this variables
% IncludeTX = 'Yes';
% IncludeChannel = 'Yes';
% IncludeRX = 'Yes';


%create signal to be transmitted

%create TX line signal

%apply channel

%receive the signal



%%
% Estimate the perfomance





