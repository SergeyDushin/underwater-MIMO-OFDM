% =====================================================================    
% Hydroacoustic MIMO-OFDM implemetation
%
% Classes:
% transmitter - class for transmitter side 
% receiver - class for receiver side
% channel - class for channel
%
% Ojects:
% TX, Ch, RX - object for link parts 
% =====================================================================    
% Sergey Dushin, nihsuds@gmail.com
% (c) 2021 by V.A. Trapeznikov Institute of control sciences of RAS
% www.ipu.ru
% =====================================================================
  
% paths
clear; close all;
% addpath('./External-libs/Watermark-dataset/BCH1/mat');
% addpath('./External-libs/Watermark-dataset/matlab');
% addpath('./External-libs/Nissel-FBMC-OFDM');
% addpath('./External-libs/Nissel-FBMC-OFDM/Theory');


% ======================================================================
% Simulation parameters   
%=======================================================================
% Simulation time and input data
SimulationTime=10;
DataSampleRate=100000;
DataLength=SimulationTime*DataSampleRate;

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
RxEqualizationMethod='OneTap'
RxInputDataSampleRate=DataSampleRate;
RxPreambulaLength=16;
RxCyclicPrefixLength=0.0001;
RxZeroGuardTimeLength=0;



% channel parameters
% TBD

% validate parameters 
% TBD

% Here is different modulation modes
% TX-->analizer
% TX-->Channel-->analizer
% TX-->RX-->analizer
% TX-->Channel --> analizer

% ======================================================================
% Create data and objects   
%=======================================================================
% data
data=data_generator(DataLength); 

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
%Ch=channel(rx_parameters);

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


% ======================================================================
% Run simulation   
%=======================================================================

%create transmit signal


%apply channel


%receive signal


%estimate the perfomznce





