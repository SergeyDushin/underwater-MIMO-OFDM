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

% paths
clear; close all;
addpath('./External-libs/Watermark-dataset/BCH1/mat');
addpath('./External-libs/Watermark-dataset/matlab');
addpath('./External-libs/Nissel-FBMC-OFDM');
addpath('./External-libs/Nissel-FBMC-OFDM/Theory');


% ======================================================================
% Simulation parameters   
%=======================================================================
% comon and data parameters
data_length=2048;


% transmitter parameters

% receiver prarameters

% validate parameters 


% ======================================================================
% Create data and objects   
%=======================================================================
% data
data=data_generator(data_length); 

% transmitter
TX=transmitter(tx_parameters);

% channel
Ch=channel(rx_parameters);

% receiver
RX=receiver(rx_parameters);


% ======================================================================
% Run simulation   
%=======================================================================

%create transmit signal


%apply channel


%receive signal


%estimate the perfomznce





