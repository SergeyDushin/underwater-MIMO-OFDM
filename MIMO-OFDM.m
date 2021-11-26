% =====================================================================    
% Данная программа реализует эксперимент, описанный в докладе ДКАШФ на EnT 2020
% Сравнивается разные типы многочастотной модуляции в применении к
% гидроакустическим каналам. В том числе рассматриваются: 
%       1) OFDM with CP (worst spectral behaviour)
%       2) FBMC (best spectral properties, complex orthogonality is replaced by real orthogonality)
%       3) WOLA (windowed OFDM. The windowing is done at TX and RX)
%       4) FOFDM (filtered OFDM, sinc + Hann, filtering at TX and RX)
%       5) UFMC (filtered OFDM, subband-wise filtering, Dolph-Chebyshev window, cyclic prefix and zero padding are supported, filtering at TX and RX)
% В программе используются сторонние библиотеки от Ronald Nissel, rnissel@nt.tuwien.ac.at
% и датасет Watermark в качестве канала связи https://ieeexplore.ieee.org/document/7932436
% ===================================================================== 


% Пути
clear; close all;
addpath('./External-libs/Watermark-dataset/BCH1/mat');
addpath('./External-libs/Watermark-dataset/matlab');
addpath('./External-libs/Nissel-FBMC-OFDM');
addpath('./External-libs/Nissel-FBMC-OFDM/Theory');

% Параметры для модели канала BCH1                               
% Загрузка импульсной характеристики BCH1 первый элемент
load BCH1_001.mat; 



% Параметры симуляции исходя из полосы пропускания 5 кГц для BCH1    
QAM_ModulationOrder          =4;                                      % QPSK для нашего случая.
Number_of_carriers = 256;
Subcarrier_Spacing = 20;
SamplingRate = Number_of_carriers*Subcarrier_Spacing*4;
dt = 1/SamplingRate;
Simulation_SNR_OFDM_dB = [5:2.5:25];                        % SNR for OFDM in dB. The average transmit power of all methods is the same! However, the SNR might be different due to filtering (in FOFDM and UFMC) or because a different bandwidth is used (different subcarrier spacing or different number of subcarriers).
NumberOfSymbolsInTime= 40; % количетво передаваемых фреймов


f_central =35000;

% Параметры исследуемых многочастотных схем

% FBMC parameters
FBMC_NumberOfSubcarriers     = Number_of_carriers;                                      % Number of subcarriers
FBMC_NumberOfSymbolsInTime   = NumberOfSymbolsInTime;                                      % Number FBMC symbols in time
FBMC_SubcarrierSpacing       = Subcarrier_Spacing;                                    % Subcarrier spacing (Hz)
FBMC_PrototypeFilter         = 'Hermite-QAM';                          % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM.
FBMC_OverlappingFactor       = 4;                                       % Overlapping factor, 2,3,4,...

% OFDM parameters
OFDM_NumberOfSubcarriers     =Number_of_carriers;                                      % Number of subcarriers 
OFDM_NumberOfSymbolsInTime   = NumberOfSymbolsInTime;                                      % Number OFDM symbols in time
OFDM_SubcarrierSpacing       = Subcarrier_Spacing;                                    % Subcarrier spacing (Hz)
OFDM_CyclicPrefixLength      = 1/(14*OFDM_SubcarrierSpacing);           % Length of the cyclic prefix (s)

% WOLA parameters
WOLA_NumberOfSubcarriers     = Number_of_carriers;                                      % Number of subcarriers                                   
WOLA_NumberOfSymbolsInTime   = NumberOfSymbolsInTime;                                      % Number WOLA symbols in time  
WOLA_SubcarrierSpacing       = Subcarrier_Spacing;                                    % Subcarrier spacing (Hz)
WOLA_CyclicPrefixLength      = 0;                                       % Length of the cyclic prefix (s) to combat the channel
WOLA_WindowLengthTX          = 1/(14*2*WOLA_SubcarrierSpacing);         % Length of the window overlapping (s) at the transmitter 
WOLA_WindowLengthRX          = 1/(14*2*WOLA_SubcarrierSpacing);         % Length of the window overlapping (s) at the receiver

% FOFDM parameters
FOFDM_NumberOfSubcarriers     = Number_of_carriers;                                     % Number of subcarriers
FOFDM_NumberOfSymbolsInTime   = NumberOfSymbolsInTime;                                     % Number FOFDM symbols in time                        
FOFDM_SubcarrierSpacing       = Subcarrier_Spacing;                                   % Subcarrier spacing (Hz)
FOFDM_CyclicPrefixLength      = 0;                                      % Length of the cyclic prefix (s) to combat the channel
FOFDM_FilterLengthTX          = 0.2*1/(FOFDM_SubcarrierSpacing);        % Length of the transmit filter (s)
FOFDM_FilterLengthRX          = 0.2*1/(FOFDM_SubcarrierSpacing);        % Length of the receive filter (s) 
FOFDM_FilterCylicPrefixLength = 1/(14*FOFDM_SubcarrierSpacing);         % Length of the additional cyclic prefix (s) to combat ISI and ICI due to the filtering

% UFMC parameters
UFMC_NumberOfSubcarriers     = Number_of_carriers;                                      % Number of subcarriers
UFMC_NumberOfSymbolsInTime   = NumberOfSymbolsInTime;                                      % Number UFMC symbols in time
UFMC_SubcarrierSpacing       = Subcarrier_Spacing;                                    % Subcarrier spacing (Hz)
UFMC_CyclicPrefixLength      = 0;                                       % Length of the cyclic prefix (s) to combat the channel. If zero padding is used, this length reprents the zero guard length instead of the CP length.
UFMC_FilterLengthTX          = 1/14*1/(UFMC_SubcarrierSpacing);         % Length of the transmit filter (s)
UFMC_FilterLengthRX          = 1/14*1/(UFMC_SubcarrierSpacing);         % Length of the receive filter (s)
UFMC_FilterCylicPrefixLength = 1/(14*UFMC_SubcarrierSpacing);           % Length of the additional cyclic prefix (or zero guard symbol if ZP is used) in seconds (s). Needed to combat ISI and ICI due to the filtering. However, small ICI and ISI is perfectly feasibly.
UFMC_ZeroPaddingInsteadOfCP  = true;                                    % TRUE for Zero Padding (ZP) and FALSE for a conventional Cyclic Prefix (CP). Note that a CP delivers nicer plots of the power spectral density because there are no zero crossing.



%% Generate " +Modulation\" Objects
% FBMC Object
FBMC = Modulation.FBMC(...
    FBMC_NumberOfSubcarriers,...                                        % Number of subcarriers
    FBMC_NumberOfSymbolsInTime,...                                      % Number FBMC symbols in time
    FBMC_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz).  Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    FBMC_PrototypeFilter,...                                            % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM. The data rate of QAM is reduced by a factor of two compared to OQAM, but robustness in doubly-selective channels is inceased
    FBMC_OverlappingFactor, ...                                         % Overlapping factor (also determines oversampling in the frequency domain)                                   
    0, ...                                                              % Initial phase shift
    true ...                                                            % Polyphase implementation
    );
FBMC_BlockOverlapTime = (FBMC.PrototypeFilter.OverlappingFactor-1/2)*FBMC.PHY.TimeSpacing;

% OFDM Object
OFDM = Modulation.OFDM(...
    OFDM_NumberOfSubcarriers,...                                        % Number of subcarriers
    OFDM_NumberOfSymbolsInTime,...                                      % Number OFDM symbols in time                                                 
    OFDM_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz) 
    SamplingRate,...                                                    % Sampling rate (Samples/s)                                       
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    OFDM_CyclicPrefixLength, ...                                        % Length of the cyclic prefix (s)                 
    FBMC_BlockOverlapTime ...                                           % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    );

% Windowed OFDM (WOLA)
WOLA = Modulation.WOLA(...
    WOLA_NumberOfSubcarriers,...                                        % Number subcarriers
    WOLA_NumberOfSymbolsInTime,...                                      % Number WOLA symbols in time 
    WOLA_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    0, ...                                                              % Length of the cyclic prefix (s)
    FBMC_BlockOverlapTime-WOLA_WindowLengthTX/2, ...                    % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    WOLA_WindowLengthTX, ...                                            % Length of the window overlapping (s) at the transmitter 
    WOLA_WindowLengthRX ...                                             % Length of the window overlapping (s) at the receiver
    );

% FOFDM (Filtered OFDM)
FOFDM = Modulation.FOFDM(...
    FOFDM_NumberOfSubcarriers,...                                       % Number of subcarriers
    FOFDM_NumberOfSymbolsInTime,...                                     % Number FOFDM symbols in time                
    FOFDM_SubcarrierSpacing,...                                         % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    0, ...                                                              % Length of the cyclic prefix (s)
    FBMC_BlockOverlapTime-FOFDM_FilterLengthTX/2, ...                   % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission                    
    FOFDM_FilterLengthTX, ...                                           % Length of the transmit filter (s)
    FOFDM_FilterLengthRX, ...                                           % Length of the receive filter (s) 
    FOFDM_FilterCylicPrefixLength ...                                   % Length of the additional cyclic prefix (s).  Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
);

% UFMC (Subband Filtered OFDM)
UFMC = Modulation.UFMC(...
    UFMC_NumberOfSubcarriers,...                                        % Number of subcarriers
    UFMC_NumberOfSymbolsInTime,...                                      % Number UFMC symbols in time
    UFMC_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    0, ...                                                              % Length of the cyclic prefix (s). If zero padding is used, this length reprents the zero guard length instead of the CP length
    FBMC_BlockOverlapTime-UFMC_FilterLengthTX/2, ...                    % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    UFMC_FilterLengthTX, ...                                            % Length of the transmit filter (s)
    UFMC_FilterLengthRX, ...                                            % Length of the receive filter (s)
    UFMC_FilterCylicPrefixLength, ...                                   % Length of the additional cyclic prefix (or zero guard symbol if ZP is used) in seconds (s). Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
    UFMC_ZeroPaddingInsteadOfCP ...                                     % TRUE for Zero Padding (ZP) or FALSE for a conventional Cyclic Prefix (CP)
);

% Number of samples
N_FBMC  = FBMC.Nr.SamplesTotal;
N_OFDM  = OFDM.Nr.SamplesTotal;
N_WOLA  = WOLA.Nr.SamplesTotal;
N_FOFDM = FOFDM.Nr.SamplesTotal;
N_UFMC  = UFMC.Nr.SamplesTotal;
N       = max([N_FBMC N_OFDM N_WOLA N_FOFDM N_UFMC]);


% PAM and QAM Object
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
if strcmp(FBMC.Method(end-3),'O')
    % FBMC-OQAM transmission, only real valued data symbols
    PAMorQAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
else
    % FBMC-QAM transmission,  complex valued data symbols
    PAMorQAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
end


% Pre-allocate Transmit Power
Ps_FBMC   = zeros(N_FBMC,1);
Ps_OFDM   = zeros(N_OFDM,1);
Ps_WOLA   = zeros(N_WOLA,1);
Ps_FOFDM  = zeros(N_FOFDM,1);
Ps_UFMC   = zeros(N_UFMC,1);

% Pre-allocate Power Spectral Density
PSD_FBMC  = zeros(N_FBMC,1);
PSD_OFDM  = zeros(N_OFDM,1);
PSD_WOLA  = zeros(N_WOLA,1);
PSD_FOFDM = zeros(N_FOFDM,1);
PSD_UFMC  = zeros(N_UFMC,1);

%% Симуляция
% Binary data
BinaryDataStream_FBMC  = randi([0 1], FBMC.Nr.Subcarriers  * FBMC.Nr.MCSymbols  * log2(PAMorQAM.ModulationOrder),1);
BinaryDataStream_OFDM  = randi([0 1], OFDM.Nr.Subcarriers  * OFDM.Nr.MCSymbols  * log2(QAM.ModulationOrder),1);
BinaryDataStream_WOLA  = randi([0 1], WOLA.Nr.Subcarriers  * WOLA.Nr.MCSymbols  * log2(QAM.ModulationOrder),1);
BinaryDataStream_FOFDM = randi([0 1], FOFDM.Nr.Subcarriers * FOFDM.Nr.MCSymbols * log2(QAM.ModulationOrder),1);
BinaryDataStream_UFMC  = randi([0 1], FOFDM.Nr.Subcarriers * FOFDM.Nr.MCSymbols * log2(QAM.ModulationOrder),1);

% Transmitted data symbols (map binary data to symbol)
x_FBMC  = reshape( PAMorQAM.Bit2Symbol(BinaryDataStream_FBMC)  , FBMC.Nr.Subcarriers  , FBMC.Nr.MCSymbols);
x_OFDM  = reshape(      QAM.Bit2Symbol(BinaryDataStream_OFDM)  , OFDM.Nr.Subcarriers  , OFDM.Nr.MCSymbols);
x_WOLA  = reshape(      QAM.Bit2Symbol(BinaryDataStream_WOLA)  , WOLA.Nr.Subcarriers  , WOLA.Nr.MCSymbols);
x_FOFDM = reshape(      QAM.Bit2Symbol(BinaryDataStream_FOFDM) , FOFDM.Nr.Subcarriers , FOFDM.Nr.MCSymbols);
x_UFMC  = reshape(      QAM.Bit2Symbol(BinaryDataStream_UFMC)  , UFMC.Nr.Subcarriers  , UFMC.Nr.MCSymbols);

% Transmitted signal in the time domain
s_FBMC  =  FBMC.Modulation( x_FBMC );
s_OFDM  =  OFDM.Modulation( x_OFDM );
s_WOLA  =  WOLA.Modulation( x_WOLA );
s_FOFDM = FOFDM.Modulation( x_FOFDM );
s_UFMC  =  UFMC.Modulation( x_UFMC );

% Свертка сигнала с нестационарным каналом в основной полосе 
% Используется модифицированный интсрумет из Вотермарк replayfilter_bb
r_FBMC_noNoise  = replayfilter_bb(s_FBMC, SamplingRate, h, fs_t, fs_tau, f_central);
r_OFDM_noNoise  = replayfilter_bb(s_OFDM, SamplingRate, h, fs_t, fs_tau, f_central);
r_WOLA_noNoise  = replayfilter_bb(s_WOLA, SamplingRate, h, fs_t, fs_tau, f_central);
r_FOFDM_noNoise = replayfilter_bb(s_FOFDM, SamplingRate, h, fs_t, fs_tau, f_central);
r_UFMC_noNoise  = replayfilter_bb(s_UFMC, SamplingRate, h, fs_t, fs_tau, f_central);



% Calculate the transmitted power over time
Ps_FBMC  = Ps_FBMC  + abs(s_FBMC).^2;
Ps_OFDM  = Ps_OFDM  + abs(s_OFDM).^2;
Ps_WOLA  = Ps_WOLA  + abs(s_WOLA).^2;
Ps_FOFDM = Ps_FOFDM + abs(s_FOFDM).^2;
Ps_UFMC  = Ps_UFMC  + abs(s_UFMC).^2;

% Calculat the power spectral density
PSD_FBMC  = PSD_FBMC  + abs(fft(s_FBMC)/sqrt(N_FBMC)).^2;
PSD_OFDM  = PSD_OFDM  + abs(fft(s_OFDM)/sqrt(N_OFDM)).^2;
PSD_WOLA  = PSD_WOLA  + abs(fft(s_WOLA)/sqrt(N_WOLA)).^2;
PSD_FOFDM = PSD_FOFDM + abs(fft(s_FOFDM)/sqrt(N_FOFDM)).^2;
PSD_UFMC  = PSD_UFMC  + abs(fft(s_UFMC)/sqrt(N_UFMC)).^2;   



%R_vecH = E(h(:)*h(:)'); % математическое ожидание автокорреляции 

    
    for i_SNR = 1:length(Simulation_SNR_OFDM_dB)
        % Add noise
        SNR_OFDM_dB = Simulation_SNR_OFDM_dB(i_SNR);
        Pn_time     = 1/OFDM.GetSymbolNoisePower(1)*10^(-SNR_OFDM_dB/10);
        noise       = sqrt(Pn_time/2)*(randn(N,1)+1j*randn(N,1));
            
        r_FBMC  = r_FBMC_noNoise  + noise(1:N_FBMC);
        % Во всех сигналах кроме FBMC почему то появлялся лишний отсчет
        % после прохождения канала
        r_OFDM  = r_OFDM_noNoise(1:N_OFDM)  + noise(1:N_OFDM);
        r_WOLA  = r_WOLA_noNoise(1:N_WOLA)  + noise(1:N_WOLA);
        r_FOFDM = r_FOFDM_noNoise(1:N_FOFDM) + noise(1:N_FOFDM);
        r_UFMC  = r_UFMC_noNoise(1:N_UFMC)  + noise(1:N_UFMC);
       
        % Демодуляция сигналов
        y_FBMC  =  FBMC.Demodulation(r_FBMC);
        y_OFDM  =  OFDM.Demodulation(r_OFDM);
        y_WOLA  =  WOLA.Demodulation(r_WOLA);
        y_FOFDM = FOFDM.Demodulation(r_FOFDM);
        y_UFMC  =  UFMC.Demodulation(r_UFMC);
        
        
%         %Параметры эквалайзера
%         PIL=20;
%         P=1;
%         M=1;
%         S=1;
%         F=1;
%         T=1;
%         Equalizer_length =3; 
%         set_adaptive_alg_param;     % set the parameters for the adaptive algorithm
%         set_plot_param;
%         
       % эквалайзер 
%        y_Equalized_FBMC  = Simple_Equalazer(x_FBMC, y_FBMC,'E-FFT',PIL, LMS.G, LMS.E,crr_grd,which_plots,Equalizer_length);
%        y_Equalized_OFDM  = Simple_Equalazer( x_OFDM  , y_OFDM,'E-FFT',PIL, LMS.G, LMS.E,crr_grd,which_plots,Equalizer_length);
%        y_Equalized_WOLA  = Simple_Equalazer( x_WOLA  , y_WOLA,'E-FFT',PIL, LMS.G, LMS.E,crr_grd,which_plots,Equalizer_length);
%        y_Equalized_FOFDM = Simple_Equalazer( x_FOFDM , y_FOFDM,'E-FFT',PIL, LMS.G, LMS.E,crr_grd,which_plots,Equalizer_length);
%        y_Equalized_UFMC  = Simple_Equalazer( x_UFMC  , y_UFMC,'E-FFT',PIL, LMS.G, LMS.E,crr_grd,which_plots,Equalizer_length);
%         
        
%         % Простой эквалайзер
       y_Equalized_FBMC  = Simple_Equalazer(x_FBMC, y_FBMC,h);
       y_Equalized_OFDM  = Simple_Equalazer( x_OFDM  , y_OFDM,h);
       y_Equalized_WOLA  = Simple_Equalazer( x_WOLA  , y_WOLA,h);
       y_Equalized_FOFDM = Simple_Equalazer( x_FOFDM , y_FOFDM,h);
       y_Equalized_UFMC  = Simple_Equalazer( x_UFMC  , y_UFMC,h);


  
        
        % Detect symbols (quantization and demapping to bits)
        DetectedBitStream_Equalized_FBMC  = PAMorQAM.Symbol2Bit(y_Equalized_FBMC(:));
        DetectedBitStream_Equalized_OFDM  = QAM.Symbol2Bit(y_Equalized_OFDM(:));
        DetectedBitStream_Equalized_WOLA  = QAM.Symbol2Bit(y_Equalized_WOLA(:));
        DetectedBitStream_Equalized_FOFDM = QAM.Symbol2Bit(y_Equalized_FOFDM(:));
        DetectedBitStream_Equalized_UFMC  = QAM.Symbol2Bit(y_Equalized_UFMC(:));
       
        % Calculate the BER
        BER_FBMC_Equalized(i_SNR,1)   = mean( BinaryDataStream_FBMC~=DetectedBitStream_Equalized_FBMC );
        BER_OFDM_Equalized(i_SNR,1)   = mean( BinaryDataStream_OFDM~=DetectedBitStream_Equalized_OFDM );
        BER_WOLA_Equalized(i_SNR,1)   = mean( BinaryDataStream_WOLA~=DetectedBitStream_Equalized_WOLA );
        BER_FOFDM_Equalized(i_SNR,1)  = mean( BinaryDataStream_FOFDM~=DetectedBitStream_Equalized_FOFDM );
        BER_UFMC_Equalized(i_SNR,1)   = mean( BinaryDataStream_UFMC~=DetectedBitStream_Equalized_UFMC ); 
        
        
        % Расчет SIR
     
%     [PS_FBMC,PI_FBMC] = FBMC.GetSignalAndInterferencePowerOQAM(...
%                 R_vecH,...                                                  % Let the received signal be r=H*s with H representing a time-variant convolution matrix. Then "VectorizedChannelCorrelationMatrix" represents the expectation E{{H(:)*H(:)'}. We can obtain such matrix by ChannelModel.GetCorrelationMatrix
%                 eye(FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols),...              % Correlation matrix of the vectorized data symbols
%                 0,...                                                       % Time offset in samples (to improve the SIR)
%                 round(FBMC.Nr.Subcarriers/2),...                            % Subcarrier position for which the SIR is calculated.
%                 round(FBMC.Nr.MCSymbols/2)...                               % FBMC symbol position in time for which the SIR is calculated.
%                 );                          
%      
%     [PS_OFDM,PI_OFDM] = OFDM.GetSignalAndInterferencePowerQAM(...
%                 R_vecH,...                                                  % Let the received signal be r=H*s with H representing a time-variant convolution matrix. Then "VectorizedChannelCorrelationMatrix" represents the expectation E{{H(:)*H(:)'}. We can obtain such matrix by ChannelModel.GetCorrelationMatrix
%                 eye(OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols),...              % Correlation matrix of the vectorized data symbols
%                 0,...                                                       % Time offset in samples (to improve the SIR)
%                 round(OFDM.Nr.Subcarriers/2),...                            % Subcarrier position for which the SIR is calculated.
%                 round(OFDM.Nr.MCSymbols/2)...                               % OFDM symbol position in time for which the SIR is calculated.
%                 );
%             
%             
%             
%         SIR_FBMC_Object         = 10*log10( PS_FBMC      /  PI_FBMC);
%         SIR_OFDM_Object         = 10*log10( PS_OFDM      /  PI_OFDM);
%         
%         
%         
%         SIR_FBMC(i_SNR,1)   = mean( BinaryDataStream_FBMC~=DetectedBitStream_Equalized_FBMC );
%         SIR_OFDM(i_SNR,1)   = mean( BinaryDataStream_OFDM~=DetectedBitStream_Equalized_OFDM );
%         SIR_WOLA(i_SNR,1)   = mean( BinaryDataStream_WOLA~=DetectedBitStream_Equalized_WOLA );
%         SIR_FOFDM(i_SNR,1)  = mean( BinaryDataStream_FOFDM~=DetectedBitStream_Equalized_FOFDM );
%         SIR_UFMC(i_SNR,1)   = mean( BinaryDataStream_UFMC~=DetectedBitStream_Equalized_UFMC );
    
    end
   

% Define colors for different modulation schemes
ColorFBMC   = [0 0 1]*0.5;
ColorOFDM   = [1 0 0];
ColorWOLA   = [1 0 1];
ColorFOFDM  = [1 1 0]*0.7;
ColorUFMC   = [0 1 0]*0.5;

%% Plot Stuff
% Plot BER
figure(1);
semilogy( Simulation_SNR_OFDM_dB , mean(BER_FBMC_Equalized,2)  ,'x-','color',ColorFBMC);  hold on;
semilogy( Simulation_SNR_OFDM_dB , mean(BER_OFDM_Equalized,2)  ,'o-','color',ColorOFDM);  hold on;
semilogy( Simulation_SNR_OFDM_dB , mean(BER_WOLA_Equalized,2)  ,'s-','color',ColorWOLA);  hold on;
semilogy( Simulation_SNR_OFDM_dB , mean(BER_FOFDM_Equalized,2) ,'d-','color',ColorFOFDM); hold on;
semilogy( Simulation_SNR_OFDM_dB , mean(BER_UFMC_Equalized,2)  ,'*-','color',ColorUFMC);  hold on;
ylabel('Bit Error Ratio');
xlabel('Signal-to-Noise Ratio at the receiver [dB]');
legend({'FBMC','OFDM','WOLA','f-OFDM','UFMC'});

% отношение SIR
%figure(2);


% Plot spectrum
figure(3);
frequency_FBMC  = (0:N_FBMC  -1) / (N_FBMC *dt);
frequency_OFDM  = (0:N_OFDM  -1) / (N_OFDM *dt);
frequency_WOLA  = (0:N_WOLA  -1) / (N_WOLA *dt);
frequency_FOFDM = (0:N_FOFDM -1) / (N_FOFDM*dt);
frequency_UFMC  = (0:N_UFMC  -1) / (N_UFMC *dt);
plot( frequency_FBMC /1e6 , 10*log10(PSD_FBMC/max(PSD_FBMC)) ,'color',ColorFBMC);  hold on;
plot( frequency_OFDM /1e6 , 10*log10(PSD_OFDM/max(PSD_FBMC)) ,'color',ColorOFDM);  hold on;
plot( frequency_WOLA /1e6 , 10*log10(PSD_WOLA/max(PSD_FBMC)) ,'color',ColorWOLA);  hold on;
plot( frequency_FOFDM/1e6 , 10*log10(PSD_FOFDM/max(PSD_FBMC)),'color',ColorFOFDM); hold on;
plot( frequency_UFMC /1e6 , 10*log10(PSD_UFMC/max(PSD_FBMC)) ,'color',ColorUFMC);  hold on;
ylabel('Power Spectral Density [dB]');
xlabel('Frequency [MHz]');
ylim([-100 0]);
xlim([0 max([frequency_FBMC(round(end/2)) frequency_OFDM(round(end/2)) frequency_WOLA(round(end/2)) frequency_FOFDM(round(end/2)) frequency_UFMC(round(end/2))]/1e6)]);
legend({'FBMC','OFDM','WOLA','f-OFDM','UFMC'});
title('Simulation');


%% impulse response
figure (4);

tau=[0:length(h(1,:))-1]/fs_tau;
plot(tau, real(h(20,:)));
xlabel('Time delay [с]', 'fontsize', 16), ylabel('Instantaneous impulse response BCH1 (ch.1)', 'fontsize', 16);




