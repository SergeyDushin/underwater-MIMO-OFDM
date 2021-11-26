classdef BCH1 < handle
    % =====================================================================        
    % Этот MATLAB класс загружает BCH1 эксперимент  
    % The channel parameters are initialized by the class
    % constructor. The convolution is performed by the method
    % ".Convolution(.)". A new channel realization is obtained by 
    % ".NewRealization;" 
    
    
    properties (SetAccess = private)
        PHY
        Nr
        Implementation
        ImpulseResponse
    end
    
    methods
        function obj = BCH1(...
                SamplingRate,...                                   % Sampling rate (Samples/s)
                SamplesTotal...                                   % Number of total samples   
                )
            
            % Initialize parameters
            obj.PHY.SamplingRate                    = SamplingRate;
            obj.Nr.SamplesTotal                     = SamplesTotal;
            obj.PHY.dt = 1/obj.PHY.SamplingRate;

            obj.NewRealization;
               
        end
        
        
        
        function NewRealization( obj )
            % obtain new realization of the channel

                    obj.ImpulseResponse = zeros(1,size(obj.Implementation.PowerDelayProfileNormalized,1), obj.Nr.txAntennas, obj.Nr.txAntennas  );
                    obj.ImpulseResponse(1,:) = (1/sqrt(2)*sqrt(obj.Implementation.PowerDelayProfileNormalized).*(randn(size(obj.Implementation.PowerDelayProfileNormalized))+1j*randn(size(obj.Implementation.PowerDelayProfileNormalized)))).';
                             
        end 
        
        
        
        
        
        
        
        function convolvedSignal = Convolution(obj, signal)
            % Convolution of the signal with the time variant impulse response
            
            N = size(signal,1);
                convolvedSignal = zeros( N);
                ConvolutionMatrix   = obj.GetConvolutionMatrix;      
                 convolvedSignal(:) = convolvedSignal(:) + ConvolutionMatrix(1:N,1:N)*signal(:);                  
        end
        
        
        function ConvolutionMatrix = GetConvolutionMatrix(obj)
            % returns the time-variant convolution matrix
            
                 ImpulseResponseTemp = obj.ImpulseResponse(:,obj.Implementation.IndexDelayTaps);
                 ConvolutionMatrix = sparse(obj.Implementation.MappingConvolutionMatrixFast(:,2),obj.Implementation.MappingConvolutionMatrixFast(:,1),ImpulseResponseTemp(obj.Implementation.CancelElementsConvolutionMatrixFast),obj.Nr.SamplesTotal,obj.Nr.SamplesTotal);

        end
        
        
        
        
        
        function ChannelTransferFunction = GetTransferFunction(obj, TimePos, FFTSize, ActiveSubcarrier)
            % Calculates the channel transfer function at certain time
            % positions. The first argument is a vector and represent the time
            % positions. The second argument represents the FFT size of the
            % underlying modulation format. The third argument (if
            % specified) determines which subcarriers are active (only 
            % these frequency positions are returned).


            ChannelTransferFunction = zeros( FFTSize, length(TimePos));
            if obj.PHY.MaximumDopplerShift==0
                TimePos = ones(length(TimePos),1);
            end

                    ImpulseTemp = [obj.ImpulseResponse(TimePos).';zeros(FFTSize-size(obj.ImpulseResponse,2),length(TimePos))];
                    ChannelTransferFunction(:,:) = fft(ImpulseTemp,[],1);
            if exist('ActiveSubcarrier','var')
                ChannelTransferFunction = ChannelTransferFunction(ActiveSubcarrier,:,:,:);
            end
        end     
        
        

        
    end
end
