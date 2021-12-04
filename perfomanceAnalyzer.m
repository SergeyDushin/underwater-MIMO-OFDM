classdef perfomanceAnalyzer
    %PERFOMANC EANALYZER for analyse the result
    % Settings description are here
    % Target description is here
    % Output is here
    % =====================================================================    
    % Sergey Dushin, nihsuds@gmail.com
    % (c) 2021 by V.A. Trapeznikov Institute of control sciences of RAS
    % www.ipu.ru
    %
    % =====================================================================
    % MIT license
    % =====================================================================
    % version 0.1 (basic version, just plotting In and Out Signals)
    %======================================================================
    
    properties
       AnalyzerSettings
    end
    
    methods
        function obj = perfomanceAnalyzer(varargin)
            %PERFOMANCEANALYZER Construct an instance of this class
            obj.AnalyzerSettings.PaMode = varargin{1};

        end
        
        function AnalyzeData(obj,LinkInputSignal,LinkOutputSignal,fs,SimTime)
            % choose the set of analysis instruments and call corresponding function 
            if strcmp(obj.AnalyzerSettings.PaMode, 'DataGeneratorOut')
                obj.AnalyzeGeneratorOut(LinkInputSignal, LinkOutputSignal,fs,SimTime);
            elseif strcmp(obj.AnalyzerSettings.PaMode, 'TransmitterOut')
                obj.AnalyzeTransmitterOut(LinkInputSignal, LinkOutputSignal,fs,SimTime);
            elseif strcmp(obj.AnalyzerSettings.PaMode, 'ReceiverOut')
                obj.AnalyzeReceiverOut(LinkInputSignal, LinkOutputSignal,fs,SimTime);
            else 
                error('unknown analyzer target');
            end
        end
    
        function AnalyzeGeneratorOut(obj,LinkInData, LinkOutData, fs,SimTime)
            figure(1);
            n=length(LinkOutData);
            t=(0:n-1)*(SimTime/n);
            plot(t,LinkOutData);
            xlabel('sample number');
            ylabel('Generator Output');
            
            
            spec=abs(fft(LinkOutData));
            f = (0:n-1)*(fs/n);
            
            figure(2);
            plot(f,spec);
            xlabel('Frequency');
            ylabel('Power');
            
        end

        function AnalyzeTransmitterOut(obj,LinkInData, LinkOutData,fs,SimTime)
            
            %TBD
        end

        function AnalyzeReceiverOut(obj,LinkInData, LinkOutData,fs,SimTime)
            
            %TBD
            
            % here is BER estimaton
            % here is SIR estimation
        end    
    
    end
end
