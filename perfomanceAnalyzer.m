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
    % version 0.1
    %======================================================================
    
    properties
       AnalyzerSettings
    end
    
    methods
        function obj = perfomanceAnalyzer(varargin)
            %PERFOMANCEANALYZER Construct an instance of this class
            obj.AnalyzerSettings.PaMode = varargin{1};

        end
        
        function AnalyzeData(obj,InputData)
            % choose the set of analysis instruments and call corresponding function 
            if strcmp(obj.AnalyzerSettings.PaMode, 'DataGeneratorOut')
                obj.AnalyzeGeneratorOut(InputData);
            elseif strcmp(obj.AnalyzerSettings.PaMode, 'TransmitterOut')
                obj.AnalyzeTransmitterOut(InputData);
            elseif strcmp(obj.AnalyzerSettings.PaMode, 'ReceiverOut')
                obj.AnalyzeReceiverOut(InputData);
            else 
                error('unknown analyzer target');
            end
        end
    
        function AnalyzeGeneratorOut(obj,data)
            plot(data);
        end

        function AnalyzeTransmitterOut(obj,data)
            plot(data);
        end

        function AnalyzeReceiverOut(obj,data)
            plot(data);
        end    
    
    end
end
