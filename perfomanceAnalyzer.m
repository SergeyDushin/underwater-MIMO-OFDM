classdef perfomanceAnalyzer
    %PERFOMANCEANALYZER Summary of this class goes here
    %   Detailed explanation goes here
    
    % It is not implemented yet
    
    properties
       AnalyzerSettings
    end
    
    methods
        function obj = perfomanceAnalyzer(varargin)
            %PERFOMANCEANALYZER Construct an instance of this class
            obj.AnalyzerSettings.PaMode = varargin{1};

        end
        
        function AnalyzeData(obj,input_data)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
           plot(input_data);
        end
    end
end

