classdef dataGenerator
    % DATA GENERATOR 
    %   Class for generation of data to be transmitted
    
    % version 0.1
    
    properties (SetAccess = private)
        GeneratorSettings
    end
    
    methods
        function obj = dataGenerator(varargin)
            % DATA GENERATOR constructor
            obj.GeneratorSettings.GenerationTime = varargin{1}; % in seconds
            obj.GeneratorSettings.DataSampleRate = varargin{2}; % samples per second
            obj.GeneratorSettings.waveform = varargin{3};
        end
        
        function Data = generateData(obj)
            % Generate the data using defined settings
            Len=obj.GeneratorSettings.GenerationTime*obj.GeneratorSettings.DataSampleRate;
            switch obj.GeneratorSettings.waveform
                case 'wgn'
                    Data = wgn(1,Len,0);
                case 'sin'
                    Data = sin(1,Len,0);
                otherwise
                    error('dataGenerator: Defined waveform is anavalable');
            end
        end
    end
end
