classdef dataGenerator
    % DATA GENERATOR 
    % Settings description are here
    % Output is here
    % =====================================================================    
    % Sergey Dushin, s.dushin@sdlab-co.ru
    % (c) 2021 by SDLab-co
    % www.sdlab-co.ru
    %
    % =====================================================================
    % MIT license
    % =====================================================================
    %version 0.1
    %======================================================================
    
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
        
        function [Data, fs] = generateData(obj)
            % Generate the data using defined settings
            Len=obj.GeneratorSettings.GenerationTime*obj.GeneratorSettings.DataSampleRate;
            fs=obj.GeneratorSettings.DataSampleRate;
            switch obj.GeneratorSettings.waveform
                case 'wgn'
                    Data = wgn(1,Len,0);
                case 'sin'
                    %Data = sin(1,Len,0);
                otherwise
                    error('dataGenerator: Defined waveform is anavalable');
            end
        end
    end
end

