classdef channel
    %CHANNEL Summary of this class goes here
    %   Detailed explanation goes here
    
    % It is not implemented yet
    
    properties
        Property1
    end
    
    methods
        function obj = channel(varargin)
            %CHANNEL Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = varargin{1};
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

