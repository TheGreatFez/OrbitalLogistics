classdef Body
    %Assumed to be in a circular orbit
    %   Detailed explanation goes here
    
    properties
        SMA
        MeanAnomalyAtEpoch = 0; % In Universal Time (UT), time
        parentMU                % m^3/s^2
        MeanAngMotion
        OrbitalSpeed
        OrbitalPeriod
    end
    
    methods
        function obj = Body(SMA,parentMU,MA)
            %BODY Construct an instance of this class
            %   SMA is the Semi Major Axis of the Body's orbit
            %   parentMU is the Gravitational Parameter (MU) of the parent body
            %   MA is the Mean Anomaly at Epoch of the Body
            obj.SMA           = SMA;
            obj.parentMU      = parentMU;
            obj.MeanAngMotion = sqrt(obj.parentMU/(obj.SMA^3));
            obj.MeanAnomalyAtEpoch = MA;
            obj.OrbitalSpeed  = sqrt(obj.parentMU/obj.SMA);
            obj.OrbitalPeriod = 2*pi*sqrt((obj.SMA)^3/obj.parentMU);
        end
        
        function TrueAnomaly = ThetaCalc(obj,UT)
            %ThetaCalc: Calculate the True Anomaly at a given UT
            TrueAnomaly = (UT)*obj.MeanAngMotion*180/pi + obj.MeanAnomalyAtEpoch;
        end
        
        function Position = PositionCalc(obj,UT)
            %PositionCalc: Calculate the Position at a given UT
            X = obj.SMA*cosd(obj.ThetaCalc(UT));
            Y = obj.SMA*sind(obj.ThetaCalc(UT));
            Position = [X;Y];
        end
        
        
    end
end

