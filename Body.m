classdef Body
    %Assumed to be in a circular orbit
    %   Detailed explanation goes here
    
    properties
        SMA
        MeanAnomalyAtEpoch = 0; % In Universal Time (UT), time
        parentMU                % m^3/s^2
        MeanAngMotion           % deg/s
        OrbitalSpeed            % m/s
        OrbitalPeriod           % sec
        Radius                  % Radius of the body, m
        SoI                     % Sphere of Influence, m
        Mass                    % kg
    end
    
    methods
        function obj = Body(SMA,parentMU,MA,Radius,Mass)
            %BODY Construct an instance of this class
            %   SMA is the Semi Major Axis of the Body's orbit
            %   parentMU is the Gravitational Parameter (MU) of the parent body
            %   MA is the Mean Anomaly at Epoch of the Body
            G                 = 6.67408E-11;
            obj.SMA           = SMA;
            obj.parentMU      = parentMU;
            obj.MeanAngMotion = sqrt(obj.parentMU/(obj.SMA^3));
            obj.MeanAnomalyAtEpoch = MA;
            obj.OrbitalSpeed  = sqrt(obj.parentMU/obj.SMA);
            obj.OrbitalPeriod = 2*pi*sqrt((obj.SMA)^3/obj.parentMU);
            obj.Radius        = Radius;
            obj.Mass          = Mass;
            obj.SoI           = obj.SMA*(obj.Mass/(obj.parentMU/G))^(2/5);
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

