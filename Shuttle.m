classdef Shuttle
    %SHUTTLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SMA                     % m
        ECC                     % non-dimensional
        MeanAnomalyAtEpoch  = 0;% deg
        parentMU                % m^3/s^2
        MeanAngMotion           % rad/s
        OrbitalPeriod           % sec
        Apoapsis                % m
        Periapsis               % m
    end
    
    methods
        function obj = Shuttle(OriginBody,TargetBody,UT)
            %SHUTTLE Construct an instance of this class
            %   Detailed explanation goes here
            obj.SMA           = (OriginBody.SMA + TargetBody.SMA)/2;
            obj.ECC           = abs(OriginBody.SMA - TargetBody.SMA)/(OriginBody.SMA + TargetBody.SMA);
            obj.parentMU      = OriginBody.parentMU;
            obj.MeanAngMotion = sqrt(obj.parentMU/(obj.SMA^3));
            obj.MeanAnomalyAtEpoch = 0 - obj.MeanAngMotion*UT*180/pi;
            obj.OrbitalPeriod = 2*pi*sqrt((obj.SMA)^3/obj.parentMU);
            
        end
        function EccAnomaly = BisectionSolver_EccAno(obj,MA)
            EA = pi;
            Min = 0;
            Max = 2*pi;
            tol = 0.1*pi/180;
            
            for i=1:10
                MA_test = EA - obj.ECC*sin(EA);
                result = MA - MA_test;
                if abs(result) < tol
                    break
                end
                
                if result > 0
                    Min = EA;
                    EA = (Max + Min)/2;
                else
                    Max = EA;
                    EA = (Max + Min)/2;
                end
            end
            
            EccAnomaly = EA;
                
        end
        
        function TrueAnomaly = ThetaCalc(obj,UT)
            MeanAnomaly_rad = obj.MeanAnomalyAtEpoch*pi/180 + obj.MeanAngMotion*UT;
            EccAnomaly = obj.BisectionSolver_EccAno(MeanAnomaly_rad);
            A = sqrt((1 + obj.ECC)/(1 - obj.ECC));
            TrueAnomaly = (2*atan(A*tan(EccAnomaly/2)))*180/pi;
        end
        
        function Position = PositionCalc(obj,UT)
            TA = obj.ThetaCalc(UT);
            A = obj.SMA*(1 - obj.ECC^2);
            B = 1 + obj.ECC*cosd(TA);
            R = A/B;
            X = R*cosd(TA);
            Y = R*sind(TA);
            Position = [X;Y];
        end
    end
end

