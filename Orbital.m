classdef Orbital
    %ORBITAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        SMA
        MeanAngMotion
        Period
        MAatEpoch = 0;  % Degrees
        Ecc = 0;
        Inc = 0;        % Degrees
        AoP = 0;        % Degrees
        LAN = 0;        % Degrees
        parentMU = 1.32712440018E20; % m^3/s^2
    end
    
    methods
        function obj = Orbital(Name,SMA,OrbitParams)
            %ORBITAL Construct an instance of this class
            %   OrbitParams is a structure with 
            obj.Name = Name;
            obj.SMA  = SMA;
            obj.MeanAngMotion = sqrt(obj.parentMU/obj.SMA^3);
            obj.Period = 2*pi*sqrt(obj.SMA^3/obj.parentMU);
            if exist('OrbitParams','var') 
                Fields = fieldnames(OrbitParams);
                for i=1:length(Fields)
                    obj.(Fields{i}) = OrbitParams.(Fields{i});
                end
            end
        end
        
        function OrbitalSpeed = VisVia(obj,SMA,R)
            OrbitalSpeed = sqrt(obj.parentMU*(2/R - 1/SMA));
        end
        
        function EccAnomaly = EA_Solver(obj,MA)
            EA = pi;
            Min = 0;
            Max = 2*pi;
            tol = 0.001*pi/180;
            
            for i=1:30
                MA_test = EA - obj.Ecc*sin(EA);
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
            disp(['Iterations ' num2str(i)])
            
            EccAnomaly = EA;
            
        end
        
        function TrueAnomaly = ThetaCalc(obj,UT)
            MeanAnomaly_rad = obj.MAatEpoch*pi/180 + obj.MeanAngMotion*UT;
            EccAnomaly = obj.EA_Solver(MeanAnomaly_rad);
            A = sqrt((1 + obj.Ecc)/(1 - obj.Ecc));
            TrueAnomaly = (2*atan(A*tan(EccAnomaly/2)))*180/pi;
        end
        
        function Position = PositionCalc(obj,UT)
            posAngle = obj.ThetaCalc(UT) + obj.AoP;
            A = obj.SMA*(1 - obj.Ecc^2);
            B = 1 + obj.Ecc*cosd(posAngle);
            R = A/B;
            X = R*cosd(posAngle);
            Y = R*sind(posAngle);
            Position = [X;Y;0];
        end
    end
end

