classdef Shuttle
    %SHUTTLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SMA                     % m
        ECC                     % non-dimensional
        MeanAnomalyAtEpoch      % deg
        ArgumentOfPeriapsis     % deg
        parentMU                % m^3/s^2
        MeanAngMotion           % rad/s
        OrbitalPeriod           % sec
        Apoapsis                % m
        Periapsis               % m
        DeltaV                  % m/s
    end
    
    methods
        function obj = Shuttle(OriginBody,TargetBody,UT)
            %SHUTTLE Construct an instance of this class
            %   Detailed explanation goes here
            obj.parentMU      = OriginBody.parentMU;
            
            % Determine the Intercept Orbit
            thetaOB    = OriginBody.ThetaCalc(UT);
            thetaTB    = TargetBody.ThetaCalc(UT);
            thetaStart = thetaOB;
            
            currentSMA = (OriginBody.SMA + TargetBody.SMA)/2;
            currentEcc = abs(OriginBody.SMA - TargetBody.SMA)/(OriginBody.SMA + TargetBody.SMA);
            currentMeanMotion = sqrt(obj.parentMU/currentSMA^3);
            
            % Calculate current Orbital Speed assuming hohmann transfer and
            % determine Min and Max speeds for solver
            currentSpeedStart = obj.VisVia(currentSMA,OriginBody.SMA);
            if OriginBody.SMA > TargetBody.SMA
                MaxSpeed = currentSpeedStart;
                MinSpeed = 0;
                obj.ArgumentOfPeriapsis = thetaOB + 180;
                if obj.ArgumentOfPeriapsis > 360
                    obj.ArgumentOfPeriapsis = obj.ArgumentOfPeriapsis - 360;
                end
            else
                MaxSpeed = sqrt(2*obj.parentMU/OriginBody.SMA) - 1;
                MinSpeed = currentSpeedStart;
                obj.ArgumentOfPeriapsis = thetaOB;
            end
            InstanceSpeed1 = NaN;
            InstanceSpeed2 = NaN;
            % Solver for intercept orbit
            for Instance = 1:2 % First or Second instance of hitting the Body's Orbit
                currentSpeed = currentSpeedStart;
                for i = 1:10%50
                    % Determine where Shuttle will intercept the orbit of the
                    % target body
                    specificOE = currentSpeed^2/2 - obj.parentMU/OriginBody.SMA;
                    currentSMA = -obj.parentMU/(2*specificOE);
                    currentMeanMotion = sqrt(obj.parentMU/(currentSMA^3));
                    if OriginBody.SMA > TargetBody.SMA
                        currentApoapsis  = OriginBody.SMA;
                        currentPeriapsis = 2*currentSMA - currentApoapsis;
                        currentEcc       = (currentApoapsis - currentPeriapsis)/(currentApoapsis + currentPeriapsis);
                    else
                        currentPeriapsis = OriginBody.SMA;
                        currentApoapsis  = 2*currentSMA - currentPeriapsis;
                        currentEcc       = (currentApoapsis - currentPeriapsis)/(currentApoapsis + currentPeriapsis);
                    end
                    thetaIntercept = acosd(max(-1,min(1,((currentSMA*(1-currentEcc^2)/TargetBody.SMA)-1)/currentEcc)));
                    
                    if Instance == 2
                        thetaIntercept = 360 - thetaIntercept;
                    end
                    % Determine when Shuttle will intercept the orbit of the
                    % target body
                    EccAnomaly = 2*atan(sqrt((1-currentEcc)/(1+currentEcc))*tand(thetaIntercept/2));
                    MeanAnomaly = EccAnomaly - currentEcc*sin(EccAnomaly);
                    timeDelta = MeanAnomaly/currentMeanMotion;
                    interceptUT = timeDelta + UT;
                    
                    % Determine distance b/w Shuttle and Target Body when they
                    % intercept
                    targetThetaCheck = TargetBody.ThetaCalc(interceptUT);
                    shuttleThetaCheck = thetaIntercept + obj.ArgumentOfPeriapsis;
                    if shuttleThetaCheck > 360
                        shuttleThetaCheck = shuttleThetaCheck - 360;
                    end
                    
                    % Determine if the Shuttle has intercepted the Target Body
                    diffThetaCheck = abs(targetThetaCheck - shuttleThetaCheck);
                    if diffThetaCheck > 180
                        diffThetaCheck = 360 - diffThetaCheck;
                    end
                    distThetaCheck = 2*TargetBody.SMA*sind(diffThetaCheck/2);
                    
                    if distThetaCheck < TargetBody.SoI
                        if Instance == 1
                            InstanceSpeed1 = currentSpeed;
                        else
                            InstanceSpeed2 = currentSpeed;
                        end
                        break;
                    end
                    if Instance == 1
                        if shuttleThetaCheck < targetThetaCheck
                            MinSpeed = currentSpeed;
                            currentSpeed = (MinSpeed + MaxSpeed)/2;
                        else
                            if currentSpeed == currentSpeedStart
                                break;
                            else
                                MaxSpeed = currentSpeed;
                                currentSpeed = (MinSpeed + MaxSpeed)/2;
                            end
                        end
                    else
                        if shuttleThetaCheck > targetThetaCheck
                            MinSpeed = currentSpeed;
                            currentSpeed = (MinSpeed + MaxSpeed)/2;
                        else
                            if currentSpeed == currentSpeedStart
                                break;
                            else
                                MaxSpeed = currentSpeed;
                                currentSpeed = (MinSpeed + MaxSpeed)/2;
                            end
                        end
                    end
                end
                finalSpeed = min([InstanceSpeed1,InstanceSpeed2]);
            end
            
            specificOE = finalSpeed^2/2 - obj.parentMU/OriginBody.SMA;
            obj.SMA           = -obj.parentMU/(2*specificOE);
            obj.MeanAngMotion = sqrt(obj.parentMU/(obj.SMA^3));
            if OriginBody.SMA > TargetBody.SMA
                obj.Apoapsis  = OriginBody.SMA;
                obj.Periapsis = 2*obj.SMA - obj.Apoapsis;
                obj.ECC       = (obj.Apoapsis - obj.Periapsis)/(obj.Apoapsis + obj.Periapsis);
                obj.MeanAnomalyAtEpoch = 180 - obj.MeanAngMotion*UT*180/pi;
            else
                obj.Periapsis = OriginBody.SMA;
                obj.Apoapsis  = 2*obj.SMA - obj.Periapsis;
                obj.ECC       = (obj.Apoapsis - obj.Periapsis)/(obj.Apoapsis + obj.Periapsis);
                obj.MeanAnomalyAtEpoch = 0 - obj.MeanAngMotion*UT*180/pi;
            end
            obj.OrbitalPeriod = 2*pi*sqrt((obj.SMA)^3/obj.parentMU);
        end
        
        function OrbitalSpeed = VisVia(obj,SMA,R)
            OrbitalSpeed = sqrt(obj.parentMU*(2/R - 1/SMA));
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
            posAngle = obj.ThetaCalc(UT) + obj.ArgumentOfPeriapsis;
            A = obj.SMA*(1 - obj.ECC^2);
            B = 1 + obj.ECC*cosd(posAngle);
            R = A/B;
            X = R*cosd(posAngle);
            Y = R*sind(posAngle);
            Position = [X;Y];
        end
    end
end

