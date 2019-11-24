%% Initialize
clear
T0 = 0; % sec
AU = 1.495978707E11; % m
Sun_MU = 1.32712440018E20; % m^3/
Earth = Body(AU,Sun_MU,0);
Mars  = Body(1.5237*AU,Sun_MU,44.3453);
ShuttleTest = Shuttle(Earth,Mars,0);

UT = linspace(0,ShuttleTest.OrbitalPeriod/2);

% Graph settings
h = figure(1);
gifName = 'Test.gif';
clf
xlim([-2*AU, 2*AU])
ylim([-2*AU, 2*AU])
axis equal
grid on
legend on

SunPath   = makeBodyPath('yellow','Sun');
addpoints(SunPath,0,0);
EarthPath = makeBodyPath('blue','Earth');
MarsPath = makeBodyPath('red','Mars');

ShuttleTestPath = makeShuttlePath('green','Shuttle');

%% Sim Run
for i=1:length(UT)
    EarthPosition = Earth.PositionCalc(UT(i));
    MarsPosition = Mars.PositionCalc(UT(i));
    ShuttleTestPosition = ShuttleTest.PositionCalc(UT(i));
    
    addpoints(EarthPath,EarthPosition(1),EarthPosition(2));
    addpoints(MarsPath,MarsPosition(1),MarsPosition(2));
    addpoints(ShuttleTestPath,ShuttleTestPosition(1),ShuttleTestPosition(2));
    
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,gifName,'gif','DelayTime',0.1,'Loopcount',inf);
    else
        imwrite(imind,cm,gifName,'gif','DelayTime',0.1,'WriteMode','append');
    end
end

%% Functions
function BodyPath = makeBodyPath(Color,Name)

BodyPath = animatedline();
BodyPath.Color = Color;
BodyPath.DisplayName = Name;
BodyPath.MaximumNumPoints = 1;
BodyPath.Marker = 'o';
BodyPath.MarkerSize = 10;
BodyPath.MarkerFaceColor = Color;

end

function ShuttlePath = makeShuttlePath(Color,Name)

ShuttlePath = animatedline();
ShuttlePath.Color = Color;
ShuttlePath.DisplayName = Name;
ShuttlePath.MaximumNumPoints = 10;

end
