%% Final Project AAE 568 MATLAB Simulation Code

clc; 
clear all; 
close all;
run("createTrajectory.m");

% Conversion constants
d2r = pi/180;      % degrees to radians
r2d = 180/pi;      % radians to degrees
gToMps = 9.80665;  % Gs to Meters per Sec
toMilli = 1000;    % meters to millimeters
fromMilli = 0.001; % millimeters to meters 
ft2m = 0.305;      % feet to meters
hr2sec = 3600;     % hours to seconds
gravity = 9.80665; % meters/sec/sec gravity acceleration constant

%% Problem Setup

% Configuration 1
timeStep   = 0.05; % Seconds
updateTime = 0;

% % Configuration 2
% timeStep   = 0.1; % Seconds
% updateTime = 0;
% endSim     = false; % Termination flag
% satUpdateFreq = 2; % Sec
% imuUpdateFreq = 4; % Sec
% satUpdateTime = 0.5; % Next time to update
% imuUpdateTime = 1; % Next time to update

%% Simulation Initialization Variables

%
% Flags, timing, indices
%
endSim  = false; % Termination flag
i       = 1    ; % Loop iteration number
simTime = 0    ; % Current simulation time 

%
% Kalman Characteristics TODO: Should fix these values 
%
q    = 3; % Process uncertainty
Q    = q^2*eye(12);
ripH = 1   ; % IMU Measurement uncertainty (horizontal)
ripV = 1   ; % IMU Measurement uncertainty (vertcal)
ripE = .158; % IMU Measurement uncertainty (euler angles)
ripR = 1.58; % IMU Measurement uncertainty (euler rates)
rgpH = 0.5 ; % GPS Measurement uncertainty (horizontal)
rgpV = 0.3 ; % GPS measurement uncertainty (vertical)
rgv  = 0.55; % GPS Velocity uncertainty (all directions)

% Covariance of the IMU (x,y,z,u,v,w)
Rimu = [ ripH^2,      0,      0,     0,     0,     0,      0,      0,      0,      0,      0,      0;
              0, ripH^2,      0,     0,     0,     0,      0,      0,      0,      0,      0,      0;
              0,      0, ripV^2,     0,     0,     0,      0,      0,      0,      0,      0,      0;
              0,      0,      0, 0.1^2,     0,     0,      0,      0,      0,      0,      0,      0;
              0,      0,      0,     0, 0.1^2,     0,      0,      0,      0,      0,      0,      0;
              0,      0,      0,     0,     0, 0.1^2,      0,      0,      0,      0,      0,      0;
              0,      0,      0,     0,     0,     0, ripE^2,      0,      0,      0,      0,      0;
              0,      0,      0,     0,     0,     0,      0, ripE^2,      0,      0,      0,      0;
              0,      0,      0,     0,     0,     0,      0,      0, ripE^2,      0,      0,      0;
              0,      0,      0,     0,     0,     0,      0,      0,      0, ripR^2,      0,      0;
              0,      0,      0,     0,     0,     0,      0,      0,      0,      0, ripR^2,      0;
              0,      0,      0,     0,     0,     0,      0,      0,      0,      0,      0, ripR^2 ];

% Covariance of the GPS (x,y,z,u,v,w)
Rgps = [ rgpH^2,      0,      0,     0,     0,     0, 0, 0, 0, 0, 0, 0;
              0, rgpH^2,      0,     0,     0,     0, 0, 0, 0, 0, 0, 0;
              0,      0, rgpV^2,     0,     0,     0, 0, 0, 0, 0, 0, 0;
              0,      0,      0, rgv^2,     0,     0, 0, 0, 0, 0, 0, 0;
              0,      0,      0,     0, rgv^2,     0, 0, 0, 0, 0, 0, 0;
              0,      0,      0,     0,     0, rgv^2, 0, 0, 0, 0, 0, 0;
              0,      0,      0,     0,     0,     0, 1, 0, 0, 0, 0, 0;
              0,      0,      0,     0,     0,     0, 0, 1, 0, 0, 0, 0;
              0,      0,      0,     0,     0,     0, 0, 0, 1, 0, 0, 0;
              0,      0,      0,     0,     0,     0, 0, 0, 0, 1, 0, 0;
              0,      0,      0,     0,     0,     0, 0, 0, 0, 0, 1, 0;
              0,      0,      0,     0,     0,     0, 0, 0, 0, 0, 0, 1];

% Total initial covariance
P = Rimu.*Rgps;

%% Initial Conditions

%
% Reference frames
%
% Geodetic Frame Position (LLA) [user defined / from running "createTrajectory.m"]
geoPosLat(i) = initialPositionLLA(1); % degrees
geoPosLon(i) = initialPositionLLA(2); % degrees
geoPosAlt(i) = initialPositionLLA(3); % meters
geoPosLLA(:,i) = [geoPosLat; geoPosLon; geoPosAlt];

% ECEF Frame Position [rotated from Geodetic]
[ecefPosX(i), ecefPosY(i), ecefPosZ(i)] = geodetic2ecef(wgs84Ellipsoid, geoPosLat, geoPosLon, geoPosAlt);
ecefPos(:,i) = [ecefPosX; ecefPosY; ecefPosZ]; % meters

% NED Frame Position [rotated from ECEF with reference at takeoff point 
% i.e. Geodetic initial LLA]
[nedPosN(i), nedPosE(i), nedPosD(i)] = ecef2ned(ecefPosX, ecefPosY, ecefPosZ, geoPosLat(1), geoPosLon(1), geoPosAlt(1), wgs84Ellipsoid);
nedPos(:,i) = [nedPosN; nedPosE; nedPosD]; % meters

% Body Frame Components
% User defined
bodyRoll(i)  = 0.0; % deg
bodyPitch(i) = 0.0; % deg
bodyYaw(i)   = 0.0; % deg
bodyEulerAng(:,i) = [bodyRoll; bodyPitch; bodyYaw];

bodyRollRate(i)  = 0.0; % deg/sec
bodyPitchRate(i) = 0.0; % deg/sec
bodyYawRate(i)   = 0.0; % deg/sec
bodyEulerAngRate(:,i) = [bodyRollRate; bodyPitchRate; bodyYawRate];

bodyVelX(i) = spec.speedMetersPerSec; % meters/second
bodyVelY(i) = 0.0; % meters/second
bodyVelZ(i) = 0.0; % meters/second
bodyVel(:,i) = [bodyVelX; bodyVelY; bodyVelZ];

bodyAccX(i) = 0.0;     % meters/second/second
bodyAccY(i) = 0.0;     % meters/second/second
bodyAccZ(i) = gravity; % meters/second/second
bodyAcc(:,i) = [bodyAccX; bodyAccY; bodyAccZ];

%% Initial dynamics 

wpTo   = 2; % Waypoint index you're traveling to
wpFrom = 1; % Waypoint index you're traveling from

% Generic heading command. The vector points to the next waypoint from
% where you currently are (at this point you're at the first waypoint)
unitDirHeadingNED = transpose((waypointsNED(wpTo,:) - waypointsNED(wpFrom,:)) ./ norm(waypointsNED(wpTo,:) - waypointsNED(wpFrom,:)));

% Velocity in NED 
nedVel(:,i) = unitDirHeadingNED .* spec.speedMetersPerSec;

% Kalman filter outputs and track for performance analysis
filterEst(:,i) = [nedPos(:,i); nedVel(:,i); bodyEulerAng(:,i); bodyEulerAngRate(:,i)];
filterCovPosX(:,i)      = P(1,1);
filterCovPosY(:,i)      = P(2,2);
filterCovPosZ(:,i)      = P(3,3);
filterCovVelX(:,i)      = P(4,4);
filterCovVelY(:,i)      = P(5,5);
filterCovVelZ(:,i)      = P(6,6);
filterCovRoll(:,i)      = P(7,7);
filterCovPitch(:,i)     = P(8,8);
filterCovYaw(:,i)       = P(9,9);
filterCovRollRate(:,i)  = P(7,7);
filterCovPitchRate(:,i) = P(8,8);
filterCovYawRate(:,i)   = P(9,9);
x = filterEst(:,i);

%% Create IMU and GPS sensors

% GPS sensor properties
satUpdateRate = 1;   % Sec
satUpdateTime = 0; % Next time to update
gps = gpsSensor('SampleRate', 1/satUpdateRate,'ReferenceLocation',geoPosLLA');

% IMU sensor properties
imuUpdateRate = 0.1; % Sec
imuUpdateTime = 0.1; % Next time to update
paramsAccel = accelparams('MeasurementRange',19.6,'Resolution',0.598e-3,'ConstantBias',0.49,'NoiseDensity',3920e-6,'TemperatureBias',0.294,'TemperatureScaleFactor',0.02,'AxesMisalignment',2);
paramsGyro = gyroparams('MeasurementRange',6.2832,'Resolution',1.332e-04,'ConstantBias',0.0349,'NoiseDensity',8.727e-4,'TemperatureBias',0.349,'TemperatureScaleFactor',0.02,'AxesMisalignment',2,'AccelerationBias',0.178e-3);
imu = imuSensor(...
     'accel-gyro', ...
     'Accelerometer', paramsAccel, ...
     'Gyroscope', paramsGyro, ...
     'SampleRate', 1/imuUpdateRate);
orient = quaternion.ones(1, 1);


%% Create true flight path taken
while ((endSim == false) && (i < 100000)) % 100,000 is to ensure while-exit criteria

    % Update time step and iteration
    simTime = simTime + timeStep;
    i = i + 1;

    % Simulate measurement from the sensor - calculate velocity and
    % position readings in the body and rotate to NED
    if ( simTime >= imuUpdateTime )

        % Time delta
        dt = min(simTime-updateTime, timeStep);

        %
        % Process IMU measurement
        %

        % Accelerometer produces body accel, gyro produces body Euler rates
        [accelImuMeas, eulerRateImuMeas] = imu( bodyAcc(:,i-1)', bodyEulerAngRate(:,i-1)'.*d2r, orient );
        % Euler rate to degrees, integrate to get Euler angles
        eulerRateImuMeas = eulerRateImuMeas .* r2d;
        eulerAngImuMeas = bodyEulerAng(:,i-1) + eulerRateImuMeas' .* timeStep;

        % Integrate body acceleration once to get body velocity, then
        % rotate to NED and integrate again to get new NED position
        bodyVelImuMeas = bodyVel(:,i-1) + accelImuMeas' .* timeStep;
        [nedVelImuMeas(:,1), nedVelImuMeas(:,2), nedVelImuMeas(:,3)] = body2ned(bodyVelImuMeas(1), bodyVelImuMeas(2), bodyVelImuMeas(3), bodyEulerAng(1,i-1), bodyEulerAng(2,i-1), bodyEulerAng(3,i-1));
        nedPosImuMeas = nedPos(:,i-1) + nedVelImuMeas' .* timeStep;

        %
        % Kalman update
        %
        
        % Measurement vector (12x1)
        y = [nedPosImuMeas; nedVelImuMeas'; eulerAngImuMeas; eulerRateImuMeas'];

        u = [0; 0; 0; nedVel(:,i-1); 0; 0; 0; bodyEulerAngRate(:,i-1)];

        F = [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ];

        B = [ 0, 0, 0,  0,  0,  0, 0, 0, 0,  0,  0,  0;
              0, 0, 0,  0,  0,  0, 0, 0, 0,  0,  0,  0;
              0, 0, 0,  0,  0,  0, 0, 0, 0,  0,  0,  0;
              0, 0, 0, dt,  0,  0, 0, 0, 0,  0,  0,  0;
              0, 0, 0,  0, dt,  0, 0, 0, 0,  0,  0,  0;
              0, 0, 0,  0,  0, dt, 0, 0, 0,  0,  0,  0;
              0, 0, 0,  0,  0,  0, 0, 0, 0,  0,  0,  0;
              0, 0, 0,  0,  0,  0, 0, 0, 0,  0,  0,  0;
              0, 0, 0,  0,  0,  0, 0, 0, 0,  0,  0,  0;
              0, 0, 0,  0,  0,  0, 0, 0, 0, dt,  0,  0;
              0, 0, 0,  0,  0,  0, 0, 0, 0,  0, dt,  0;
              0, 0, 0,  0,  0,  0, 0, 0, 0,  0,  0, dt ]; 

        H = [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ];

        [x, P] = kalmanFilter(F, x, B, u, P, H, y, Q, Rimu);
        updateTime = simTime;
        imuUpdateTime = imuUpdateTime + imuUpdateRate;
    end

    % Simulate measurement from the GPS signal
    % add measurments of X and Y with some error and noise
    if ( simTime >= satUpdateTime )

        % Time delta
        dt = min(simTime-updateTime, timeStep);

        % Get GPS measurement of position in geodetic and velocity in NED
        [geoPosLLAGpsMeas, nedVelLLAGpsMeas] = gps(nedPos(:,i-1)', nedVel(:,i-1)');
        % Rotate position to ECEF
        [ecefGpsPosMeas(1), ecefGpsPosMeas(2), ecefGpsPosMeas(3)] = geodetic2ecef(wgs84Ellipsoid, geoPosLLAGpsMeas(1), geoPosLLAGpsMeas(2), geoPosLLAGpsMeas(3));
        % Rotate position to NED
        [nedGpsPosMeas(1), nedGpsPosMeas(2), nedGpsPosMeas(3)] = ecef2ned(ecefGpsPosMeas(1), ecefGpsPosMeas(2), ecefGpsPosMeas(3), geoPosLat(1), geoPosLon(1), geoPosAlt(1), wgs84Ellipsoid);

        y = [nedGpsPosMeas'; nedVelLLAGpsMeas'; 0; 0; 0; 0; 0; 0];

        u = [0; 0; 0; nedVel(:,i-1); 0; 0; 0; 0; 0; 0];

        F = [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ];

        B = [ 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, dt,  0,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0, dt,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0,  0, dt, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0;
              0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0 ];

        H = [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

        [x, P] = kalmanFilter(F, x, B, u, P, H, y, Q, Rgps);
        updateTime = simTime;
        satUpdateTime = satUpdateTime + satUpdateRate;
    end

    % If a time step goes by and no measurement was made, propagate the
    % state and the covariance
    if ( simTime > updateTime )
        x(1:3) = x(1:3)' + x(4:6)'.*timeStep;
        x(7:9) = x(7:9)' + x(10:12)'.*timeStep;
        P = P + Q;
    end

    % Log the filter state estimate and covariance matrix
    filterEst(:,i)          = x;
    filterCovPosX(:,i)      = P(1,1);
    filterCovPosY(:,i)      = P(2,2);
    filterCovPosZ(:,i)      = P(3,3);
    filterCovVelX(:,i)      = P(4,4);
    filterCovVelY(:,i)      = P(5,5);
    filterCovVelZ(:,i)      = P(6,6);
    filterCovRoll(:,i)      = P(7,7);
    filterCovPitch(:,i)     = P(8,8);
    filterCovYaw(:,i)       = P(9,9);
    filterCovRollRate(:,i)  = P(10,10);
    filterCovPitchRate(:,i) = P(11,11);
    filterCovYawRate(:,i)   = P(12,12);
    
    % Propagate the position based on the velocity from before for the next
    % step
    nedPos(:,i) = nedPos(:,i-1) + nedVel(:,i-1) .* timeStep;
    [ecefPos(1,i), ecefPos(2,i), ecefPos(3,i)] = ned2ecef(nedPos(1,i), nedPos(2,i), nedPos(3,i), geoPosLLA(1,1), geoPosLLA(2,1), geoPosLLA(3,1), wgs84Ellipsoid);
    [geoPosLLA(1,i), geoPosLLA(2,i), geoPosLLA(3,i)] = ecef2geodetic(wgs84Ellipsoid, ecefPos(1,i), ecefPos(2,i), ecefPos(3,i));

    % Determine the acceleration for the next step
    bodyAcc(:,i) = bodyAcc(:,i-1);

    % If you've passed the waypoint of interest, head towards the next one
    display( num2str(distance(geoPosLLA(1,i-1), geoPosLLA(2,i-1), waypointsLLA(wpTo,1), waypointsLLA(wpTo,2), wgs84Ellipsoid)) );
    if ( distance(geoPosLLA(1,i-1), geoPosLLA(2,i-1), waypointsLLA(wpTo,1), waypointsLLA(wpTo,2), wgs84Ellipsoid) < (spec.speedMetersPerSec * timeStep) )
        wpFrom = wpTo;
        wpTo = min(wpTo + 1, height(waypointsLLA));
    end

    % Find direction of travel for the next iteration
    unitDirHeadingNED = (waypointsNED(wpTo,:)' - filterEst(1:3,i)) ./ norm(waypointsNED(wpTo,:)' - filterEst(1:3,i));

    % Determine the velocity for the next step
    nedVel(:,i)   = unitDirHeadingNED .* spec.speedMetersPerSec;

    % Determine Euler angle necessary to achieve this 
    reqPitch = atan2d(sqrt(nedVel(1,i)^2 + nedVel(2,i)^2), nedVel(3,i));
    if ((bodyEulerAng(2,i-1) + spec.pitchRateDegPerSec * timeStep) < reqPitch)
        bodyEulerAngRate(2,i) = spec.pitchRateDegPerSec;
        bodyEulerAng(2,i) = bodyEulerAng(2,i-1) + bodyEulerAngRate(2,i) * timeStep;
    elseif ((bodyEulerAng(2,i-1) - spec.pitchRateDegPerSec * timeStep) > reqPitch)
         bodyEulerAngRate(2,i) = -spec.pitchRateDegPerSec;
         bodyEulerAng(2,i) = bodyEulerAng(2,i-1) + bodyEulerAngRate(2,i) * timeStep;
    else
        bodyEulerAngRate(2,i) = 0;
        bodyEulerAng(2,i) = reqPitch;
    end

    reqYaw = atan2d(nedVel(1,i), -nedVel(2,i));
    if ((bodyEulerAng(3,i-1) + spec.yawRateDegPerSec * timeStep) < reqYaw)
        bodyEulerAngRate(3,i) = spec.yawRateDegPerSec;
        bodyEulerAng(3,i) = bodyEulerAng(3,i-1) + bodyEulerAngRate(3,i) * timeStep;
    elseif ((bodyEulerAng(3,i-1) - spec.yawRateDegPerSec * timeStep) > reqYaw)
        bodyEulerAngRate(3,i) = -spec.yawRateDegPerSec;
        bodyEulerAng(3,i) = bodyEulerAng(3,i-1) + bodyEulerAngRate(3,i) * timeStep;
    else
        bodyEulerAngRate(3,i) = 0;
        bodyEulerAng(3,i) = reqYaw;  
    end

    % Calculate body velocity by rotating NED velocity into body
    [bodyVel(1,i), bodyVel(2,i), bodyVel(3,i)] = ned2body(nedVel(1,i), nedVel(2,i), nedVel(3,i), bodyEulerAng(1,i), bodyEulerAng(2,i), bodyEulerAng(3,i));

    % If you passed the final waypoint, terminate simulation
    if ( wpTo == height(waypointsNED) && ...
            distance(geoPosLLA(1,i-1), geoPosLLA(2,i-1), waypointsLLA(wpTo,1), waypointsLLA(wpTo,2), wgs84Ellipsoid) < (spec.speedMetersPerSec * timeStep) )
        endSim = true;
    end

end

%% DEBUG
outData = 1:100:i;
timeVec = timeStep.*[0:1:i];
timeVec = timeVec(outData);

% True Geodetic position waypoint based
figure; hold on; grid on;
plot3(waypointsLLA(:,1), waypointsLLA(:,2), waypointsLLA(:,3), 'b');
scatter3(waypointsLLA(1,1), waypointsLLA(1,2), waypointsLLA(1,3), 'ko', 'filled');
scatter3(waypointsLLA(2:end-1,1), waypointsLLA(2:end-1,2), waypointsLLA(2:end-1,3), 'bo');
scatter3(geoPosLLA(1,outData), geoPosLLA(2,outData), geoPosLLA(3,outData), 'r.');
view(3)
title('Geodetic Position Visualization');
xlabel('Latitude (deg)');
ylabel('Longitude (deg)');
zlabel('Altitude (meters)');
legend('Planned Path', 'Starting Waypoint', 'Waypoints', 'True Position', 'Location', 'Best')

%% Plotting 

% Filter the data points so that plots are more visible
outData = 1:100:i;
timeVec = timeStep.*[0:1:i];
timeVec = timeVec(outData);

% True position waypoint based
figure; hold on; grid on;
plot3(waypointsNED(:,1), waypointsNED(:,2), -waypointsNED(:,3), 'b');
scatter3(waypointsNED(1,1), waypointsNED(1,2), -waypointsNED(1,3), 'ko', 'filled');
scatter3(waypointsNED(2:end-1,1), waypointsNED(2:end-1,2), -waypointsNED(2:end-1,3), 'bo');
scatter3(nedPos(1,outData), nedPos(2,outData), -nedPos(3,outData), 'r.');
plot3(filterEst(1,outData), filterEst(2,outData), -filterEst(3,outData), 'g');
view(3)
title('NED Position Visualization');
xlabel('xNorth (meters)');
ylabel('Negative yEast (meters)');
zlabel('Negative zDown (meters)');
legend('Planned Path', 'Starting Waypoint', 'Waypoints', 'True Position', 'KF Estimate', 'Location', 'Best')

% Filtered position for tracking - 3D
figure; hold on; grid on;
scatter3(nedPos(1,outData), nedPos(2,outData), -nedPos(3,outData), 'r.');
plot3(filterEst(1,outData), filterEst(2,outData), -filterEst(3,outData), 'b')
view(3)
legend('True Body Pos', 'KF Estimation', 'Location', 'Best')
title('Body Position Filter Estimate Comparison');
xlabel('xNorth (meters)');
ylabel('yEast (meters)');
zlabel('Negative zDown (meters)');

% Filtered position for tracking - 2D
figure; hold on; grid on;
scatter(nedPos(1,outData), nedPos(2,outData), 'r.');
plot(filterEst(1,outData), filterEst(2,outData), 'b')
legend('True Body Pos', 'KF Estimation', 'Location', 'Best')
title('Body Position Filter Estimate Comparison');
xlabel('xNorth (meters)');
ylabel('yEast (meters)');

% For axes linking
axInd = 1;

% X position compared to error in NED
figure;
subplot(3,1,1:2)
hold on; grid on;
title('North Position Comparison')
plot(timeVec, nedPos(1,outData), 'b');
plot(timeVec, filterEst(1,outData), 'r');
plot(timeVec, filterEst(1,outData) + sqrt(filterCovPosX(outData)), 'k');
plot(timeVec, filterEst(1,outData) - sqrt(filterCovPosX(outData)), 'k');
xlabel('Time (s)');
ylabel('X Pos North (m)');
legend('Truth', 'KF Est', '1\sigma Uncert', 'Location', 'Best')
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovPosX(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(filterCovPosX), 2), 'k')
legend('1\sigma Value', 'Mean', 'Location', 'Best')
ylabel('\sigma(1,1)');
ax(axInd) = gca;
axInd = axInd + 1;

% Y position compared to error in NED
figure;
subplot(3,1,1:2)
hold on; grid on;
title('East Position Comparison')
plot(timeVec, nedPos(2,outData), 'b');
plot(timeVec, filterEst(2,outData), 'r');
plot(timeVec, filterEst(2,outData) + sqrt(filterCovPosY(outData)), 'k');
plot(timeVec, filterEst(2,outData) - sqrt(filterCovPosY(outData)), 'k');
xlabel('Time (s)');
ylabel('Y Pos East (m)');
legend('Truth', 'KF Est', '1\sigma Uncert', 'Location', 'Best')
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovPosY(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovPosY)), 2), 'k')
legend('1\sigma Value', 'Mean', 'Location', 'Best')
ylabel('\sigma (2,2)');
ax(axInd) = gca;
axInd = axInd + 1;

% Z position compared to error in NED
figure;
subplot(3,1,1:2)
hold on; grid on;
title('Down Position Comparison')
plot(timeVec, nedPos(3,outData), 'b');
plot(timeVec, filterEst(3,outData), 'r');
plot(timeVec, filterEst(3,outData) + sqrt(filterCovPosZ(outData)), 'k');
plot(timeVec, filterEst(3,outData) - sqrt(filterCovPosZ(outData)), 'k');
xlabel('Time (s)');
ylabel('Z Pos Down (m)');
legend('Truth', 'KF Est', '1\sigma Uncert', 'Location', 'Best')
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovPosZ(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovPosZ)), 2), 'k')
legend('1\sigma Value', 'Mean', 'Location', 'Best')
ylabel('\sigma (3,3)');
ax(axInd) = gca;
axInd = axInd + 1;

% X velocity compared to error in NED
figure;
subplot(3,1,1:2)
hold on; grid on;
title('East Velocity Comparison')
plot(timeVec, nedVel(1,outData), 'b');
plot(timeVec, filterEst(4,outData), 'r');
plot(timeVec, filterEst(4,outData) + sqrt(filterCovVelX(outData)), 'k');
plot(timeVec, filterEst(4,outData) - sqrt(filterCovVelX(outData)), 'k');
xlabel('Time (s)');
ylabel('X Vel East (m)');
legend('Truth', 'KF Est', '1\sigma Uncert', 'Location', 'Best')
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovVelX(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovVelX)), 2), 'k')
legend('1\sigma Value', 'Mean', 'Location', 'Best')
ylabel('\sigma (5,5)');
ax(axInd) = gca;
axInd = axInd + 1;

% Y velocity compared to error in NED
figure;
subplot(3,1,1:2)
hold on; grid on;
title('East Velocity Comparison')
plot(timeVec, nedVel(2,outData), 'b');
plot(timeVec, filterEst(5,outData), 'r');
plot(timeVec, filterEst(5,outData) + sqrt(filterCovVelY(outData)), 'k');
plot(timeVec, filterEst(5,outData) - sqrt(filterCovVelY(outData)), 'k');
xlabel('Time (s)');
ylabel('Y Vel East (m)');
legend('Truth', 'KF Est', '1\sigma Uncert', 'Location', 'Best')
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovVelY(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovVelY)), 2), 'k')
legend('1\sigma Value', 'Mean', 'Location', 'Best')
ylabel('\sigma (5,5)');
ax(axInd) = gca;
axInd = axInd + 1;

% Z velocity compared to error in NED
figure;
subplot(3,1,1:2)
hold on; grid on;
title('Down Velocity Comparison')
plot(timeVec, nedVel(3,outData), 'b');
plot(timeVec, filterEst(6,outData), 'r');
plot(timeVec, filterEst(6,outData) + sqrt(filterCovVelZ(outData)), 'k');
plot(timeVec, filterEst(6,outData) - sqrt(filterCovVelZ(outData)), 'k');
xlabel('Time (s)');
ylabel('Z Vel Down (m)');
legend('Truth', 'KF Est', '1\sigma Uncert', 'Location', 'Best')
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovVelZ(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovVelZ)), 2), 'k')
legend('1\sigma Value', 'Mean', 'Location', 'Best')
ylabel('\sigma (6,6)');
ax(axInd) = gca;
axInd = axInd + 1;

% Roll angle compared to error in Body frame
figure;
subplot(3,1,1:2)
hold on; grid on;
title('Roll Angle Comparison')
plot(timeVec, bodyEulerAng(1,outData), 'b');
plot(timeVec, filterEst(7,outData), 'r');
plot(timeVec, filterEst(7,outData) + sqrt(filterCovRoll(outData)), 'k');
plot(timeVec, filterEst(7,outData) - sqrt(filterCovRoll(outData)), 'k');
xlabel('Time (s)');
ylabel('Roll Angle (deg)');
legend('Truth', 'KF Est', '1\sigma Uncert', 'Location', 'Best')
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovRoll(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovRoll)), 2), 'k')
legend('1\sigma Value', 'Mean', 'Location', 'Best')
ylabel('\sigma (7,7)');
ax(axInd) = gca;
axInd = axInd + 1;

% Pitch angle compared to error in Body frame
figure;
subplot(3,1,1:2)
hold on; grid on;
title('Pitch Angle Comparison')
plot(timeVec, bodyEulerAng(2,outData), 'b');
plot(timeVec, filterEst(8,outData), 'r');
plot(timeVec, filterEst(8,outData) + sqrt(filterCovPitch(outData)), 'k');
plot(timeVec, filterEst(8,outData) - sqrt(filterCovPitch(outData)), 'k');
xlabel('Time (s)');
ylabel('Pitch Angle (deg)');
legend('Truth', 'KF Est', '1\sigma Uncert', 'Location', 'Best')
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovPitch(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovPitch)), 2), 'k')
legend('1\sigma Value', 'Mean', 'Location', 'Best')
ylabel('\sigma (8,8)');
ax(axInd) = gca;
axInd = axInd + 1;

% Yaw angle compared to error in Body frame
figure;
subplot(3,1,1:2)
hold on; grid on;
title('Yaw Angle Comparison')
plot(timeVec, bodyEulerAng(3,outData), 'b');
plot(timeVec, filterEst(9,outData), 'r');
plot(timeVec, filterEst(9,outData) + sqrt(filterCovYaw(outData)), 'k');
plot(timeVec, filterEst(9,outData) - sqrt(filterCovYaw(outData)), 'k');
xlabel('Time (s)');
ylabel('Yaw Angle (deg)');
legend('Truth', 'KF Est', '1\sigma Uncert', 'Location', 'Best')
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3)
hold on; grid on;
plot(timeVec, sqrt(filterCovYaw(outData)), 'r')
plot([timeVec(1) timeVec(end)], repelem(mean(sqrt(filterCovYaw)), 2), 'k')
legend('1\sigma Value', 'Mean', 'Location', 'Best')
ylabel('\sigma (9,9)');
ax(axInd) = gca;
axInd = axInd + 1;

% NED Position errors in subplots
figure;
subplot(3,1,1);
hold on; grid on; grid minor;
plot(timeVec, nedPos(1,outData)-filterEst(1,outData), 'b');
title('Position Errors')
ylabel('North Pos Err (m)');
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,2);
hold on; grid on; grid minor;
plot(timeVec, nedPos(2,outData)-filterEst(2,outData), 'b');
ylabel('East Pos Err (m)');
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3);
hold on; grid on; grid minor;
plot(timeVec, nedPos(3,outData)-filterEst(3,outData), 'b');
ylabel('Down Pos Err (m)');
xlabel('Time (s)');
ax(axInd) = gca;
axInd = axInd + 1;

% NED Velocity errors in subplots
figure;
subplot(3,1,1);
hold on; grid on; grid minor;
plot(timeVec, nedVel(1,outData)-filterEst(4,outData), 'b');
title('Velocity Errors')
ylabel('North Vel Err (m)');
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,2);
hold on; grid on; grid minor;
plot(timeVec, nedVel(2,outData)-filterEst(5,outData), 'b');
ylabel('East Vel Err (m)');
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3);
hold on; grid on; grid minor;
plot(timeVec, nedVel(3,outData)-filterEst(6,outData), 'b');
ylabel('Down Vel Err (m)');
xlabel('Time (s)');
ax(axInd) = gca;
axInd = axInd + 1;

% Body Euler errors in subplots
figure;
subplot(3,1,1);
hold on; grid on; grid minor;
plot(timeVec, bodyEulerAng(1,outData)-filterEst(7,outData), 'b');
title('Euler Angle Errors')
ylabel('Roll Error (deg)');
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,2);
hold on; grid on; grid minor;
plot(timeVec, bodyEulerAng(2,outData)-filterEst(8,outData), 'b');
ylabel('Pitch Error (deg)');
ax(axInd) = gca;
axInd = axInd + 1;
subplot(3,1,3);
hold on; grid on; grid minor;
plot(timeVec, bodyEulerAng(3,outData)-filterEst(9,outData), 'b');
ylabel('Yaw Error (deg)');
xlabel('Time (s)');
ax(axInd) = gca;
axInd = axInd + 1;

linkaxes(ax, 'x')

%% Save figures that have been generated
h =  findobj('type','figure');
n = length(h);
for i = 1:1:n
    saveas(figure(i), sprintf('save_for_later_figure%i', i), 'png' )
    saveas(figure(i), sprintf('save_for_later_figure%i', i), 'fig' )
end