%SatNav Assignment - AKusmirek [20203284]

%constants and input:
prompt = "Please enter the last four digits of your student number: ";
sn = input(prompt);
c = 299792458;
meanEarthradius = 6371E3;
satellitealtitudes = 20200E3;
r = meanEarthradius + satellitealtitudes;

%Satellite Data from Receiver:
format long
addpath('D:\1. Alexander\2. University College Dublin\6. Applications of Space Science\Assignment 3 Satellite Navigation')
SatData = Receiver(sn);
disp("Satellite data received...")
disp("Calculating...")
%SatData now generated for the four satellites

%Converting SatData to Cartesian coordinate system:
%Change spherical degrees to radians
D = [SatData(1,1), SatData(1,2); SatData(2,1), SatData(2,2); SatData(3,1), SatData(3,2); SatData(4,1), SatData(4,2)];
R = (pi/180).*D;
%Note: for larger problems, this can be done with a "for" loop
%Create array of Cartesian satellite coordinates
xyz1 = [r.*cos(R(1,1)).*cos(R(1,2)), r.*cos(R(1,1)).*sin(R(1,2)), r.*sin(R(1,1))];
xyz2 = [r.*cos(R(2,1)).*cos(R(2,2)), r.*cos(R(2,1)).*sin(R(2,2)), r.*sin(R(2,1))];
xyz3 = [r.*cos(R(3,1)).*cos(R(3,2)), r.*cos(R(3,1)).*sin(R(3,2)), r.*sin(R(3,1))];
xyz4 = [r.*cos(R(4,1)).*cos(R(4,2)), r.*cos(R(4,1)).*sin(R(4,2)), r.*sin(R(4,1))];
SatCart = [xyz1; xyz2; xyz3; xyz4];
disp("Satellite positions determined (Cartesian).")

%Calculating pseudoranges from times of flight:
tof = [SatData(1,3); SatData(2,3); SatData(3,3); SatData(4,3)]; %measured time of flight (seconds)
pr = tof.*c; %pseudoranges (metres)

%Array of converted satellite position values [x, y, z, d]:
SatPos = [SatCart, pr];

%Newton's method:
mu = zeros(4,1); %initial vector of variables
F = zeros(4,1); %initial vector of functions
J = zeros(4,4); %initial Jacobian
tol = 1e-4;
disp("Determining receiver location using Newton's method.")
disp("Calculating Jacobians...")

for n = 1:20
    %Calculate function values and Jacobian for current iteration:
    for i = 1:4
        d = (mu(1)-SatPos(i,1)).^2 + (mu(2)-SatPos(i,2)).^2 + (mu(3)-SatPos(i,3)).^2 - (SatPos(i,4)-mu(4)).^2;
        F(i) = d;
        J1 = [(2*mu(1)-2*SatPos(1,1)),(2*mu(2)-2*SatPos(1,2)),(2*mu(3)-2*SatPos(1,3)),(-2*mu(4)+2*SatPos(1,4))];
        J2 = [(2*mu(1)-2*SatPos(2,1)),(2*mu(2)-2*SatPos(2,2)),(2*mu(3)-2*SatPos(2,3)),(-2*mu(4)+2*SatPos(2,4))];
        J3 = [(2*mu(1)-2*SatPos(3,1)),(2*mu(2)-2*SatPos(3,2)),(2*mu(3)-2*SatPos(3,3)),(-2*mu(4)+2*SatPos(3,4))];
        J4 = [(2*mu(1)-2*SatPos(4,1)),(2*mu(2)-2*SatPos(4,2)),(2*mu(3)-2*SatPos(4,3)),(-2*mu(4)+2*SatPos(4,4))];
        J = [J1; J2; J3; J4];
    end
    
    %Iterate through Newton's method:
    mu = mu - J\F;
end

%Location of receiver in Cartesian coordinates
disp("Receiver location determined (Cartesian).")
recx = mu(1);
recy = mu(2);
recz = mu(3);
reccprerr = mu(4);
recclockoffset = (reccprerr/c)*10.^6; %microseconds

%Converting receiver location to spherical coordinates and degrees
longr = atan2(recy, recx);
latr = atan2(recz, sqrt(recx.^2 + recy.^2));
long = longr.*(180/pi);
lat = latr.*(180/pi);
rr = sqrt(recx.^2 + recy.^2 + recz.^2);
relev = rr - meanEarthradius;

%Code output:
fprintf('Calculated Latitude and Longitude: %f, %f degrees\n', lat, long)
fprintf('Calculated Elevation: %f metres\n', relev)
fprintf('Calculated Receiver Clock Offset: %f microseconds\n', recclockoffset)