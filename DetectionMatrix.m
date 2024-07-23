clearvars;

load("EscapeFunctions.mat");

%% parameters to sweep

% lengths [cm]
Lx = 5;
Ly = 5;
Lz = 5;

% voxels
Nx = 250;
Ny = 250;
Nz = 250;

% depth levels (number of z voxels)
depths = 50;
z_low = linspace(0, Lz-Lz/depths, depths);
z_upp = linspace(Lz/depths, Lz, depths);
z_mid = (z_low+z_upp)/2;

% rings of radial distances [cm]
r_low = 0:0.2:2;
r_upp = [0.2:0.2:2 Inf];
r = [0.1:0.2:1.9 2.25];

% centers of voxels
x = linspace(-(Lx-Lx/Nx)/2, (Lx-Lx/Nx)/2, Nx);
y = linspace(-(Ly-Ly/Ny)/2, (Ly-Ly/Ny)/2, Ny);
z = 0.5*(Lz/Nz):(Lz/Nz):(Lz - 0.5*Lz/Nz);

% coordinates
[X, Y, Z] = ndgrid(x, y, z);
R = sqrt(X.^2 + Y.^2);

% scattering [cm^-1]
muss = [0.1 1 10 20:20:100];

% absorption [cm^-1]
muas = [0.01 0.1 0.2:0.2:1 5];

for i = 1:6
    for j = 1:6
        if ~exist("Detections" + num2str(i) + num2str(j) + ".mat", "file")
            disp("mus = " + num2str(muss(i)));
            disp("mua = " + num2str(muas(j)));

            E = griddedInterpolant({r, z_mid}, squeeze(escapeFunction(:,:,i,j)).', "linear", "nearest");
            detectionMatrix = E(R, Z);

            save("Detections" + num2str(i) + num2str(j) + ".mat", "detectionMatrix", "-v7.3");
        end
    end
end