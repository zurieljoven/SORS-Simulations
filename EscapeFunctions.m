%% EscapeFunctions.m -- Zuriel Joven
% Generate escape function tensor E = f(z, r, mu_a, mu_s)
% Must be run in MCmatlab-Release folder
% Generation progress is saved between simulations

% housekeeping
clearvars;

if exist("EscapeFunctions.mat", "file")
    load("EscapeFunctions.mat");
end

%% parameters to sweep

% depth levels (number of z voxels)
depths = 50;

% rings of radial distances [cm]
r_low = 0:0.2:2;
r_upp = [0.2:0.2:2 Inf];

% scattering [cm^-1]
muss = [10 20:20:100];

% absorption [cm^-1]
muas = [0.1 0.2:0.2:1];

% initialize escape function tensor E = f(z, r, mu_s, mu_a)
if ~exist("escapeFunction", "var")
    escapeFunction = ones(depths, numel(r_low), numel(muss), numel(muas));
    save("EscapeFunctions.mat", "escapeFunction", "-v7.3"); % save result to .mat file
end

%% fixed model parameters
model = MCmatlab.model;

% Geometry
model.G.nx = 50; % number of x voxels
model.G.ny = 50; % number of y voxels
model.G.nz = depths; % number of z voxels
model.G.Lx = 5; % length in x [cm]
model.G.Ly = 5; % length in y [cm]
model.G.Lz = 5; % length in z [cm]

[X, Y] = ndgrid(model.G.x, model.G.y); % store coordinates

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % media properties, see end of file
model.G.geomFunc = @geometryDefinition; % media distribution, see end of file

% Monte Carlo parameters
model.MC.matchedInterfaces = true; % assume all RI equal
model.MC.boundaryType = 1; % all boundaries escaping

model.MC.calcNormalizedFluenceRate = false;

model.MC.useLightCollector = true;
model.MC.lightCollector.x = 0; % x position [cm]
model.MC.lightCollector.y = 0; % y position [cm]
model.MC.lightCollector.z = -0.001; % z position [cm]

model.MC.lightCollector.theta = 0; % polar angle [rad]
model.MC.lightCollector.phi = pi/2; % azimuthal angle [rad]

model.MC.lightCollector.f = Inf; % fiber collector
model.MC.lightCollector.diam = 0.1; % diameter of collector aperature [cm]
model.MC.lightCollector.fieldSize = .1;
model.MC.lightCollector.NA = 0.22; % fiber numerical aperature

model.MC.lightCollector.res = 1; % fiber collector = 1 px resolution

%model.MC.depositionCriteria.onlyCollected = true; % record photon packets collected by detector

model.MC.silentMode = true;
model.MC.useGPU = true;

%% first, try with 5e6 photons

model.MC.nPhotonsRequested = 5e6;

for z = 1:depths
    for ring = 1:numel(r_low)
        % source distribution = ring of emitters (1s) at single z-level
        model.MC.sourceDistribution = zeros([size(X) depths]);
        model.MC.sourceDistribution(:,:,z) = (X.^2 + Y.^2 >= r_low(ring)^2).*(X.^2 + Y.^2 < r_upp(ring)^2);

        for j = 1:numel(muss)
            for k = 1:numel(muas)
                if escapeFunction(z, ring, j, k) == 1
                    model.G.mediaPropParams = {muss(j), muas(k)}; % pass into mediaPropertiesFunc function
                    disp("Simulations Left: " + num2str(nnz(escapeFunction == 1) + 1));
                    disp("Depth = " + num2str(z));
                    disp("Ring = " + num2str(ring));
                    disp("mus = " + num2str(muss(j)));
                    disp("mua = " + num2str(muas(k)));
                    model = runMonteCarlo(model); % run simulation
                    if model.MC.lightCollector.image > 0
                        escapeFunction(z, ring, j, k) = model.MC.lightCollector.image; % extract % of incident light that hits detector
                    else
                        escapeFunction(z, ring, j, k) = -1;
                    end
                    save("EscapeFunctions.mat", "escapeFunction", "-v7.3"); % save result to .mat file
                end
            end
        end
    end
end

%% if 0 photons collected, retry with 5e7 photons

model.MC.nPhotonsRequested = 5e7;

for z = 1:depths
    for ring = 1:numel(r_low)
        % source distribution = ring of emitters (1s) at single z-level
        model.MC.sourceDistribution = zeros([size(X) depths]);
        model.MC.sourceDistribution(:,:,z) = (X.^2 + Y.^2 >= r_low(ring)^2).*(X.^2 + Y.^2 < r_upp(ring)^2);

        for j = 1:numel(muss)
            for k = 1:numel(muas)
                if escapeFunction(z, ring, j, k) == -1
                    model.G.mediaPropParams = {muss(j), muas(k)}; % pass into mediaPropertiesFunc function
                    disp("Redo Simulations Left: " + num2str(nnz(escapeFunction < 0) + 1));
                    disp("Depth = " + num2str(z));
                    disp("Ring = " + num2str(ring));
                    disp("mus = " + num2str(muss(j)));
                    disp("mua = " + num2str(muas(k)));
                    model = runMonteCarlo(model); % run simulation
                    escapeFunction(z, ring, j, k) = model.MC.lightCollector.image; % extract % of incident light that hits detector
                    save("EscapeFunctions.mat", "escapeFunction", "-v7.3"); % save result to .mat file
                end
            end
        end
    end
end

%% Geometry function -- homogeneous tissue
function M = geometryDefinition(X, Y, Z, parameters)
    M = ones(size(X));
end

%% Media function
function mediaProperties = mediaPropertiesFunc(parameters)
    mediaProperties = MCmatlab.mediumProperties;

    j = 1;
    mediaProperties(j).name = "tissue";
    mediaProperties(j).mus = parameters{1};
    mediaProperties(j).mua = parameters{2};
    mediaProperties(j).g = 0.9;
    mediaProperties(j).n = 1.4;
    mediaProperties(j).VHC = 4.19;
    mediaProperties(j).TC = 5.8e-3;
end