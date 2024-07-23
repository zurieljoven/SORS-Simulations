% housekeeping
clearvars;

%% parameters to sweep

% voxel resolutions
nx = 500;
ny = 250;
nz = 250;

% scattering [cm^-1]
muss = [10 20:20:100];

% absorption [cm^-1]
muas = [0.1 0.2:0.2:1];

%% fixed model parameters
model = MCmatlab.model;

% Geometry
model.G.nx = nx; % number of x voxels
model.G.ny = ny; % number of y voxels
model.G.nz = nz; % number of z voxels
model.G.Lx = 10; % length in x [cm]
model.G.Ly = 5; % length in y [cm]
model.G.Lz = 5; % length in z [cm]

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % media properties, see end of file
model.G.geomFunc = @geometryDefinition; % media distribution, see end of file

% Monte Carlo parameters
model.MC.matchedInterfaces = true; % assume all RI equal
model.MC.boundaryType = 1; % all boundaries escaping

model.MC.useLightCollector = false;

model.MC.lightSource.sourceType   = 4; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr = 1; % Radial focal plane intensity distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth = .03; % [cm] Radial focal plane 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.lightSource.angularIntensityDistribution.radialDistr = 1; % Radial angular intensity distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.lightSource.angularIntensityDistribution.radialWidth = 0; % [rad] Radial angular 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom. For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength*1e-9/(pi*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth*1e-2))
model.MC.lightSource.xFocus       = 0; % [cm] x position of focus
model.MC.lightSource.yFocus       = 0; % [cm] y position of focus
model.MC.lightSource.zFocus       = 0; % [cm] z position of focus
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

model.MC.nPhotonsRequested = 5e8;

model.MC.silentMode = true;
model.MC.useGPU = true;

 
for i = 1:numel(muss)
    for j = 1:numel(muas)
        if ~exist("EmissionMatrices/Emissions" + num2str(i) + num2str(j) + ".mat", "file")
            model.G.mediaPropParams = {muss(i), muas(j)};

            disp("mus = " + num2str(muss(i)));
            disp("mua = " + num2str(muas(j)));

            model = runMonteCarlo(model);

            normalizedFluence = model.MC.normalizedFluenceRate;
            save("EmissionMatrices/Emissions" + num2str(i) + num2str(j) + ".mat", "normalizedFluence", "-v7.3");
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
