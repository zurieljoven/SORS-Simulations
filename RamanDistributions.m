clearvars;

%% parameters

% spatial offsets [cm]
offsets = 0:0.02:2;

% rings of radial distances [cm]
r_low = 0:0.2:2;
r_upp = [0.2:0.2:2 Inf];

% scattering [cm^-1]
muss = [10 20:20:100];

% absorption [cm^-1]
muas = [0.1 0.2:0.2:1];

% z layers
depths = 250;

% voxels
nx = 250;
ny = 250;
Nx = 500;
Ny = 250;

for i = 1:numel(muss)
    for j = 1:numel(muas)
        disp("mus = " + num2str(muss(i)));
        disp("mua = " + num2str(muas(j)));
        if i*j > 1
            clear detectionMatrix normalizedFluence;
        end
        load("EmissionMatrices/Emissions" + num2str(i) + num2str(j) + ".mat");
        load("DetectionMatrices/Detections" + num2str(i) + num2str(j) + ".mat");
        for deltaX = 0:(numel(offsets)-1)
            leftX = (Nx+1)/2 - (nx-1)/2 + deltaX;
            rightX = (Nx+1)/2 + (nx-1)/2 + deltaX;
            leftY = (Ny+1)/2 - (ny-1)/2;
            rightY = (Ny+1)/2 + (ny-1)/2;
            
            if exist("RamanDistribution/Material" + num2str(i) + num2str(j), "dir") == 0
                mkdir("RamanDistribution/Material" + num2str(i) + num2str(j));
            end

            if exist("RamanDistribution/Material" + num2str(i) + num2str(j) + "/CollectedRaman" + num2str(deltaX+1) + ".mat", "file") == 0
                CollectedRamanDistribution = normalizedFluence(leftX:rightX,leftY:rightY,:).*detectionMatrix;

                borderingIntensity = 1 - sum(CollectedRamanDistribution(2:end-1,2:end-1,2:end-1), "all")/sum(CollectedRamanDistribution, "all");
                borderingNonzeroVoxels = 1 - nnz(CollectedRamanDistribution(2:end-1,2:end-1,2:end-1))/nnz(CollectedRamanDistribution);

                save("RamanDistribution/Material" + num2str(i) + num2str(j) + "/CollectedRaman" + num2str(deltaX+1) + ".mat", "CollectedRamanDistribution", "borderingIntensity", "borderingNonzeroVoxels", "-v7.3");
            end
        end
    end
end

