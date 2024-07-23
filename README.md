# SORS-Simulations

Monte Carlo Code (run using MCmatlab v4.4.5):
- `Emissions.m` (generates laser excitation distribution in media)
- `EscapeFunctions.m` (generates escape function cylindrical distribution in media)

Raw Data and Processing Code:
- `EscapeFunctions.mat` (resulting data from `EscapeFunctions.m`)
- `Emissions__.mat` (resulting data from `Emissions.m`)
-- on Figshare: 10.6084/m9.figshare.26354497
- `DetectionMatrix.m` (generates Cartesian escape function Cartesian distribution from `EscapeFunctions.mat`)
- `Detections__.mat` (resulting data from `DetectionMatrix.m`)
-- on Figshare: 10.6084/m9.figshare.26354497
- `RamanDistributions.m` (generates collected Raman intensity distribution in media from emissions and escape value distributions for various spatial offsets)
-- Note: resulting data is too large to share here (~150 GB)
