%% readRawCal data
VolumeImages = h5read('./RawCalibration.h5','/VolumeImages');
SurfaceImagesAllPassSets = h5read('./RawCalibration.h5','/SurfaceImagesAllPassSets');
SurfaceImagesPairedPassSets=  h5read('./RawCalibration.h5','/SurfaceImagesPairedPassSets');



% VolumeImages = squeeze(VolumeImages.real + 1i .* VolumeImages.imag);
% SurfaceImagesAllPassSets = squeeze(SurfaceImagesAllPassSets.real + 1i .* SurfaceImagesAllPassSets.imag);
% SurfaceImagesPairedPassSets = squeeze(SurfaceImagesPairedPassSets.real + 1i .* SurfaceImagesPairedPassSets.imag);
info = h5info('ref.h5');

% Display all datasets
for i = 1:length(info.Datasets)
    disp(info.Datasets(i).Name)
end

tmp=h5read('ref.h5','/ChannelsWithMaxSnr');