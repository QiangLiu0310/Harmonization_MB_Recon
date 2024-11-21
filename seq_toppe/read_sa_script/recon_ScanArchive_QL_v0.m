fullpath='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_03_13_ge/Exam16107/Series7/'
tmp=strcat(fullpath,'ScanArchive_LONGWOOD30MR2_20240313_201603515.h5');
pfile = fullfile(tmp);
archive = GERecon('Archive.Load', pfile);

% Scan parameters
xRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
yRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres - 1;
stop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.stop_rcv;
start = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.start_rcv;
nChannels = stop - start + 1;
phs = archive.Passes; % number of phases

% Keep track of the current pass
pass = 1;
zRes = archive.SlicesPerPass(pass);

% Allocate K-space
% kspace = complex(zeros(xRes, yRes, nChannels, zRes));
kspace = complex(zeros(xRes,  nChannels, yRes, zRes));

table=[]
% Loop through each control, sorting frames if applicable
for i = 1:archive.ControlCount

    % Retrieve the next control packet along with its associated data (if
    % applicable) from the archive. Frame data is returned in the
    % control.Data field and is organized as: ReadoutSize x NumChannels x
    % NumFrames where NumFrames is the number of frames corresponding to
    % this control packet. Each programmable packet corresponds to a single
    % frame. Thus, for this example, the frames dimension of control.Data
    % will always have a size of 1
    control = GERecon('Archive.Next', archive);

    % Sort only programmable packets in range
    %     if(control.opcode == 1 && ...
    %        control.viewNum > 0 && ...
    %        control.viewNum <= yRes && ...
    %        control.sliceNum < zRes)

    % Each programmable packet contains a single frame of data. Squeeze
    % off the singular frames dimension and copy the data into the
    % kSpace matrix
    if isfield(control, 'Data')
        kspace(:,:,:,i) = squeeze(control.Data);
    else
        fprintf('%d\n', i);
    end

    %     elseif(control.opcode == 0) % end of pass and/or scan
    %
    %         % Reconstruct this pass
    %         for slice = 1:zRes
    %             for channel = 1:nChannels
    %                 channelImages(:,:,channel) = GERecon('Transform', kspace(:,:,channel,slice));
    %             end
    %
    %             % Apply Channel Combination
    %             combinedImage = GERecon('SumOfSquares', channelImages);
    %
    %             % Create Magnitude Image
    %             magnitudeImage = abs(combinedImage);
    %
    %             % Get Info
    %             info = GERecon('Archive.Info', archive, pass, slice);
    %
    %             % Orient the image
    %             finalImage = GERecon('Orient', magnitudeImage, info.Orientation);
    %
    %             % meic complex images
    %             ksp(:, :, channel, slice, pass) = kspace;
    %             ima(:, :, channel, slice, pass) = ifftshift(ifft(ifft(ifftshift(ifftshift(kspace,1),2),[],1),[],2),1);

    % Display
    %meic imagesc(finalImage);
    %figure(111);hold on
    %im(ima);
    %title(['Pass: ' num2str(pass) 'Slice: ' num2str(info.Number)]);
    %pause(0.5);

    % Save DICOMs
    %filename = ['DICOMs/image' num2str(info.Number) '.dcm'];
    %GERecon('Dicom.Write', filename, finalImage, info.Number, info.Orientation, info.Corners);
end

% Move the next pass and clear out kspace if not last pass
%         if(pass < archive.Passes)
%             pass = pass + 1;
%             zRes = archive.SlicesPerPass(pass);
%             kspace = complex(zeros(xRes, yRes, nChannels, zRes));
%         end

%     end
% end

% Close the archive to release it from memory
GERecon('Archive.Close', archive);