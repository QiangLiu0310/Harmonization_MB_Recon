
addpath(genpath('/data/pnl/home/ql087/orchestra-sdk-2.1-1.matlab'))
addpath(genpath('/data/pnl/home/ql087/arrShow-develop'))
addpath(genpath('/data/pnl/home/ql087/Bruker_2022'))

% read the 2-shot GE raw data for calibration

fullpath='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_03_28_ge_sub2/Exam16351/Series14/'
tmp=strcat(fullpath,'ScanArchive_LONGWOOD30MR2_20240329_214020589.h5');
pfile = fullfile(tmp);
archive = GERecon('Archive.Load', pfile);

% Scan parameters
xRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
yRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres - 1;
stop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.stop_rcv;
start = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.start_rcv;
nChannels = stop - start + 1;
phs = archive.Passes; % number of phases

pass = 1;
zRes = archive.SlicesPerPass(pass);

kspace = complex(zeros(xRes,  nChannels, yRes, zRes));

for i = 1:archive.ControlCount

    control = GERecon('Archive.Next', archive);

    if isfield(control, 'Data')
        kspace(:,:,1:2:end,i) = squeeze(control.Data);
    else
        fprintf('%d\n', i);
    end

end
GERecon('Archive.Close', archive);
kspace=single(kspace);
kspace=kspace(:,:,:,1:2:end);
save('kspace_ref.mat','kspace','-v7.3')

