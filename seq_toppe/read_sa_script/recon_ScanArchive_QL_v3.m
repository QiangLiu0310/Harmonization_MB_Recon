
addpath(genpath('/data/pnl/home/ql087/orchestra-sdk-2.1-1.matlab'))
addpath(genpath('/data/pnl/home/ql087/arrShow-develop'))
addpath(genpath('/data/pnl/home/ql087/Bruker_2022'))

% read the 2-shot GE raw data for calibration

fullpath='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series2/'
tmp=strcat(fullpath,'ScanArchive_LONGWOOD30MR2_20240409_210858458.h5');
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
        kspace(:,:,1:1:end,i) = squeeze(control.Data);
    else
        fprintf('%d\n', i);
    end

end
GERecon('Archive.Close', archive);
kspace=single(kspace);
tmp=[1:2:90 91:size(kspace,4)];
kspace=kspace(:,:,:,tmp);
save('kspace.mat','kspace','-v7.3')

