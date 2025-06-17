
addpath(genpath('/data/pnl/home/ql087/orchestra-sdk-2.1-1.matlab'))
addpath(genpath('/data/pnl/home/ql087/arrShow-develop'))
addpath(genpath('/data/pnl/home/ql087/Bruker_2022'))

% read the 2-shot GE raw data for calibration

fullpath='/data/pnlx/home/ql087/data_bwh/GE_magnus/Exam235/Series2/'
tmp=strcat(fullpath,'ScanArchive_BWHMAGNUS007MR_20250606_171516753.h5');
pfile = fullfile(tmp);
% archive = GERecon('Archive.Load', pfile);

archive = GERecon('Archive.Load', pfile);
SA_head =  archive.DownloadData.rdb_hdr_rec;
mydiff_head.MB = SA_head.rdb_hdr_mb_factor;
mydiff_head.total_slice = SA_head.rdb_hdr_nslices;
mydiff_head.pass = SA_head.rdb_hdr_npasses;
mydiff_head.diffdir= SA_head.rdb_hdr_numdifdirs;
mydiff_head.logic_slice = mydiff_head.total_slice / mydiff_head.pass;
mydiff_head.b0pass= mydiff_head.pass-mydiff_head.diffdir;
mydiff_head.tensor= SA_head.rdb_hdr_user11;
mydiff_head.ms= SA_head.rdb_hdr_ileaves;
mydiff_head.ExamNumber= archive.ExamNumber;
mydiff_head.pcref_start= SA_head.rdb_hdr_pcref_start;
mydiff_head.pcref_stop= SA_head.rdb_hdr_pcref_stop;
mydiff_head.noise_cal = '';
mydiff_head.asset_cal = '';
mydiff_head.raw_cal = '';
mydiff_head.slthick= archive.DownloadData.rdb_hdr_image.slthick;
mydiff_head.slspace= archive.DownloadData.rdb_hdr_image.scanspacing;
mydiff_head.fov= archive.DownloadData.rdb_hdr_image.dfov;
mydiff_head.dx= archive.DownloadData.rdb_hdr_image.dim_X;
mydiff_head.dy= archive.DownloadData.rdb_hdr_image.dim_Y;
mydiff_head.day= SA_head.rdb_hdr_da_yres -1;
mydiff_head.ctr_R= archive.DownloadData.rdb_hdr_image.ctr_R;
mydiff_head.ctr_A= archive.DownloadData.rdb_hdr_image.ctr_A;
mydiff_head.ctr_S= archive.DownloadData.rdb_hdr_image.ctr_S;
mydiff_head.tableposition= archive.DownloadData.rdb_hdr_series.tablePosition;
mydiff_head.im_size= SA_head.rdb_hdr_im_size;
mydiff_head.RawHeader = SA_head;
mydiff_head.se_no = archive.SeriesNumber;
mydiff_head.se_desc= archive.DownloadData.rdb_hdr_series.se_desc;
mydiff_head.bval = SA_head.rdb_hdr_bvalstab(1);
mydiff_head.sliceorder1 = archive.DownloadData.rdb_hdr_series.start_ras;
mydiff_head.sliceorder2 = archive.DownloadData.rdb_hdr_series.end_ras;
mydiff_head.pepolar= archive.DownloadData.rdb_hdr_image.ihpepolar;
mydiff_head.phasefov=SA_head.rdb_hdr_phase_scale;

for i = 1:mydiff_head.logic_slice
    info = GERecon('Archive.Info', archive, 1, i);
    corners(i)= info.Corners;
    a(i,2) =i;
    a(i,1)= corners(i).UpperLeft(3);
    b= sortrows(a,1);
end

mydiff_head.Orientation = info.Orientation;
% now need sort the cornerpoint in order
for i = 1:mydiff_head.logic_slice
    mydiff_head.corners(i)= corners(b(i,2));
end

mydiff_head.geo_slice= mydiff_head.logic_slice;

% multibandCalHandle = GERecon('Epi.Hyperband.LoadCalibration', mydiff_head.raw_cal);
multibandCalHandle = GERecon('Asset.LoadCalibration', 'Asset-ID171440746-Repetition0000.h5')

for i = 1:mydiff_head.logic_slice
    surfaceImages(:,:,:,i) = GERecon('Calibration.Data',  multibandCalHandle, 'SurfaceImageSpaceAllPassSets', ...
        [mydiff_head.dx; mydiff_head.dy],  mydiff_head.corners(i), mydiff_head.tableposition);
end


