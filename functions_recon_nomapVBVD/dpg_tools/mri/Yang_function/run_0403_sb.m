%% sb
clc;clear;close all;
% [~, directory] = system('pwd');
% directory = strtrim(directory);
% directory_k=[directory '/Raw_K/'];
directory='/rfanfs/pnl-zorro/home/yang/data_paper2/data_0403/';
directory_k='/rfanfs/pnl-zorro/home/yang/data_paper2/data_0403/Raw_K/';
%% 
sb_TE1  = 'sb_Part1_TE79';
% sb_TE2  = 'sb_Part1_TE108';
% sb_TE3  = 'sb_Part1_TE137';

k_sb_TE1 = mapVBVD(   findfile(directory_k,sb_TE1)  );
% k_sb_TE2 = mapVBVD(   findfile(directory_k,sb_TE2)  );
% k_sb_TE3 = mapVBVD(   findfile(directory_k,sb_TE3)  );

% k_sb_TE1_PA = mapVBVD(   findfile_PA(directory_k,sb_TE1)  );
% k_sb_TE2_PA = mapVBVD(   findfile_PA(directory_k,sb_TE2)  );
% k_sb_TE3_PA = mapVBVD(   findfile_PA(directory_k,sb_TE3)  );
% 
[img_sb_TE1]  = recon_std_EPI_QL_debug1(k_sb_TE1);
% [img_sb_TE2]  = recon_std_EPI(k_sb_TE2);
% [img_sb_TE3]  = recon_std_EPI(k_sb_TE3);
% 
% [img_sb_TE2_PA]  = recon_std_EPI(k_sb_TE2_PA);
% [img_sb_TE3_PA]  = recon_std_EPI(k_sb_TE3_PA);
% 
% 
% img3_TE1 = img_sb_TE1;
% img3_TE2 = img_sb_TE2;
% img3_TE3 = img_sb_TE3;
% 
% img3_TE1_PA_B0 = img_sb_TE1_PA;
% img3_TE2_PA_B0 = img_sb_TE2_PA;
% img3_TE3_PA_B0 = img_sb_TE3_PA;
% 
% save([directory_k 'sb.mat']);
%% 
% dir=[ '/rfanfs/pnl-zorro/home/yang/data_paper2/0215/data_0215/NING_21_02_15-17_18_51-STD-1_3_12_2_1107_5_2_43_66109/__20210215_172113_815000'];
% 
% seq1='SB_STD_PART1_TE79_0011';
% path_seq1=[dir  '/' seq1];
% stdSB_P1_TE1S=reshape(readDicom_SS(path_seq1),[100 100 20 76]);
% 
% seq2='SB_STD_PART2_TE79_0017';
% path_seq2=[dir  '/' seq2];
% stdSB_P2_TE1S=reshape(readDicom_SS(path_seq2),[100 100 20 76]);
% 
% seq3='SB_STD_PART3_TE79_0023';
% path_seq3=[dir  '/' seq3];
% stdSB_P3_TE1S=reshape(readDicom_SS(path_seq3),[100 100 20 76]);
% 
% 
% seq11='SB_STD_PART1_TE108_0012';
% path_seq11=[dir  '/' seq11];
% stdSB_P1_TE2S=reshape(readDicom_SS(path_seq11),[100 100 20 76]);
% 
% seq22='SB_STD_PART2_TE108_0018';
% path_seq22=[dir  '/' seq22];
% stdSB_P2_TE2S=reshape(readDicom_SS(path_seq22),[100 100 20 76]);
% 
% seq33='SB_STD_PART3_TE108_0024';
% path_seq33=[dir  '/' seq33];
% stdSB_P3_TE2S=reshape(readDicom_SS(path_seq33),[100 100 20 76]);
% 
% seq111='SB_STD_PART1_TE137_0013';
% path_seq111=[dir  '/' seq111];
% stdSB_P1_TE3S=reshape(readDicom_SS(path_seq111),[100 100 20 76]);
% 
% seq222='SB_STD_PART2_TE137_0019';
% path_seq222=[dir  '/' seq222];
% stdSB_P2_TE3S=reshape(readDicom_SS(path_seq222),[100 100 20 76]);
% 
% seq333='SB_STD_PART3_TE137_0025';
% path_seq333=[dir  '/' seq333];
% stdSB_P3_TE3S=reshape(readDicom_SS(path_seq333),[100 100 20 76]);


