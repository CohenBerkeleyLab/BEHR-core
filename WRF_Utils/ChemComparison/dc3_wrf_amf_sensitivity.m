function [ amfs_dc3, amfs_wrf, amfs_dc3_ft_wrf_bl, amfs_dc3_mt_wrf_bl_ut, amfs_dc3_rl_wrf_bl_ft, amfs_wrf_ft_dc3_bl, amfs_wrf_mt_dc3_bl_ut ] = dc3_wrf_amf_sensitivity( data_file, skytype )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

us_clon = -95;
us_clat = 37.5;
us_month = 6;

data = load(data_file);
data = data.profs;

warn_state=warning('off');
amfs_dc3 = amfSensitivityTest_par(data.dc3_prof, data.dc3_pres, us_clon, us_clat, us_month, skytype);
amfs_wrf = amfSensitivityTest_par(data.wrf_prof, data.wrf_pres, us_clon, us_clat, us_month, skytype);
amfs_dc3_mt_wrf_bl_ut = amfSensitivityTest_par(data.dc3_mt_wrf_bl_ut, data.dc3_pres, us_clon, us_clat, us_month, skytype);
amfs_dc3_rl_wrf_bl_ft = amfSensitivityTest_par(data.dc3_rl_wrf_bl_ft, data.dc3_pres, us_clon, us_clat, us_month, skytype);
amfs_dc3_ft_wrf_bl = amfSensitivityTest_par(data.dc3_ft_wrf_bl, data.dc3_pres, us_clon, us_clat, us_month, skytype);
amfs_wrf_ft_dc3_bl = amfSensitivityTest_par(data.wrf_ft_dc3_bl, data.wrf_pres, us_clon, us_clat, us_month, skytype);
amfs_wrf_mt_dc3_bl_ut = amfSensitivityTest_par(data.wrf_mt_dc3_bl_ut, data.wrf_pres, us_clon, us_clat, us_month, skytype);
warning(warn_state);
end

