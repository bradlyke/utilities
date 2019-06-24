from astropy.io import fits
import numpy as np

def multiwave_maker(num_rec):
    colarr = np.array(['SDSS_NAME-U18',
                        'RA-float64',
                        'DEC-float64',
                        'PLATE-int32',
                        'MJD-int32',
                        'FIBERID-int16',
                        'MY_CLASS_PQN-U6',
                        'MY_CLASS_ISA-U6',
                        'VI_SET-U1',
                        'IS_QSO_QN-int16',
                        'Z_QN-float64',
                        'RANDOM_SELECT-int16',
                        'Z_10K-float64',
                        'Z_CONF_10K-int16',
                        'PIPE_CORR_10K-int16',
                        'IS_QSO_10K-int16',
                        'PRIM_REC-int16',
                        'THING_ID-int64',
                        'Z_VI-float64',
                        'Z_CONF-int16',
                        'CLASS_PERSON-int16',
                        'Z_DR12Q-float64',
                        'IS_QSO_DR12Q-int16',
                        'Z_DR7Q_SCH-float64',
                        'IS_QSO_DR7Q-int16',
                        'Z_DR7Q_HW-float64',
                        'IS_QSO_FINAL-int16',
                        'Z-float64',
                        'SOURCE_Z-U12',
                        'Z_PIPE-float64',
                        'ZWARNING-int32',
                        'OBJID-U19',
                        'Z_PCA-float64',
                        'ZWARN_PCA-int64',
                        'DELTACHI2_PCA-float64',
                        'Z_HALPHA-float64',
                        'ZWARN_HALPHA-int64',
                        'DELTACHI2_HALPHA-float64',
                        'Z_HBETA-float64',
                        'ZWARN_HBETA-int64',
                        'DELTACHI2_HBETA-float64',
                        'Z_MGII-float64',
                        'ZWARN_MGII-int64',
                        'DELTACHI2_MGII-float64',
                        'Z_CIII-float64',
                        'ZWARN_CIII-int64',
                        'DELTACHI2_CIII-float64',
                        'Z_CIV-float64',
                        'ZWARN_CIV-int64',
                        'DELTACHI2_CIV-float64',
                        'Z_LYA-float64',
                        'ZWARN_LYA-int64',
                        'DELTACHI2_LYA-float64',
                        'Z_DLA-5float64',
                        'NHI_DLA-5float64',
                        'CONF_DLA-5float64',
                        'BAL_PROB-float32',
                        'BI_CIV-float64',
                        'ERR_BI_CIV-float64',
                        'AI_CIV-float64',
                        'ERR_AI_CIV-float64',
                        'BI_SIV-float64',
                        'ERR_BI_SIV-float64',
                        'AI_SIV-float64',
                        'ERR_AI_SIV-float64',
                        'BOSS_TARGET1-int64',
                        'EBOSS_TARGET0-int64',
                        'EBOSS_TARGET1-int64',
                        'EBOSS_TARGET2-int64',
                        'ANCILLARY_TARGET1-int64',
                        'ANCILLARY_TARGET2-int64',
                        'NSPEC_SDSS-int32',
                        'NSPEC_BOSS-int32',
                        'NSPEC-int32',
                        'PLATE_DUPLICATE-99int32',
                        'MJD_DUPLICATE-99int32',
                        'FIBERID_DUPLICATE-99int16',
                        'SPECTRO_DUPLICATE-99int16',
                        'SKYVERSION-B',
                        'RUN_NUMBER-int32',
                        'RERUN_NUMBER-U3',
                        'CAMCOL_NUMBER-int32',
                        'FIELD_NUMBER-int32',
                        'ID_NUMBER-int32',
                        'LAMBDA_EFF-float64',
                        'ZOFFSET-float64',
                        'XFOCAL-float64',
                        'YFOCAL-float64',
                        'CHUNK-U14',
                        'TILE-int32',
                        'PLATESN2-float64',
                        'PSFFLUX-5float32',
                        'PSFFLUX_IVAR-5float64',
                        'PSFMAG-5float32',
                        'PSFMAGERR-5float64',
                        'EXTINCTION-5float32',
                        'M_I-float64',
                        'SN_MEDIAN_ALL-float64',
                        'GALEX_MATCHED-int16',
                        'FUV-float64',
                        'FUV_IVAR-float64',
                        'NUV-float64',
                        'NUV_IVAR-float64',
                        'UKIDSS_MATCHED-int16',
                        'YFLUX-float64',
                        'YFLUX_ERR-float64',
                        'JFLUX-float64',
                        'JFLUX_ERR-float64',
                        'HFLUX-float64',
                        'HFLUX_ERR-float64',
                        'KFLUX-float64',
                        'KFLUX_ERR-float64',
                        'W1_FLUX-float32',
                        'W1_FLUX_IVAR-float32',
                        'W1_CHI2-float32',
                        'W1_SNR-float32',
                        'W1_SRC_FRAC-float32',
                        'W1_EXT_FLUX-float32',
                        'W1_EXT_FRAC-float32',
                        'W1_NPIX-int16',
                        'W1_NEXP-int16',
                        'W1_FIT_NEXP-float32',
                        'W2_FLUX-float32',
                        'W2_FLUX_IVAR-float32',
                        'W2_CHI2-float32',
                        'W2_SNR-float32',
                        'W2_SRC_FRAC-float32',
                        'W2_EXT_FLUX-float32',
                        'W2_EXT_FRAC-float32',
                        'W2_NPIX-int16',
                        'W2_NEXP-int16',
                        'W2_FIT_NEXP-float32',
                        'W3_FLUX-float32',
                        'W3_FLUX_IVAR-float32',
                        'W3_CHI2-float32',
                        'W3_SNR-float32',
                        'W3_SRC_FRAC-float32',
                        'W3_EXT_FLUX-float32',
                        'W3_EXT_FRAC-float32',
                        'W3_NPIX-int16',
                        'W3_NEXP-int16',
                        'W3_FIT_NEXP-float32',
                        'W4_FLUX-float32',
                        'W4_FLUX_IVAR-float32',
                        'W4_CHI2-float32',
                        'W4_SNR-float32',
                        'W4_SRC_FRAC-float32',
                        'W4_EXT_FLUX-float32',
                        'W4_EXT_FRAC-float32',
                        'W4_NPIX-int16',
                        'W4_NEXP-int16',
                        'W4_FIT_NEXP-float32',
                        'FIRST_MATCHED-int16',
                        'FIRST_FLUX-float64',
                        'FIRST_SNR-float64',
                        'SDSS2FIRST_SEP-float64',
                        'JMAG-float64',
                        'JMAG_ERR-float64',
                        'JSNR-float64',
                        'JRDFLAG-int32',
                        'HMAG-float64',
                        'HMAG_ERR-float64',
                        'HSNR-float64',
                        'HRDFLAG-int32',
                        'KMAG-float64',
                        'KMAG_ERR-float64',
                        'KSNR-float64',
                        'KRDFLAG-int32',
                        'SDSS2MASS_SEP-float64',
                        'RASS_COUNTS-float64',
                        'RASS_COUNTS_SNR-float64',
                        'SDSS2ROSAT_SEP-float64',
                        'FLUX_0.2_2.0keV-float32',
                        'FLUX_0.2_2.0keV_ERR-float32',
                        'FLUX_2.0_12.0keV-float32',
                        'FLUX_2.0_12.0keV_ERR-float32',
                        'FLUX_0.2_12.0keV-float32',
                        'FLUX_0.2_12.0keV_ERR-float32',
                        'LUM_0.2_12.0keV-float32',
                        'SDSS2XMM_SEP-float64',
                        'GAIA_MATCHED-int16',
                        'GAIA_DESIGNATION-U28',
                        'GAIA_RA-float64',
                        'GAIA_DEC-float64',
                        'GAIA_RA_ERR-float64',
                        'GAIA_DEC_ERR-float64',
                        'GAIA_PARALLAX-float64',
                        'GAIA_PARALLAX_ERR-float64',
                        'GAIA_PM_RA-float64',
                        'GAIA_PM_RA_ERR-float64',
                        'GAIA_PM_DEC-float64',
                        'GAIA_PM_DEC_ERR-float64',
                        'GAIA_G_FLUX-float64',
                        'GAIA_G_FLUX_ERR-float64',
                        'GAIA_G_MAG-float64',
                        'GAIA_BP_FLUX-float64',
                        'GAIA_BP_FLUX_ERR-float64',
                        'GAIA_BP_MAG-float64',
                        'GAIA_RP_FLUX-float64',
                        'GAIA_RP_FLUX_ERR-float64',
                        'GAIA_RP_MAG-float64',
                        'SDSS2GAIA_SEP-float64'])
    col_names = np.array([])
    col_forms = np.array([])
    for i in range(len(colarr)):
        tt = colarr[i].split('-')
        col_names = np.append(col_names,tt[0])
        col_forms = np.append(col_forms,tt[1])

    estruct = np.zeros(num_rec,dtype={'names':col_names,'formats':col_forms})

    return estruct

def supercat_maker(num_rec):
    colarr = np.array(['SDSS_NAME-U18',
                        'RA-float64',
                        'DEC-float64',
                        'PLATE-int32',
                        'MJD-int32',
                        'FIBERID-int16',
                        'MY_CLASS_PQN-U6',
                        'MY_CLASS_ISA-U6',
                        'VI_SET-U1',
                        'IS_QSO_QN-int16',
                        'Z_QN-float64',
                        'RANDOM_SELECT-int16',
                        'Z_10K-float64',
                        'Z_CONF_10K-int16',
                        'PIPE_CORR_10K-int16',
                        'IS_QSO_10K-int16',
                        'PRIM_REC-int16',
                        'THING_ID-int64',
                        'Z_VI-float64',
                        'Z_CONF-int16',
                        'CLASS_PERSON-int16',
                        'Z_DR12Q-float64',
                        'IS_QSO_DR12Q-int16',
                        'Z_DR7Q_SCH-float64',
                        'IS_QSO_DR7Q-int16',
                        'Z_DR7Q_HW-float64',
                        'IS_QSO_FINAL-int16',
                        'Z-float64',
                        'SOURCE_Z-U12',
                        'Z_PIPE-float64',
                        'ZWARNING-int32',
                        'OBJID-U19',
                        'Z_PCA-float64',
                        'ZWARN_PCA-int64',
                        'DELTACHI2_PCA-float64',
                        'Z_HALPHA-float64',
                        'ZWARN_HALPHA-int64',
                        'DELTACHI2_HALPHA-float64',
                        'Z_HBETA-float64',
                        'ZWARN_HBETA-int64',
                        'DELTACHI2_HBETA-float64',
                        'Z_MGII-float64',
                        'ZWARN_MGII-int64',
                        'DELTACHI2_MGII-float64',
                        'Z_CIII-float64',
                        'ZWARN_CIII-int64',
                        'DELTACHI2_CIII-float64',
                        'Z_CIV-float64',
                        'ZWARN_CIV-int64',
                        'DELTACHI2_CIV-float64',
                        'Z_LYA-float64',
                        'ZWARN_LYA-int64',
                        'DELTACHI2_LYA-float64',
                        'Z_DLA-5float64',
                        'NHI_DLA-5float64',
                        'CONF_DLA-5float64',
                        'BAL_PROB-float32',
                        'BI_CIV-float64',
                        'ERR_BI_CIV-float64',
                        'AI_CIV-float64',
                        'ERR_AI_CIV-float64',
                        'BI_SIV-float64',
                        'ERR_BI_SIV-float64',
                        'AI_SIV-float64',
                        'ERR_AI_SIV-float64',
                        'BOSS_TARGET1-int64',
                        'EBOSS_TARGET0-int64',
                        'EBOSS_TARGET1-int64',
                        'EBOSS_TARGET2-int64',
                        'ANCILLARY_TARGET1-int64',
                        'ANCILLARY_TARGET2-int64',
                        'NSPEC_SDSS-int32',
                        'NSPEC_BOSS-int32',
                        'NSPEC-int32',
                        'PLATE_DUPLICATE-99int32',
                        'MJD_DUPLICATE-99int32',
                        'FIBERID_DUPLICATE-99int16',
                        'SPECTRO_DUPLICATE-99int16',
                        'SKYVERSION-B',
                        'RUN_NUMBER-int32',
                        'RERUN_NUMBER-U3',
                        'CAMCOL_NUMBER-int32',
                        'FIELD_NUMBER-int32',
                        'ID_NUMBER-int32',
                        'LAMBDA_EFF-float64',
                        'ZOFFSET-float64',
                        'XFOCAL-float64',
                        'YFOCAL-float64',
                        'CHUNK-U14',
                        'TILE-int32',
                        'PLATESN2-float64',
                        'PSFFLUX-5float32',
                        'PSFFLUX_IVAR-5float64',
                        'PSFMAG-5float32',
                        'PSFMAGERR-5float64',
                        'EXTINCTION-5float32',
                        'M_I-float64',
                        'SN_MEDIAN_ALL-float64'])

    col_names = np.array([])
    col_forms = np.array([])
    for i in range(len(colarr)):
        tt = colarr[i].split('-')
        col_names = np.append(col_names,tt[0])
        col_forms = np.append(col_forms,tt[1])

    estruct = np.zeros(num_rec,dtype={'names':col_names,'formats':col_forms})

    return estruct

if __name__=='__main__':
    #this is for testing purposes only.
    num_rec = 3
    multiwave_test = multiwave_maker(num_rec)
    multiwave_hdu = fits.BinTableHDU.from_columns(multiwave_test,name='CATALOG')
    multiwave_hdu.writeto('multiwave_test_cat.fits')

    supercat_test = supercat_maker(num_rec)
    supercat_hdu = fits.BinTableHDU.from_columns(supercat_test,name='CATALOG')
    supercat_hdu.writeto('supercat_test.fits')
