import pandas as pd
import astropy.coordinates as coord
from astroquery.vizier import Vizier
import astropy.units as u

catalogue_names=['I/345/gaia2','I/355/gaiadr3']

row_limit=Vizier.ROW_LIMIT = -1

def renameHeadersDR3(data):
	#Manually rename all columns to their original column names. This one line is 6849 characters long.
	data.rename(columns={'RA_ICRS':'ra','DE_ICRS':'dec','Source':'source_id','e_RA_ICRS':'ra_error','e_DE_ICRS':'dec_error','Plx':'parallax','e_Plx':'parallax_error','PM':'pm','pmRA':'pmra','e_pmRA':'pmra_error','pmDE':'pmdec','e_pmDE':'pmdec_error','RUWE':'ruwe','FG':'phot_g_mean_flux','e_FG':'phot_g_mean_flux_error','Gmag':'phot_g_mean_mag','FBP':'phot_bp_mean_flux','e_FBP':'phot_bp_mean_flux_error','BPmag':'phot_bp_mean_mag','FRP':'phot_rp_mean_flux','e_FRP':'phot_rp_mean_flux_error','RPmag':'phot_rp_mean_mag','BP-RP':'bp_rp','RV':'radial_velocity','e_RV':'radial_velocity_error','Vbroad':'vbroad','GRVSmag':'grvs_mag','QSO':'in_qso_candidates','Gal':'in_galaxy_candidates','NSS':'in_galaxy_candidates','XPcont':'has_xp_continuous','XPsamp':'has_xp_sampled','RVS':'has_rvs','EpochPh':'has_epoch_photometry','EpochRV':'has_epoch_photometry','MCMCGSP':'has_mcmc_gspphot','MCMCMSC':'has_mcmc_msc','And':'in_andromeda_survey','Teff':'teff_gspphot','logg':'logg_gspphot','__Fe_H_':'mh_gspphot','Dist':'distance_gspphot','A0':'distance_gspphot','HIP':'hip_original_ext_source_id','PS1':'ps1_original_ext_source_id','SDSS13':'sdss13_ext_source_id','SKYM2':'skym2_original_ext_source_id','TYC2':'tyc2_original_ext_source_id','URAT1':'urat1_original_ext_source_id','ALLWISE':'allwise_original_ext_source_id','APASS9':'apass9_original_ext_source_id','GSC23':'gsc23_original_ext_source_id','RAVE5':'gsc23_original_ext_source_id','_2MASS':'twomass_original_ext_source_id','RAVE6':'rave6_original_ext_source_id','RAJ2000':'ra2000','DEJ2000':'dec2000','DR3Name':'designation','SolID':'solution_id','RandomI':'random_index','RPlx':'parallax_over_error','RADEcor':'ra_dec_corr','RAPlxcor':'ra_parallax_corr','RApmRAcor':'ra_pmra_corr','RApmDEcor':'ra_pmdec_corr','DEPlxcor':'dec_parallax_corr','DEpmRAcor':'dec_pmra_corr','DEpmDEcor':'dec_pmdec_corr','PlxpmRAcor':'parallax_pmra_corr','PlxpmDEcor':'parallax_pmdec_corr','pmRApmDEcor':'pmra_pmdec_corr','NAL':'astrometric_n_obs_al','NAC':'astrometric_n_obs_ac','NgAL':'astrometric_n_good_obs_al','NbAL':'astrometric_n_bad_obs_al','gofAL':'astrometric_gof_al','chi2AL':'astrometric_chi2_al','epsi':'astrometric_chi2_al','sepsi':'astrometric_chi2_al','Solved':'astrometric_params_solved','APF':'astrometric_primary_flag','nueff':'nu_eff_used_in_astrometry','pscol':'pseudocolour','e_pscol':'pseudocolour_error','RApscolCorr':'ra_pseudocolour_corr','DEpscolCorr':'dec_pseudocolour_corr','PlxpscolCorr':'parallax_pseudocolour_corr','pmRApscolCorr':'pmra_pseudocolour_corr','pmDEpscolCorr':'pmdec_pseudocolour_corr','MatchObsA':'astrometric_matched_transits','Nper':'visibility_periods_used','amax':'astrometric_sigma5d_max','MatchObs':'matched_transits','NewMatchObs':'new_matched_transits','MatchObsrm':'matched_transits_removed','IPDgofha':'ipd_gof_harmonic_amplitude','IPDgofhp':'ipd_gof_harmonic_phase','IPDfmp':'ipd_frac_multi_peak','IPDfow':'ipd_frac_odd_win','SDSk1':'scan_direction_strength_k1','SDSk2':'scan_direction_strength_k2','SDSk3':'scan_direction_strength_k3','SDSk4':'scan_direction_strength_k4','SDMk1':'scan_direction_mean_k1','SDMk2':'scan_direction_mean_k2','SDMk3':'scan_direction_mean_k3','SDMk4':'scan_direction_mean_k4','Dup':'duplicated_source','o_Gmag':'phot_g_n_obs','RFG':'phot_g_mean_flux_over_error','e_Gmag':'phot_g_mean_mag_error','o_BPmag':'phot_bp_n_obs','RFBP':'phot_bp_mean_flux_over_error','e_BPmag':'phot_bp_mean_mag_error','o_RPmag':'phot_rp_n_obs','RFRP':'phot_rp_mean_flux_over_error','e_RPmag':'phot_rp_mean_mag_error','E_BP_RP_':'phot_bp_rp_excess_factor','NBPcont':'phot_bp_n_contaminated_transits','NBPblend':'phot_bp_n_blended_transits','NRPcont':'phot_rp_n_contaminated_transits','NRPblend':'phot_rp_n_blended_transits','Mode':'phot_proc_mode','BP-G':'bp_g','G-RP':'g_rp','n_RV':'rv_method_used','o_RV':'rv_nb_transits','o_RVd':'rv_nb_deblended_transits','RVNper':'rv_visibility_periods_used','RVS_N':'rv_expected_sig_to_noise','RVgof':'rv_renormalised_gof','RVchi2':'rv_chisq_pvalue','RVTdur':'rv_time_duration','RVamp':'rv_amplitude_robust','RVtempTeff':'rv_template_teff','RVtemplogg':'rv_template_logg','RVtemp_Fe_H_':'rv_template_fe_h','Vatmparam':'rv_atm_param_origin','e_Vbroad':'vbroad_error','o_Vbroad':'vbroad_nb_transits','e_GRVSmag':'grvs_mag_error','o_GRVSmag':'grvs_mag_nb_transits','RVSS_N':'rvs_spec_sig_to_noise','GLON':'l','GLAT':'b','ELON':'ecl_lon','ELAT':'ecl_lat','PQSO':'classprob_dsc_combmod_quasar','PGal':'classprob_dsc_combmod_galaxy','PSS':'classprob_dsc_combmod_star','b_Teff':'teff_gspphot_lower','B_Teff':'teff_gspphot_upper','b_logg':'logg_gspphot_lower','B_logg':'logg_gspphot_upper','b__Fe_H_':'mh_gspphot_lower','B__Fe_H_':'mh_gspphot_upper','b_Dist':'distance_gspphot_lower','B_Dist':'distance_gspphot_upper','b_A0':'azero_gspphot_lower','B_A0':'azero_gspphot_upper','AG':'ag_bspphot','b_AG':'ag_gspphot_lower','B_AG':'ag_gspphot_upper','E_BP-RP_':'ebpminrp_gspphot','b_E_BP-RP_':'ebpminrp_gspphot_lower','B_E_BP-RP_':'ebpminrp_gspphot_upper','Lib':'libname_gspphot','dHIP':'hip_angular_distance','nHIP':'hip_number_of_neighbours','f_HIP':'hip_xm_flag','PS1coid':'clean_panstarrs1_oid','dPS1':'ps1_angular_distance','nPS1':'ps1_number_of_neighbours','mPS1':'ps1_number_of_mates','f_PS1':'ps1_xm_flag','SDSS13coid':'clean_sdssdr13_oid','dSDSS13':'sdss13_angular_distance','nSDSS13':'sdss13_number_of_neighbours','mSDSS13':'sdss13_number_of_mates','f_SDSS13':'sdss13_xm_flag','dSKYM2':'skym2_angular_distance','nSKYM2':'skym2_number_of_neighbours','mSKYM2':'skym2_number_of_mates','f_SKYM2':'skym2_xm_flag','dTYC2':'tyc2_angular_distance','f_TYC2':'tyc2_xm_flag','TYC2moid':'tycho2tdsc_merge_oid','nTYC2':'tyc2_number_of_neighbours','dURAT1':'urat1_angular_distance','f_URAT1':'urat1_xm_flag','URAT1oid':'urat1_oid','nURAT1':'urat1_number_of_neighbours','mURAT1':'urat1_number_of_mates','dAllWISE':'allwise_angular_distance','f_AllWISE':'allwise_xm_flag','AllWISEoid':'allwise_oid','nAllWISE':'allwise_number_of_neighbours','mAllWISE':'allwise_number_of_mates','APASS9coid':'clean_apassdr9_oid','dAPASS9':'apass9_angular_distance','nAPASS9':'apass9_number_of_neighbours','mAPASS9':'apass9_number_of_mates','f_APASS9':'apass9_xm_flag','dGSC23':'gsc23_angular_distance','f_GSC23':'gsc23_xm_flag','GSC23coid':'clean_gsc23_oid','nGSC23':'gsc23_number_of_neighbours','mGSC23':'gsc23_number_of_mates','dRAVE5':'rave5_angular_distance','f_RAVE5':'rave5_xm_flag','RAVE5coid':'clean_ravedr5_oid','nRAVE5':'rave5_number_of_neighbours','d2MASS':'twomass_angular_distance','f_2MASS':'twomass_xm_flag','_2MASScoid':'clean_tmass_psc_xsc_oid','n2MASS':'twomass_number_of_neighbours','m2MASS':'twomass_number_of_mates','dRAVE6':'rave6_angular_distance','f_RAVE6':'rave6_xm_flag','RAVE6oid':'ravedr6_oid','nRAVE6':'rave6_number_of_neighbours','e_RAJ2000':'ra2000_error','e_DEJ2000':'dec2000_error','RADEcorJ2000':'ra2000_de2000_corr'}, inplace=True)
	return data
	
#Returns all Gaia data from a given data release for an object, given its Gaia designation
def GaiaQueryDesignation(Designation,catalogue='dr3'):
	#Raises an exception if multiple designations are given, added since GaiaGetDesignation can return a tuple which then breaks if immediately passed back into this function without parsing first.
	if isinstance(Designation,list):
		raise Exception('multiple designations given')
	#'GAIA DR...' prefix in object designations is not present in returned tables, but is required for Vizier query. Therefore need to create two separate parameters, one with the prefix (Search Designation) and one without (Designation)
	if isinstance(Designation,str):
		Designation=Designation.upper()
	catalogue=catalogue.lower()
	if catalogue!='dr2' and catalogue!='dr3':
		raise Exception('invalid catalogue')
		exit()

	if 'GAIA DR3 ' in str(Designation) or 'GAIA DR2 ' in str(Designation):
		SearchDesignation=Designation
		if 'GAIA DR3' in str(Designation):
			catalogue='dr3'
		elif 'GAIA DR2' in str(Designation):
			catalogue='dr2'
		Designation=Designation[9:]
	else:
		if catalogue=='dr2':
			SearchDesignation='GAIA DR2 '+str(Designation)
		elif catalogue=='dr3':
			SearchDesignation='GAIA DR3 '+str(Designation)

	data=[]

	#Sends Vizier query for given designation. Returns multiple objects within a radius, but as you are searching with a GAIA designation it is assumed that you already know the exact object that want and so all others are removed.
	v=Vizier(columns=['**'],column_filters={"Source":"=="+str(Designation)},row_limit=row_limit)
	
	if catalogue=='dr2':
		data.append(v.query_object(SearchDesignation,catalog=catalogue_names[0]))
		if len(data[0])==0:
			print('[Gaia: GaiaQueryDesignation] Error: no Gaia data found for given designation')
			return None
		data=data[0][0].to_pandas()
		data=renameHeadersDR3(data)
		return data
	elif catalogue=='dr3':
		data.append(v.query_object(SearchDesignation,catalog=catalogue_names[1],))
		if len(data[0])==0:
			print('[Gaia: GaiaQueryDesignation] Error: no Gaia data found for given designation')
			return None
		data=data[0][0].to_pandas()
		data=renameHeadersDR3(data)
		return data

#Returns all Gaia data from a given data release for an object, given its coordinates
def GaiaQueryCoords(ObjRa,ObjDec,radius=3,catalogue='dr3'):
	catalogue=catalogue.lower()
	
	if not isinstance(radius,int) and not isinstance(radius,float):
		raise Exception('search radius must be a float or an integer')
	
	if not isinstance(ObjRa,float):
		try:
			ObjRa=float(ObjRa)
		except:
			raise Exception('object RA must be a float or an integer')
	if not isinstance(ObjDec,float):
		try:
			ObjDec=float(ObjDec)
		except:
			raise Exception('object DEC must be a float or an integer')
	
	if catalogue!='dr2' and catalogue!='dr3':
		raise Exception('invalid catalogue')
	
	data=[]
	
	#Get all columns
	v=Vizier(columns=['**'],row_limit=row_limit)
	
	#Sends Vizier query, can return multiple objects in given search radius
	if catalogue=='dr2':
		data.append(v.query_region(coord.SkyCoord(ra=ObjRa,dec=ObjDec,unit=(u.deg,u.deg),frame='icrs'),width=radius*u.arcsec,catalog=catalogue_names[0]))
		if len(data[0])==0:
			print('[Gaia: GaiaQueryCoords] Error: no Gaia data found at given coordinates')
			return None
		data=data[0][0].to_pandas()
		data=data.reset_index(drop=True)
		data=renameHeadersDR3(data)
		return data
	if catalogue=='dr3':
		data.append(v.query_region(coord.SkyCoord(ra=ObjRa,dec=ObjDec,unit=(u.deg,u.deg),frame='icrs'),width=radius*u.arcsec,catalog=catalogue_names[1]))
		if len(data[0])==0:
			print('[Gaia: GaiaQueryCoords] Error: no Gaia data found at given coordinates')
			return None
		data=data[0][0].to_pandas()
		data=data.reset_index(drop=True)
		data=renameHeadersDR3(data)
		return data

#Returns an objects Gaia designation given its coordinates. Designations are given the prefix 'GAIA DR2' or 'GAIA DR3' depending on the catalogue given, as this is the format used by Vizier.
def GaiaGetDesignation(ObjRa,ObjDec,radius=3,catalogue='dr3'):
	data=GaiaQueryCoords(ObjRa,ObjDec,radius,catalogue)
	
	catalogue=catalogue.lower()
	
	Designation=data['source_id']
	#Returns designation if only one object is found
	if len(data)==1:
		if catalogue=='dr2':
			Designation='GAIA DR2 '+str(Designation[0])
		elif catalogue=='dr3':
			Designation='GAIA DR3 '+str(Designation[0])
	#Returns designation as a tuple if multiple objects are found
	if len(data)>1:
		print('[Gaia: GaiaGetDesignation] Note: '+len(data), 'objects detected within search radius.')
		Designation=Designation.to_list()
	if not isinstance(data,pd.DataFrame):
		return None
	return Designation

#Returns an object's coordinates in a given data release given a Gaia designation, returns them as tuple [RA,DEC]
def GaiaGetCoords(Designation,catalogue='dr3'):
	data=GaiaQueryDesignation(Designation,catalogue)
	if not isinstance(data,pd.DataFrame) or data.empty:
		return None
	
	#If used correctly this should never trigger as each object should have a unique Gaia designation
	if len(data)>1:
		raise Exception('too many objects returned by search')
	
	objRa=data['ra']
	objDec=data['dec']
	Coords=[objRa[0],objDec[0]]
	
	return Coords
	
def GaiaGetPhotometryDesignation(Designation,catalogue='dr3'):
	data=GaiaQueryDesignation(Designation,catalogue)
	if not isinstance(data,pd.DataFrame) or data.empty:
		return None
	
	photometry=data[['ra','dec','source_id','phot_g_mean_mag','phot_g_mean_mag_error','phot_bp_mean_mag','phot_bp_mean_mag_error','phot_rp_mean_mag','phot_rp_mean_mag_error']].copy()
	
	return photometry
	
def GaiaGetPhotometryCoords(ObjRa,ObjDec,radius=3):
	data=GaiaQueryCoords(ObjRa,ObjDec,radius=radius)
	
	if not isinstance(data,pd.DataFrame) or data.empty:
		return None
	
	photometry=data[['ra','dec','source_id','phot_g_mean_mag','phot_g_mean_mag_error','phot_bp_mean_mag','phot_bp_mean_mag_error','phot_rp_mean_mag','phot_rp_mean_mag_error']].copy()

	return photometry

if __name__=='__main__':
	print(GaiaQueryDesignation(6050296829033196032))
	
	