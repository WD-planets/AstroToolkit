from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
import requests
from io import BytesIO
import pandas as pd
import warnings
import xml.etree.ElementTree as ET
import itertools

from ..Misc.ProperMotionCorrection import PMCorrection

# ignores a pandas warning about some columns being duplicated in returned data and therefore omitted
warnings.simplefilter(action='ignore', category=UserWarning)

# ensure that the row limit in returned data is infinite
row_limit = -1
Vizier.ROW_LIMIT = -1

'''
returns the list of currently supported dataquery/photquery surveys
'''
def get_survey_list():
	survey_list=['gaia','panstarrs','skymapper','galex','rosat','sdss','wise','twomass','erosita']
	return survey_list

'''
returns the times of each supported surveys as a dictionary
'''
def get_survey_times():
	survey_time_dict={'gaia':[2016,0],'panstarrs':[2012,0],'skymapper':[2016,0],'galex':[2007,0],'rosat':[1991,0],'sdss':[2017,0],'wise':[2010,5],'twomass':[1999,0],'ztf':[2019,0],'erosita':[2022,0],'atlas':[2021,0],'gaia_lc':[2016,0],'asassn':[2015,0],'crts':[2000,0]}
	return survey_time_dict

'''
renames all Gaia headers to their original column names (e.g plx --> parallax)
'''
def renameHeadersDR3(data):
	data.rename(columns={'RA_ICRS':'ra','DE_ICRS':'dec','Source':'source_id','e_RA_ICRS':'ra_error','e_DE_ICRS':'dec_error','Plx':'parallax','e_Plx':'parallax_error','PM':'pm','pmRA':'pmra','e_pmRA':'pmra_error','pmDE':'pmdec','e_pmDE':'pmdec_error','RUWE':'ruwe','FG':'phot_g_mean_flux','e_FG':'phot_g_mean_flux_error','Gmag':'phot_g_mean_mag','FBP':'phot_bp_mean_flux','e_FBP':'phot_bp_mean_flux_error','BPmag':'phot_bp_mean_mag','FRP':'phot_rp_mean_flux','e_FRP':'phot_rp_mean_flux_error','RPmag':'phot_rp_mean_mag','BP-RP':'bp_rp','RV':'radial_velocity','e_RV':'radial_velocity_error','Vbroad':'vbroad','GRVSmag':'grvs_mag','QSO':'in_qso_candidates','Gal':'in_galaxy_candidates','NSS':'in_galaxy_candidates','XPcont':'has_xp_continuous','XPsamp':'has_xp_sampled','RVS':'has_rvs','EpochPh':'has_epoch_photometry','EpochRV':'has_epoch_photometry','MCMCGSP':'has_mcmc_gspphot','MCMCMSC':'has_mcmc_msc','And':'in_andromeda_survey','Teff':'teff_gspphot','logg':'logg_gspphot','__Fe_H_':'mh_gspphot','Dist':'distance_gspphot','A0':'distance_gspphot','HIP':'hip_original_ext_source_id','PS1':'ps1_original_ext_source_id','SDSS13':'sdss13_ext_source_id','SKYM2':'skym2_original_ext_source_id','TYC2':'tyc2_original_ext_source_id','URAT1':'urat1_original_ext_source_id','ALLWISE':'allwise_original_ext_source_id','APASS9':'apass9_original_ext_source_id','GSC23':'gsc23_original_ext_source_id','RAVE5':'gsc23_original_ext_source_id','_2MASS':'twomass_original_ext_source_id','RAVE6':'rave6_original_ext_source_id','RAJ2000':'ra2000','DEJ2000':'dec2000','DR3Name':'designation','SolID':'solution_id','RandomI':'random_index','RPlx':'parallax_over_error','RADEcor':'ra_dec_corr','RAPlxcor':'ra_parallax_corr','RApmRAcor':'ra_pmra_corr','RApmDEcor':'ra_pmdec_corr','DEPlxcor':'dec_parallax_corr','DEpmRAcor':'dec_pmra_corr','DEpmDEcor':'dec_pmdec_corr','PlxpmRAcor':'parallax_pmra_corr','PlxpmDEcor':'parallax_pmdec_corr','pmRApmDEcor':'pmra_pmdec_corr','NAL':'astrometric_n_obs_al','NAC':'astrometric_n_obs_ac','NgAL':'astrometric_n_good_obs_al','NbAL':'astrometric_n_bad_obs_al','gofAL':'astrometric_gof_al','chi2AL':'astrometric_chi2_al','epsi':'astrometric_excess_noise','sepsi':'astrometric_excess_noise_sig','Solved':'astrometric_params_solved','APF':'astrometric_primary_flag','nueff':'nu_eff_used_in_astrometry','pscol':'pseudocolour','e_pscol':'pseudocolour_error','RApscolCorr':'ra_pseudocolour_corr','DEpscolCorr':'dec_pseudocolour_corr','PlxpscolCorr':'parallax_pseudocolour_corr','pmRApscolCorr':'pmra_pseudocolour_corr','pmDEpscolCorr':'pmdec_pseudocolour_corr','MatchObsA':'astrometric_matched_transits','Nper':'visibility_periods_used','amax':'astrometric_sigma5d_max','MatchObs':'matched_transits','NewMatchObs':'new_matched_transits','MatchObsrm':'matched_transits_removed','IPDgofha':'ipd_gof_harmonic_amplitude','IPDgofhp':'ipd_gof_harmonic_phase','IPDfmp':'ipd_frac_multi_peak','IPDfow':'ipd_frac_odd_win','SDSk1':'scan_direction_strength_k1','SDSk2':'scan_direction_strength_k2','SDSk3':'scan_direction_strength_k3','SDSk4':'scan_direction_strength_k4','SDMk1':'scan_direction_mean_k1','SDMk2':'scan_direction_mean_k2','SDMk3':'scan_direction_mean_k3','SDMk4':'scan_direction_mean_k4','Dup':'duplicated_source','o_Gmag':'phot_g_n_obs','RFG':'phot_g_mean_flux_over_error','e_Gmag':'phot_g_mean_mag_error','o_BPmag':'phot_bp_n_obs','RFBP':'phot_bp_mean_flux_over_error','e_BPmag':'phot_bp_mean_mag_error','o_RPmag':'phot_rp_n_obs','RFRP':'phot_rp_mean_flux_over_error','e_RPmag':'phot_rp_mean_mag_error','E_BP_RP_':'phot_bp_rp_excess_factor','NBPcont':'phot_bp_n_contaminated_transits','NBPblend':'phot_bp_n_blended_transits','NRPcont':'phot_rp_n_contaminated_transits','NRPblend':'phot_rp_n_blended_transits','Mode':'phot_proc_mode','BP-G':'bp_g','G-RP':'g_rp','n_RV':'rv_method_used','o_RV':'rv_nb_transits','o_RVd':'rv_nb_deblended_transits','RVNper':'rv_visibility_periods_used','RVS_N':'rv_expected_sig_to_noise','RVgof':'rv_renormalised_gof','RVchi2':'rv_chisq_pvalue','RVTdur':'rv_time_duration','RVamp':'rv_amplitude_robust','RVtempTeff':'rv_template_teff','RVtemplogg':'rv_template_logg','RVtemp_Fe_H_':'rv_template_fe_h','Vatmparam':'rv_atm_param_origin','e_Vbroad':'vbroad_error','o_Vbroad':'vbroad_nb_transits','e_GRVSmag':'grvs_mag_error','o_GRVSmag':'grvs_mag_nb_transits','RVSS_N':'rvs_spec_sig_to_noise','GLON':'l','GLAT':'b','ELON':'ecl_lon','ELAT':'ecl_lat','PQSO':'classprob_dsc_combmod_quasar','PGal':'classprob_dsc_combmod_galaxy','PSS':'classprob_dsc_combmod_star','b_Teff':'teff_gspphot_lower','B_Teff':'teff_gspphot_upper','b_logg':'logg_gspphot_lower','B_logg':'logg_gspphot_upper','b__Fe_H_':'mh_gspphot_lower','B__Fe_H_':'mh_gspphot_upper','b_Dist':'distance_gspphot_lower','B_Dist':'distance_gspphot_upper','b_A0':'azero_gspphot_lower','B_A0':'azero_gspphot_upper','AG':'ag_bspphot','b_AG':'ag_gspphot_lower','B_AG':'ag_gspphot_upper','E_BP-RP_':'ebpminrp_gspphot','b_E_BP-RP_':'ebpminrp_gspphot_lower','B_E_BP-RP_':'ebpminrp_gspphot_upper','Lib':'libname_gspphot','dHIP':'hip_angular_distance','nHIP':'hip_number_of_neighbours','f_HIP':'hip_xm_flag','PS1coid':'clean_panstarrs1_oid','dPS1':'ps1_angular_distance','nPS1':'ps1_number_of_neighbours','mPS1':'ps1_number_of_mates','f_PS1':'ps1_xm_flag','SDSS13coid':'clean_sdssdr13_oid','dSDSS13':'sdss13_angular_distance','nSDSS13':'sdss13_number_of_neighbours','mSDSS13':'sdss13_number_of_mates','f_SDSS13':'sdss13_xm_flag','dSKYM2':'skym2_angular_distance','nSKYM2':'skym2_number_of_neighbours','mSKYM2':'skym2_number_of_mates','f_SKYM2':'skym2_xm_flag','dTYC2':'tyc2_angular_distance','f_TYC2':'tyc2_xm_flag','TYC2moid':'tycho2tdsc_merge_oid','nTYC2':'tyc2_number_of_neighbours','dURAT1':'urat1_angular_distance','f_URAT1':'urat1_xm_flag','URAT1oid':'urat1_oid','nURAT1':'urat1_number_of_neighbours','mURAT1':'urat1_number_of_mates','dAllWISE':'allwise_angular_distance','f_AllWISE':'allwise_xm_flag','AllWISEoid':'allwise_oid','nAllWISE':'allwise_number_of_neighbours','mAllWISE':'allwise_number_of_mates','APASS9coid':'clean_apassdr9_oid','dAPASS9':'apass9_angular_distance','nAPASS9':'apass9_number_of_neighbours','mAPASS9':'apass9_number_of_mates','f_APASS9':'apass9_xm_flag','dGSC23':'gsc23_angular_distance','f_GSC23':'gsc23_xm_flag','GSC23coid':'clean_gsc23_oid','nGSC23':'gsc23_number_of_neighbours','mGSC23':'gsc23_number_of_mates','dRAVE5':'rave5_angular_distance','f_RAVE5':'rave5_xm_flag','RAVE5coid':'clean_ravedr5_oid','nRAVE5':'rave5_number_of_neighbours','d2MASS':'twomass_angular_distance','f_2MASS':'twomass_xm_flag','_2MASScoid':'clean_tmass_psc_xsc_oid','n2MASS':'twomass_number_of_neighbours','m2MASS':'twomass_number_of_mates','dRAVE6':'rave6_angular_distance','f_RAVE6':'rave6_xm_flag','RAVE6oid':'ravedr6_oid','nRAVE6':'rave6_number_of_neighbours','e_RAJ2000':'ra2000_error','e_DEJ2000':'dec2000_error','RADEcorJ2000':'ra2000_de2000_corr'}, inplace=True)
	return data

'''
returns Gaia data if input is Gaia source_id
'''
def source_query(source):
	# get all columns, only get specific source requested
	v=Vizier(columns=['**'],column_filters={"Source":"=="+str(source)},row_limit=row_limit)
	
	designation = f'GAIA DR3 {source}'

	data=v.query_object(designation,catalog='I/355/gaiadr3')
	if len(data)==0:
		print('Note: gaia source query returned no data.')
		return None
	
	data=data[0].to_pandas()
	data=renameHeadersDR3(data)
	
	if not data.empty:
		return data	
	else:
		return None

'''
returns survey data if input is a position [ra,dec]
'''
def vizier_query(catalogue_name,ra,dec,radius=3,survey=None):
	data=[]
	
	v=Vizier(columns=['**'],row_limit=row_limit)
	
	data.append(v.query_region(coord.SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),frame='icrs'),width=radius*u.arcsec,catalog=catalogue_name))
	if len(data[0])==0:
		print(f'Note: {survey} query returned no data.')
		return None
	
	data=data[0][0].to_pandas()
	data=data.reset_index(drop=True)
	if catalogue_name=='I/355/gaiadr3':
		data=renameHeadersDR3(data)
	
	if not data.empty:
		return data	
	else:
		return None

'''
returns survey data for surveys that require a separate URL (i.e. don't have data on Vizier)
'''
def url_query(survey,ra,dec,radius=3):
	# convert radius to degrees
	radius=radius/3600	

	if survey=='panstarrs':
		url=f'https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?CAT=PS1dr2OBJECTS&RA={ra}&DEC={dec}&SR={radius}&FORMAT=json'
	elif survey=='skymapper':
		url=f'https://skymapper.anu.edu.au/sm-cone/public/query?RA={ra}&DEC={dec}&SR={radius}&RESPONSEFORMAT=CSV'
	elif survey=='erosita':
		CATALOGUE='DR1_Main'
		VERBOSITY=2
		url=f'https://erosita.mpe.mpg.de/dr1/erodat/catalogue/SCS?CAT={CATALOGUE}&RA={ra}&DEC={dec}&SR={radius}&VERB={VERBOSITY}'

	# get data
	try:
		r=requests.get(url)
	except:
		print(f'Note: Experiencing issues with {survey}.')
		return None
	
	if r.status_code!=200:
		print(f'Note: Experiencing issues with {survey}.')
		return None
	
	# read returned data
	data=BytesIO(r.content)

	# pandas couldn't understand csv return, so had to use json for panstarrs
	if survey=='panstarrs':
		try:
			data=pd.read_json(data)
		except:
			print(f'Note: {survey} dataquery returned no data.')
			return None
	
	elif survey=='skymapper':
		try:
			data=pd.read_csv(data)
		except:
			print(f'Note: {survey} dataquery returned no data.')
			return None
			
	elif survey=='erosita':
		# parse xml table
		tree=ET.parse(data)
		root=tree.getroot()
		
		# move to the position of the field names
		table=root[1][1]
		
		# get a list of field names
		field_names=[]
		for i in range(0,len(table)):
			tag=table[i].tag
			if tag=='FIELD':
				field_name=table[i].attrib['name']
				field_names.append(field_name)
			elif tag=='DATA':
				data_location=i
				break			

		# get a list of lists of field values
		field_values=[]
		values_table=table[data_location][0]
		for row in values_table:
			values=[]
			for value in row:
				# try to convert each value to a float if possible
				try:
					val=float(value.text)
					values.append(val)
				except:
					values.append(value.text)

			# there are 252 potential fields (for VERB=3, 41 for VERB=2, 9 for VERB=1), but for some reason a different number of these have data for each source, i.e. data is missing and so it's impossible to map these
			# effect of this for a sample query (from eROSITA site) was 135 objects returned for VERB=1, 134 objects returned for VERB=2 (i.e. one lost to a missing column), 0 objects returned for VERB=3
			if len(values)==len(field_names):
				field_values.append(values)
			else:
				pass

		# set up data dictionary to be turned into a dataframe
		data_dict={}
		i=0
		for key in field_names:
			value_arr=[]
			for row in field_values:
				value_arr.append(row[i])
		
			data_dict[key]=value_arr
			i+=1
	
		# create dataframe from data_dict
		data=pd.DataFrame.from_dict(data_dict)

	if not data.empty:
		return data
	else:
		print(f'Note: {survey} dataquery returned no data.')
		return None

'''
maps incoming query requests to the correct query method, and corrects for their proper motion before sending the query
'''
def survey_map(survey,pos=None,source=None,radius=3):
	if pos!=None:
		ra,dec=pos[0],pos[1]
	
	elif source!=None:
		# get Gaia data of object (note that since this uses the source_query rather than Tools.dataquery, don't need to use ['data'] here as this is the data before it is appended to a dictionary)
		gaia_data=source_query(source=source)

		if isinstance(gaia_data,pd.DataFrame):
			gaia_data=pd.DataFrame.to_dict(gaia_data,orient='list')	
		else:
			return {'survey':survey,'type':'data','source':source,'pos':pos,'data':None}

		if survey=='gaia':
			data_dict={'survey':survey,'type':'data','source':source,'pos':pos,'data':gaia_data}
			return data_dict
	
		survey_times=get_survey_times()
		
		# correct for proper motion
		ra,dec,pmra,pmdec=gaia_data['ra'][0],gaia_data['dec'][0],gaia_data['pmra'][0],gaia_data['pmdec'][0]
		pos_corrected=PMCorrection(survey_times['gaia'],survey_times[survey],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	
	# send requests
	if survey=='gaia':
		data=vizier_query(catalogue_name='I/355/gaiadr3',ra=ra,dec=dec,radius=radius,survey='gaia')
	elif survey=='galex':
		data=vizier_query(catalogue_name='II/335/galex_ais',ra=ra,dec=dec,radius=radius,survey='galex')
	elif survey=='rosat':
		data=vizier_query(catalogue_name='J/A+A/588/A103',ra=ra,dec=dec,radius=radius,survey='rosat')
	elif survey=='sdss':
		data=vizier_query(catalogue_name='V/154/sdss16',ra=ra,dec=dec,radius=radius,survey='sdss')
	elif survey=='twomass':
		data=vizier_query(catalogue_name='II/246/out',ra=ra,dec=dec,radius=radius,survey='twomass')
	elif survey=='wise':
		data=vizier_query(catalogue_name='II/311/wise',ra=ra,dec=dec,radius=radius,survey='wise')
			
	elif survey=='panstarrs':
		data=url_query(survey='panstarrs',ra=ra,dec=dec,radius=radius)
	elif survey=='skymapper':
		data=url_query(survey='skymapper',ra=ra,dec=dec,radius=radius)
	elif survey=='erosita':
		data=url_query(survey='erosita',ra=ra,dec=dec,radius=radius)

	elif survey=='gaia_lc':
		data=vizier_query(catalogue_name='I/355/epphot',ra=ra,dec=dec,radius=radius,survey='gaia_lc')
	
	# return data as a dataframe if it is lightcurve data (easier to handle)
	if survey=='gaia_lc':
		return data

	if isinstance(data,pd.DataFrame):
		data=pd.DataFrame.to_dict(data,orient='list')

	# set up data dictionary
	data_dict={'survey':survey,'type':'data','source':source,'pos':pos,'data':data}

	return data_dict