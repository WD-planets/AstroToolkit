import pandas as pd 
from astropy.io import fits
from astropy.wcs import WCS
import math
import csv
import numpy as np
import copy

from ..Data.photometry import get_bulkphot_surveys
from ..Misc.file_naming import name_file

'''
creates file according to ATK's format specifications
'''
def create_file(data_copy): 
    if isinstance(data_copy,dict):
        if data_copy['type']=='image':
            file_name=name_file(data_copy,'ATKimage')
            file_type='image'
        elif data_copy['type']=='data':
            file_name=name_file(data_copy,'ATKdata')
            file_type='data'
        elif data_copy['type']=='phot':
            file_name=name_file(data_copy,'ATKphot')
            file_type='phot'
        elif data_copy['type']=='bulkphot':
            file_name=name_file(data_copy,'ATKbulkphot')
            file_type='bulkphot'
        elif data_copy['type']=='sed':
            file_name=name_file(data_copy,'ATKsed')
            file_type='sed'
        elif data_copy['type']=='spectra':
            file_name=name_file(data_copy,'ATKspectrum')
            file_type='spectrum'
    elif isinstance(data_copy,list):
        for element in data_copy:
            if element!=None:
                if element['type']=='lightcurve':
                    file_name=name_file(data_copy,'ATKlightcurve')
                    file_type='lightcurve'
    elif data_copy==None:
        print('Note: None object passed to save_file, suggests no data was found.')
        return None

    # stops input data from actually being edited
    data=copy.deepcopy(data_copy)

    # handles image files
    if file_type=='image':
        if data['pos']!=None:
            ra,dec=data['pos'][0],data['pos'][1]
            
        # makes an astropy fits hdu
        hdu=fits.PrimaryHDU(data['data']['image_data'],header=data['data']['image_header'])

        # sets up some values that ATK uses in fits header
        hdu.header['atk_survey']=data['survey']
        hdu.header['atk_source']=data['source']
        # fits header doesn't support arrays, so split pos into pos_ra and pos_dec
        if data['pos']!=None:
            hdu.header['atk_pos_ra']=data['pos'][0]
            hdu.header['atk_pos_dec']=data['pos'][1]
        else:
            hdu.header['atk_pos_ra']=None
            hdu.header['atk_pos_dec']=None
        hdu.header['atk_location_ra']=data['data']['location'][0]
        hdu.header['atk_location_dec']=data['data']['location'][1]
        hdu.header['atk_size']=data['data']['size']
        hdu.header['atk_time_year']=data['data']['image_time'][0]
        hdu.header['atk_time_month']=data['data']['image_time'][1]
            
        overlay_data=data['data']['overlay']
            
        # appends all overlay data into fits file header as columns which can then be read
        for i in range(0,len(overlay_data)):
            hdu.header[f'atk_overlay_survey_{i}']=overlay_data[i]['survey']
            hdu.header[f'atk_overlay_position_ra_{i}']=overlay_data[i]['position'][0]
            hdu.header[f'atk_overlay_position_dec_{i}']=overlay_data[i]['position'][1]
            hdu.header[f'atk_overlay_radius_{i}']=overlay_data[i]['radius']
            hdu.header[f'atk_overlay_corrected_{i}']=overlay_data[i]['corrected']
            hdu.header[f'atk_overlay_mag_{i}']=overlay_data[i]['mag']
            hdu.header[f'atk_overlay_marker_{i}']=overlay_data[i]['marker']
               
        fname=f'{file_name[:-5]}.fits'

        # create fits file
        hdu.writeto(fname,overwrite=True)
        
    # handles data or phot files
    if file_type=='data' or file_type=='phot':
        source=data['source'] 
            
        if data['source']!=None:
            ra,dec=None,None
        elif data['pos']!=None:
            ra,dec=data['pos'][0],data['pos'][1]
    
        survey=data['survey']
    
        # formats file name
        if data['type']=='data':
            fname=f'{survey}_{file_name[:-5]}.csv'
        else:
            fname=f'{survey}_{file_name[:-5]}.csv'
        
        survey_data=data['data']
        data_type=data['type']
        
        # have to split pos into pos_ra,pos_dec (separate columns) to create dataframe
        survey_data['survey']=survey
        survey_data['type']=data_type
        survey_data['source']=source
        survey_data['pos_ra']=ra
        survey_data['pos_dec']=dec

        # creates dataframe from data
        df=pd.DataFrame.from_dict(survey_data)
        df.to_csv(fname,index=False)

    # handles bulkphot files
    if file_type=='bulkphot':
        source=data['source'] 

        # get file name
        if data['source']!=None:
            ra,dec=None,None
            
        elif data[0]['pos']!=None:
            ra,dec=data['pos'][0],data['pos'][1]
 
        fname=f'{file_name[:-5]}.csv'

        data_arr=[]
        for key in data['data']:
            survey_data=data['data'][key]
            if survey_data!=None:
                data_type=data['type']                

                survey_data['survey']=key
                survey_data['type']=data_type
                survey_data['source']=source
                survey_data['pos_ra']=ra
                survey_data['pos_dec']=dec            
                
                data_arr.append(survey_data)

        # write dataframes to file, one for each survey for which data was available. each survey is separated by a new line (i.e. an empty row)
        with open(fname,'w') as file:
            for element in data_arr:
                df=pd.DataFrame.from_dict(element)
                df.to_csv(file,index=False,lineterminator='\n')
                file.write('\n')
        
    # handles sed files
    if file_type=='sed':
        source=data['source']
            
        if data['source']!=None:
            ra,dec=None,None
            
        elif data['pos']!=None:
            ra,dec=data['pos'][0],data['pos'][1]
               
        # gets file name
        fname=f'{file_name[:-5]}.csv'
            
        data_type=data['type']

        # converts sed data into dataframe which is then written to a file
        df=pd.DataFrame(data['data'])
            
        df['type']=data_type
        df['source']=source
        df['pos_ra']=ra
        df['pos_dec']=dec

        df.to_csv(fname,index=False)
          
    # handles spectra files
    if file_type=='spectrum':
        source=data['source']          

        if data['source']!=None:
            ra,dec=None,None
            
        elif data['pos']!=None:
            ra,dec=data['pos'][0],data['pos'][1]
                
        fname=f'{file_name[:-5]}.csv'
        data_type=data['type']
            
        data_dict={'wavelength':data['data']['wavelength'],'flux':data['data']['flux']}
            
        df=pd.DataFrame(data_dict)
            
        df['type']=data_type
        df['survey']=data['survey']
        df['source']=source
        df['pos_ra']=ra
        df['pos_dec']=dec

        df.to_csv(fname,index=False)

    # if data type is a list (e.g. lightcurves)
    if file_type=='lightcurve':
        source=None
        ra,dec=None,None
        for element in data:
            if element['data']!=None:
                source=element['source']
                pos=element['pos']
        
        if pos==None:
            ra,dec=None,None
        else:
            ra,dec=pos[0],pos[1]

        if source!=None:
            for i in range(0,len(data)):
                if data[i]['data']!=None:
                    del data[i]['pos']
                    data[i]['pos_ra']=None
                    data[i]['pos_dec']=None
            
        elif ra!=None and dec!=None:
            for i in range(0,len(data)):
                if data[i]['data']!=None:
                    del data[i]['pos']
                    data[i]['pos_ra']=ra
                    data[i]['pos_dec']=dec
            
        # collapse 'data' keys into flattened dict
        for i in range(0,len(data)):
            if data[i]['data']!=None:
                for key in data[i]['data']:
                    data[i][key]=data[i]['data'][key]
                del data[i]['data']

        fname=f'{file_name[:-5]}.csv'

        # same as bulkphot (see above)
        with open(fname,'w') as file:
            for element in data:
                try:
                    df=pd.DataFrame.from_dict(element)
                    df.to_csv(file,index=False,lineterminator='\n')
                    file.write('\n')
                except:
                    pass
    
    return fname

'''
reads ATK files, output is lossless (i.e. the exact same data that was used to create the file is retrieved)
'''
def read_file(file_name):
    # reads image files
    if 'ATKimage' in file_name:
        image=fits.open(file_name)[0]
        
        header=image.header
        data=image.data
        
        # gets basic parameters used to reconstruct image data from its header
        survey=header['atk_survey']
        source=header['atk_source']
        if header['atk_pos_ra']==None and header['atk_pos_dec']==None:
            pos=None
        else:
            pos=[header['atk_pos_ra'],header['atk_pos_dec']]
        location=[header['atk_location_ra'],header['atk_location_dec']]
        size=header['atk_size']
        image_time=[header['atk_time_year'],header['atk_time_month']]

        wcs=WCS(header)

        # finds the number of overlay data points by iterating through header, searching for overlay parameters. since they are indexed, this function fails when i = number of data points
        i=0
        while True:
            try:
                overlay_len_check=header[f'atk_overlay_survey_{i}']
                i+=1
            except:
                break
        overlay_len=i
        
        # reconstructs overlay data from header
        overlay_data=[]
        for i in range(0,overlay_len):
            overlay_dict={'survey':header[f'atk_overlay_survey_{i}'],'position':[header[f'atk_overlay_position_ra_{i}'],header[f'atk_overlay_position_dec_{i}']],'radius':header[f'atk_overlay_radius_{i}'],'corrected':header[f'atk_overlay_corrected_{i}'],'mag':header[f'atk_overlay_mag_{i}'],'marker':header[f'atk_overlay_marker_{i}']}
            overlay_data.append(overlay_dict)
            
        image_dict={'type':'image','survey':survey,'source':source,'pos':pos,'data':{'image_data':data,'image_header':header,'location':location,'size':size,'image_time':image_time,'wcs':wcs,'overlay':overlay_data}}

        return image_dict
    
    # reads lightcurve files
    if 'ATKlightcurve' in file_name:
        # reads multiple sets of lightcurve data into a single dataframe
        data=pd.read_csv(file_name)

        # finds which bands are included in the data
        bands=data.band.unique().tolist()
        if 'band' in bands:
            bands.remove('band')

        # splits dataframe into bands
        data_arr=[]
        for band in bands:
            data_element=data[data['band']==band]
            data_arr.append(data_element.reset_index(drop=True))
        
        # reconstructs lightcurve data. here, element = band
        final_arr=[]
        for element in data_arr:
            file_type=element['type'][0]
            source=element['source'][0]
            pos_ra=float(element['pos_ra'][0])
            pos_dec=float(element['pos_dec'][0])
            survey=element['survey'][0]
            band=element['band'][0]
            ra=element['ra'].tolist()
            dec=element['dec'].tolist()
            
            # need to handle both mjd and hjd time formats
            try:
                time=element['hjd'].tolist()
                time_ori=element['hjd_ori'].tolist()
                time_type='hjd'
            except:
                pass
            
            try:
                time=element['mjd'].tolist()
                time_ori=element['mjd_ori'].tolist()
                time_type='mjd'
            except:
                pass
            
            mag=element['mag'].tolist()
            mag_err=element['mag_err'].tolist()

            for i in range(0,len(element)):
                ra[i]=float(ra[i])
                dec[i]=float(dec[i])
                time[i]=float(time[i])
                time_ori[i]=float(time_ori[i])
                mag[i]=float(mag[i])
                mag_err[i]=float(mag_err[i])
                
            pos=[pos_ra,pos_dec]

            # writing dataframe to csv converts None to nan, so need to reverse this
            if math.isnan(pos_ra) and math.isnan(pos_dec):
                pos=None
            elif math.isnan(source):
                source=None
            
            if time_type=='hjd':
                element_dict={'type':file_type,'source':source,'pos':pos,'survey':survey,'band':band,'data':{'ra':ra,'dec':dec,'hjd':time,'hjd_ori':time_ori,'mag':mag,'mag_err':mag_err}}
            elif time_type=='mjd':
                element_dict={'type':file_type,'source':source,'pos':pos,'survey':survey,'band':band,'data':{'ra':ra,'dec':dec,'mjd':time,'mjd_ori':time_ori,'mag':mag,'mag_err':mag_err}}
            final_arr.append(element_dict)
        
        return final_arr

    # reads data and phot files
    if 'ATKdata' in file_name or 'ATKphot' in file_name:
        data=pd.read_csv(file_name)
        
        # essentially does the same as for lightcurves, but doesn't need to handle multiple tables
        file_type=data['type'][0]
        source=data['source'][0]
        pos_ra=float(data['pos_ra'][0])
        pos_dec=float(data['pos_dec'][0])
        survey=data['survey'][0]
        
        pos=[pos_ra,pos_dec]

        if math.isnan(pos_ra) and math.isnan(pos_dec):
            pos=None
        elif math.isnan(source):
            source=None

        # delete atk columns
        data.drop(['type','source','pos_ra','pos_dec','survey'],inplace=True,axis=1)

        data_dict={'survey':survey,'type':file_type,'source':source,'pos':pos,'data':data}

        return data_dict
    
    # reads bulkphot files
    if 'ATKbulkphot' in file_name:
        # this is much more complicated than for lightcurves, because the tables do not have a consistent number of columns
        with open(file_name) as file:
            # reads file line by line
            content=list(csv.reader(file,delimiter=','))
        
        # splits content into tables
        chunks=[]
        chunk=[]
        for i in range(0,len(content)):
            # if current line is a line containing data
            if content[i]:
                chunk.append(content[i])
            # if a blank line is reached, write the current chunk (i.e. a table) to the list of tables (chunks)
            else:
                chunks.append(chunk)
                chunk=[]
        
        # converts empty strings ('') into None values
        for i in range(0,len(chunks)):
            for j in range(0,len(chunks[i])):
                conv=lambda i:i or None
                chunks[i][j]=[conv(i) for i in chunks[i][j]]
        
        # write chunks (tables) into dataframes
        data_arr=[]
        for element in chunks:
            headers=element.pop(0)
            df=pd.DataFrame(element,columns=headers)
            data_arr.append(df)
            
        # reconstruct data
        source=data_arr[0]['source'][0]
        pos_ra=data_arr[0]['pos_ra'][0]
        pos_dec=data_arr[0]['pos_dec'][0]

        pos=[pos_ra,pos_dec]

        if pos_ra==None and pos_dec==None:
            pos=None
        
        survey_list=get_bulkphot_surveys()

        data_dict={}
        for df in data_arr:
            survey=df['survey'][0]
            df.drop(['type','source','pos_ra','pos_dec','survey'],inplace=True,axis=1)
            
            # convert dataframe back to dictionary and set any single-element lists to scalar values
            df_to_dict=pd.DataFrame.to_dict(df,orient='list')
            data_dict[survey]=df_to_dict

        # adds survey:None for surveys that did not return any data 
        for survey in survey_list:
            if survey not in data_dict.keys():
                data_dict[survey]=None

        bulkphot_dict={'type':'bulkphot','source':source,'pos':pos,'data':data_dict}
        
        return bulkphot_dict
    
    # reads sed files
    if 'ATKsed' in file_name:
        df=pd.read_csv(file_name)
        
        # reconstuct data
        dict_arr=[]
        for i in range(0,len(df)):
            data_dict={'survey':df['survey'][i],'wavelength':df['wavelength'][i],'flux':df['flux'][i],'rel_err':df['rel_err'][i]}
            dict_arr.append(data_dict)
        
        file_type=df['type'][0]
        source=df['source'][0]
        pos=[float(df['pos_ra'][0]),float(df['pos_dec'][0])]
        
        if math.isnan(pos[0]) and math.isnan(pos[1]):
            pos=None
        elif math.isnan(source):
            source=None

        sed_dict={'type':file_type,'source':source,'pos':pos,'data':dict_arr}
        
        return sed_dict
    
    # reads spectrum files
    if 'ATKspectrum' in file_name:
        df=pd.read_csv(file_name)
        
        # reconstruct data
        wavelength=df['wavelength'].tolist()
        flux=df['flux'].tolist()
        
        file_type=df['type'][0]
        source=df['source'][0]
        pos=[float(df['pos_ra'][0]),float(df['pos_dec'][0])]
        survey=df['survey'][0]
        
        if math.isnan(pos[0]) and math.isnan(pos[1]):
            pos=None
        elif math.isnan(source):
            source=None
            
        spectrum_dict={'type':file_type,'survey':survey,'source':source,'pos':pos,'data':{'wavelength':np.asarray(wavelength),'flux':np.asarray(flux)}}
        
        return spectrum_dict