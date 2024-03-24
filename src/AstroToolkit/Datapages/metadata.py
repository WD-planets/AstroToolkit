from bokeh.models import ColumnDataSource,DataTable,TableColumn
import pandas as pd

from ..Tools import dataquery
from ..Data.data import get_survey_list

'''
appends data points to the lists in list_dict
'''
def getdata(list_dict,sub_dict,survey,source,pos,default):
    supported_surveys=get_survey_list()    

    # splits list_dict into its individual lists
    survey_list=list_dict['survey_list']
    parameter_list=list_dict['parameter_list']
    value_list=list_dict['value_list']
    error_list=list_dict['error_list']
    note_list=list_dict['note_list']

    # if key is a survey, grabs the data for that survey
    if survey in supported_surveys:
        if source!=None:
            data=dataquery(survey=survey,source=source)['data']
            if data==None:
                return list_dict
        elif pos!=None:
            data=dataquery(survey=survey,pos=pos)['data']
            if data==None:
                return list_dict
       
    # if a key (survey) has 'default' as its value, sets the default sub_dict used to create the data table
    if default==True:
        if survey=='gaia':
            sub_dict={'parameters':['source_id','ra','dec','pmra','pmdec','parallax','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag'],
                      'errors':[None,'ra_error','dec_error','pmra_error','pmdec_error','parallax_error','phot_g_mean_mag_error','phot_bp_mean_mag_error','phot_rp_mean_mag_error'],
                      'notes':['source id','right ascension [deg]','declination [deg]','Proper motion in RA direction [mas/yr]','Proper motion in DEC direction [mas/yr]','parallax [mas]','g mag','bp mag','rp mag']
                      }
            
        if survey=='panstarrs':
            sub_dict={'parameters':['gMeanPSFMag','rMeanPSFMag','iMeanPSFMag','zMeanPSFMag','yMeanPSFMag'],
                      'errors':['gMeanPSFMagErr','rMeanPSFMagErr','iMeanPSFMagErr','zMeanPSFMagErr','yMeanPSFMagErr'],                
                      'notes':['g mag','r mag','i mag','z mag','y mag']
                      }
        
        if survey=='skymapper':
            sub_dict={'parameters':['g_psf','r_psf','i_psf','z_psf','u_psf','v_psf'],
                      'errors':['e_g_psf','e_r_psf','e_i_psf','e_z_psf','e_u_psf','e_v_psf'],                
                      'notes':['g mag','r mag','i mag','z mag','u mag','v mag']
                      }
            
        if survey=='galex':
            sub_dict={'parameters':['NUVmag','FUVmag'],
                      'errors':['e_NUVmag','e_FUVmag'],                
                      'notes':['FUV mag','NUV mag']
                      }
            
        if survey=='rosat':
            sub_dict=None
        
        if survey=='sdss':
            sub_dict={'parameters':['gPmag','rPmag','iPmag','zPmag','uPmag'],
                      'errors':['e_gPmag','e_rPmag','e_iPmag','e_zPmag','e_uPmag'],                
                      'notes':['g mag','r mag','i mag','z mag','u mag']
                      }
            
        if survey=='wise':
            sub_dict={'parameters':['W1mag','W2mag','W3mag','W4mag'],
                      'errors':['e_W1mag','e_W2mag','e_W3mag','e_W4mag'],                
                      'notes':['W1 mag','W2 mag','W3 mag','W4 mag']
                      }
                    
        if survey=='twomass':
            sub_dict={'parameters':['Jmag','Hmag','Kmag'],
                      'errors':['e_Jmag','e_Hmag','e_Kmag'],                
                      'notes':['J mag','H mag','K mag']
                      }
    
    # used to ignore any surveys for which a default sub_dict is not supported (e.g. ROSAT)
    if sub_dict==None:
        return list_dict
    
    # grabs parameters from data if the key was a survey
    if survey in supported_surveys:
        for parameter in sub_dict['parameters']:
            survey_list.append(survey)
            parameter_list.append(parameter)
            value_list.append(str(data[parameter][0]))
        for error in sub_dict['errors']:
            if error!=None:
                error_list.append(str(data[error][0]))
            else:
                error_list.append('---')
        for note in sub_dict['notes']:
            if note!=None:
                note_list.append(note)
            else:
                note_list.append('---')
                
    # otherwise, just appends the parameters and values given directly
    else:
        for parameter in sub_dict['parameters']:
            survey_list.append(survey)
            parameter_list.append(parameter)
        for value in sub_dict['values']:
            value_list.append(value)
        for error in sub_dict['errors']:
            if error!=None:
                error_list.append(error)
            else:
                error_list.append('---')
        for note in sub_dict['notes']:
            if note!=None:
                note_list.append(note)
            else:
                note_list.append('---')

    list_dict={'survey_list':survey_list,'parameter_list':parameter_list,'value_list':value_list,'error_list':error_list,'note_list':note_list}

    return list_dict

'''
cycle through keys (surveys/custom) in metadata_dict, and grabs their data before creating a table from the result
'''
def gettable(metadata_dict,pos=None,source=None):
    from ..Misc.file_naming import convertsource

    supported_surveys=get_survey_list()
    
    # get identifier if a source was given as input
    if source!=None:
        identifier=convertsource(source)
    else:
        identifier=None

    # create data lists
    survey_list=[]
    parameter_list=[]
    value_list=[]
    error_list=[]
    note_list=[]

    # add some rows to the metadata table that show information used by the toolkit
    
    survey_list.append('ATK')
    parameter_list.append('ATK source')
    value_list.append(str(source))
    error_list.append('---')
    note_list.append('ATK input source')
    
    survey_list.append('ATK')
    parameter_list.append('ATK pos')
    value_list.append(str(pos))
    error_list.append('---')
    note_list.append('ATK input pos')
    
    survey_list.append('ATK')
    parameter_list.append('ATK identifier')
    value_list.append(str(identifier))
    error_list.append('---')
    note_list.append('---')

    # combine lists into list_dict
    list_dict={'survey_list':survey_list,'parameter_list':parameter_list,'value_list':value_list,'error_list':error_list,'note_list':note_list}
    
    # get rows
    for key in metadata_dict:
        sub_dict=metadata_dict[key]
        if sub_dict=='default':
            if key not in supported_surveys:
                raise Exception('default metadata only supported for supported surveys.')
            else:
                list_dict=getdata(list_dict=list_dict,sub_dict=sub_dict,survey=key,source=source,pos=pos,default=True)
        else:
            list_dict=getdata(list_dict=list_dict,sub_dict=sub_dict,survey=key,source=source,pos=pos,default=False)

    # replace None (now 'None' as everything is converted to strings) to '---'
    survey_list=['---' if x=='None' else x for x in survey_list]
    parameter_list=['---' if x=='None' else x for x in parameter_list]
    value_list=['---' if x=='None' else x for x in value_list]
    error_list=['---' if x=='None' else x for x in error_list]
    note_list=['---' if x=='None' else x for x in note_list]

    # create table
    data=dict(
        survey=survey_list,
        parameter=parameter_list,
        value=value_list,
        error=error_list,
        notes=note_list
        )
    
    source=ColumnDataSource(data)
    columns=[
        TableColumn(field='survey',title='Survey'),  
        TableColumn(field='parameter',title='Parameter'),
        TableColumn(field='value',title='Value'),
        TableColumn(field='error',title='Error'),
        TableColumn(field='notes',title='Notes')
        ]
    
    data_table=DataTable(source=source,columns=columns,width=1200,height=400)

    return data_table