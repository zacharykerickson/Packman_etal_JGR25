import xarray
import numpy as np
import pandas as pd
import datetime
from glob import glob
from netCDF4 import Dataset

# Implementing decimal year fxn found here: https://stackoverflow.com/questions/6451655/how-to-convert-python-datetime-dates-to-decimal-float-years
# Input: array of np.datetime64 objects
# Output: same array converted into decimal year format, to be used for fitting models to a time series
def to_decimal_year(time):
    dt = pd.to_datetime(time) # converts np.datetime64 into datetime object
    dec_year = []
    for i in range(len(dt)):
        start = datetime.date(dt[i].year,1,1).toordinal()
        year_length = datetime.date(dt[i].year+1, 1, 1).toordinal() - start
        dec_year.append(dt[i].year + (float(dt[i].toordinal() - start) / year_length)) # calculates how much time has passed in the year, divides by year length, then adds to start year
    return np.array(dec_year)

# Remove trend from data by fitting a series of harmonics
# Input: array of times, data array, number of harmonics (int), whether or not to remove trend (Boolean)
# Output: data array with mean, trend and/or seasonal cycle removed if specified
def detrend_data(t, data, num_harmonics=0, remove_mean=False, remove_trend=False, nc=False):
    dyr = to_decimal_year(t) # convert time data from date-time to decimal year
    
    model = [np.ones(len(dyr)), dyr - np.mean(dyr)]
    for i in range(num_harmonics): # adding harmonics to model
        omega = 2*(i+1)*np.pi
        model.append(np.sin(omega*dyr))
        model.append(np.cos(omega*dyr))
    pmodel = np.linalg.pinv(np.array(model))
    
    # calculate coefficients / residuals (as specified by input)
    if nc == True:
        coefs = np.matmul(pmodel.transpose(), data)
    else: 
        coefs = np.matmul(data,pmodel)
    mean = coefs[0]
    trend = coefs[1]*(dyr - np.mean(dyr))
    res_coefs = coefs[2:]
    seasons_model = [] # initialize model to be populated with code below
    
    if remove_mean == True:
        detrend_data = data - mean    
    else: 
        detrend_data = data

    if remove_trend == True: # case for removing both trend and seasonal cycle
        detrend_data = detrend_data - trend
        for j in range(0,len(res_coefs),2):
            omega = (j+2)*np.pi
            detrend_data = detrend_data-(res_coefs[j]*np.sin(omega*dyr))-(res_coefs[j+1]*np.cos(omega*dyr)) 
            seasons_model.append(res_coefs[j]*np.sin(omega*dyr))
            seasons_model.append(res_coefs[j+1]*np.cos(omega*dyr))
            
    else: # case for leaving trend, removing seasonal cycle
        for j in range(0,len(res_coefs),2):
            omega = (j+2)*np.pi
            detrend_data = detrend_data-(res_coefs[j]*np.sin(omega*dyr))-(res_coefs[j+1]*np.cos(omega*dyr))
            seasons_model.append(res_coefs[j]*np.sin(omega*dyr))
            seasons_model.append(res_coefs[j+1]*np.cos(omega*dyr))
    
    seasons_model_final = np.sum(np.array(seasons_model), axis=0)
    
    results = {}
    results['data'] = detrend_data
    results['coefs'] = coefs
    results['trend'] = trend
    results['seasons_model'] = seasons_model_final
    
    return results


# Extract data
files = np.sort(glob('/Users/erickson/Documents/Data/RFROM/*.nc'))
data = xarray.open_mfdataset(files)
latitude = data['latitude'].values
longitude = data['longitude'].values
time = data['time'].values
pressure = data['mean_pressure'].values

output_nc = Dataset('trend.nc','w')
output_nc.createDimension('latitude',len(latitude))
output_nc.createDimension('longitude',len(longitude))
output_nc.createDimension('pressure',len(pressure))
output_nc.createDimension('coeffs',14)
output_nc.createVariable('latitude','f4',('latitude'))[:] = latitude
output_nc.createVariable('longitude','f4',('longitude'))[:] = longitude
coefs_var = output_nc.createVariable('coefs','f4',('coeffs','pressure','latitude','longitude'))

kinds = np.where(np.isin(data['mean_pressure'].values,[2.5,50,100,150,200,250]))[0]
for i in range(len(latitude)):
    print(latitude[i],end='')
    for j in range(len(longitude)):
        print('.',end='',flush=True)
        if np.isfinite(data['ocean_temperature'][0,0,i,j]):
            for k in kinds:
                temps = data['ocean_temperature'][:,k,i,j].values
                if np.all(np.isfinite(temps)):
                    temp_detrend = detrend_data(time,temps,num_harmonics=6,remove_trend=True)
                    coefs_var[:,k,i,j] = temp_detrend['coefs']
                else:
                    coefs_var[:,k,i,j] = np.nan
        else:
            coefs_var[:,:,i,j] = np.nan


output_nc.close()
data.close()
