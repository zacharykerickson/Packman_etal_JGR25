import xarray
import numpy as np
from glob import glob
from mhws_functions import detrend_data
from netCDF4 import Dataset

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
