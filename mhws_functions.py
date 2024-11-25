import pandas as pd
import datetime
import numpy as np

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

# Performs lagged correlation analysis on two time series using np.correlate
# Inputs: two time series of equal length (arrays), number of years to search thru for a max correlation, whether to normalize (Boolean)
# Outputs: array of all possible lags + correlation, max correlation, and lag @ which max correlation occurs
    # Interpreting output: max lag is the amount of time the first time series must be shifted forward in order to "catch up" with second time series
def lag_corr(t1, t2, search_years, norm = True):
    if norm == True: 
        a = (t2 - np.mean(t2)) / (np.std(t2)*len(t2)) # Normalizing time series
        b = (t1 - np.mean(t1)) / np.std(t1)
    else:
        a = t2 / len(t2)
        b = t1
    
    corr = np.correlate(a, b, 'full') 
    lags = np.arange(-len(a)+1, len(a)) # all possible lags over which np.correlate() operates
    
    # Max correlation (lead/lag up to specified number of years)
    weeks = search_years * 52
    median_index = np.where(lags == np.median(lags))[0][0]
    max_corr = max(corr[median_index - weeks : median_index + weeks])
    
    # Accounting for possible exceptions (if data doesn't exist at certain point)
    if np.isnan(max_corr)==False:
    	max_lag = float(lags[np.where(corr == max_corr)[0][0]])
    else:
        max_lag = float('NaN')
    
    # Organizing results
    results = {}
    results['lags'] = lags
    results['corr'] = corr
    results['max_lag'] = max_lag
    results['max_corr'] = max_corr
    
    return results