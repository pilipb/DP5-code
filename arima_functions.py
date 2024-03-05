import pandas as pd
from statsmodels.tsa.arima.model import ARIMA
from sklearn.metrics import mean_squared_error
import numpy as np
import matplotlib.pyplot as plt
import warnings

import requests
import pandas as pd
import time
from OSGridConverter import grid2latlong


# script for returning elevation from lat, long, based on open elevation data
# which in turn is based on SRTM
def get_elevation(lat, long):
    query = ('https://api.open-elevation.com/api/v1/lookup'
             f'?locations={lat},{long}')
    r = requests.get(query).json()  # json object, various ways you can extract value
    # one approach is to use pandas json functionality:
    try:
        elevation = pd.json_normalize(r, 'results')['elevation'].values[0]
    except KeyError as e:
        print(f"Unable to normalize json: {e} Likely too many requests. Sleeping for 1 second and trying again.")
        time.sleep(1)
        r = requests.get(query).json()
        elevation = pd.json_normalize(r, 'results')['elevation'].values[0]

    return elevation

# a function that checks pairs of points at a radius around the original lat,lon and returns the maximum slope
def calc_slope(lat, lon, radius):
    if isinstance(lat, str):
        lat = float(lat)
    if isinstance(lon, str):
        lon = float(lon)
    # generate a locus of points around the original lat, lon in the N,NE,E,SE,S,SW,W,NW directions
    # calculate the slope between the original point and the new point
    # return the maximum slope
    lat1 = [lat, lat + radius, lat + radius, lat + radius, lat, lat - radius, lat - radius, lat - radius]
    lon1 = [lon, lon, lon + radius, lon - radius, lon + radius, lon, lon - radius, lon + radius]

    lat2 = [lat + radius, lat + radius, lat, lat - radius, lat - radius, lat - radius, lat, lat + radius]
    lon2 = [lon, lon + radius, lon + radius, lon + radius, lon, lon - radius, lon - radius, lon - radius]

    n_points = len(lat1)

    max_slope = 0
    for i in range(n_points):
        
        upstream = get_elevation(lat1[i], lon1[i])
        downstream = get_elevation(lat2[i], lon2[i])
        dist = np.sqrt((lat1[i] - lat2[i])**2 + (lon1[i] - lon2[i])**2 ) * 111139
        slope = abs(downstream - upstream)/dist
        if slope > max_slope:
            max_slope = slope

    return max_slope

# Function to read the data from the file
def read_data(file_path):
    metadata = {}
    data = []

    with open(file_path, 'r', encoding='latin-1') as file:
        for line in file:
            if line.startswith('#'):
                if line.startswith('# DATA'):
                    continue
                elif line.startswith('# '):
                    line_split = line[2:].split(':')
                    if len(line_split) == 2:
                        key = line_split[0]
                        value = line_split[1]
                        metadata[key.strip()] = value.strip()
                    else:
                        continue
            elif line.startswith('YYYY-MM-DD'):
                continue
            else:
                data.append(line.strip().split(';'))

    return metadata, data

def read_data_from_csv(file_path):
    '''
    file in form:
    file,timestamp,2024-03-05T18:11:35
    database,id,nrfa-public-30
    database,name,UK National River Flow Archive
    station,id,53001
    station,name,Avon at Melksham
    station,gridReference,ST9029864102
    station,descriptionSummary,Velocity area station.
    station,descriptionGeneral,Velocity area station. Closed in 1980.
    station,descriptionStationHydrometry,Cableway present.
    station,descriptionCatchment,Gaugings completed in 1975/76 (upon which the recent Flood Studies high flows are based) suggest that archived dmfs appreciably underestimate the true discharge.
    dataType,id,gdf
    dataType,name,Gauged Daily Flow
    dataType,parameter,Flow
    dataType,units,m3/s
    dataType,period,day (P1D)
    dataType,measurementType,Mean
    data,first,1953-10-01
    data,last,1980-06-30
    1953-10-01,5.748
    1953-10-02,5.465
    1953-10-03,4.446
    1953-10-04,4.078
    
    
    '''


    metadata = {}
    data = []

    with open(file_path, 'r', encoding='latin-1') as file:
        for line in file:
            if line.startswith('file'):
                continue
            elif line.startswith('database'):
                continue
            elif line.startswith('station'):
                data.append(line.strip().split(','))
                continue
            elif line.startswith('dataType'):
                continue
            elif line.startswith('data'):
                continue
            else:
                data.append(line.strip().split(','))
    # print(data)
    metadata['Grid Reference'] = str(data[2][2])
    metadata['River'] = data[1][2]
    l = grid2latlong(metadata['Grid Reference'])
    metadata['Latitude (DD)'] = l.latitude
    metadata['Longitude (DD)'] = l.longitude

    return metadata, data[10:]

def process_data_from_csv(metadata, data):
    # Extract metadata
    station_info = {
        'River': metadata.get('River'),
        'Grid Reference': metadata.get('Grid Reference'),
        'Latitude (DD)': metadata.get('Latitude (DD)'),
        'Longitude (DD)': metadata.get('Longitude (DD)'),
    }

    # Extract data
    parsed_data = [{'Date': row[0], 'Value': float(row[1]) if row[1] != '--:--' else None} for row in data]
    parsed_data = pd.DataFrame(parsed_data)
    parsed_data['Value'] = parsed_data['Value'].replace(-999.0, np.nan)
    return station_info, parsed_data


# Function to process the data
def process_data(metadata, data):
    # Extract metadata
    station_info = {
        'GRDC-No.': metadata.get('GRDC-No.'),
        'River': metadata.get('River'),
        'Station': metadata.get('Station'),
        'Country': metadata.get('Country'),
        'Latitude (DD)': metadata.get('Latitude (DD)'),
        'Longitude (DD)': metadata.get('Longitude (DD)'),
        'Catchment area (km²)': metadata.get('Catchment area (km�)'),
        'Altitude (m ASL)': metadata.get('Altitude (m ASL)'),
        'Time series': metadata.get('Time series'),
        'No. of years': metadata.get('No. of years'),
        'Last update': metadata.get('Last update')
    }

    # Extract data
    parsed_data = [{'Date': row[0], 'Value': float(row[2]) if row[2] != '--:--' else None} for row in data]
    parsed_data = pd.DataFrame(parsed_data)
    parsed_data['Value'] = parsed_data['Value'].replace(-999.0, np.nan)
    return station_info, parsed_data

def downsample(data):
    data['Date'] = pd.to_datetime(data['Date'])
    data = data.set_index('Date')

    # if there is a Velocity column, sample it
    if 'Velocity' in data.columns:
        data['Sampled'] = data['Velocity'].resample('5D').mean()
    else: 
        data['Sampled'] = data['Value'].resample('5D').mean()

    return data


def split_data(data, split=0.8):

    train = data['Sampled'].iloc[:int(len(data)*split)]
    test = data['Sampled'].iloc[int(len(data)*split):]
    return train, test

def arima_model(train, order=(2,0,2), seasonal_order=(0,1,0,73)):
    warnings.filterwarnings("ignore")
    model = ARIMA(train.dropna(), order=order, seasonal_order=seasonal_order)
    model_fit = model.fit()
    # print(model_fit.summary())
    return model_fit

def forecast(model, test):
    forecast = model.forecast(steps=int(len(test)/5))
    return forecast

def evaluate_model(model, test, plot=True):
    # Get the forecasted values
    # mute warnings
    warnings.filterwarnings("ignore")
    forecast_values = model.forecast(steps=int(len(test)/5))
    # Drop NaN values
    test_clean = test.dropna()
    forecast_values_clean = forecast_values.dropna()

    # Truncate the longer series to the length of the shorter series
    if len(test_clean) > len(forecast_values_clean):
        test_clean = test_clean[:len(forecast_values_clean)]
    elif len(forecast_values_clean) > len(test_clean):
        forecast_values_clean = forecast_values_clean[:len(test_clean)]

    # Compute the mean squared error
    mse = mean_squared_error(test_clean, forecast_values_clean)
    # print("Mean Squared Error:", mse)
    rmse = np.sqrt(mse)
    # print("Root Mean Squared Error:", rmse)
    mape = np.mean(np.abs((test.dropna() - forecast_values) / test.dropna())) * 100
    # print("Mean Absolute Percentage Error:", mape)

    # plot the residuals
    if plot:
        residuals = model.resid
        plt.figure(figsize=(10,6))
        plt.plot(residuals)
        plt.title('Residuals')
        plt.show()

    return mse, rmse, mape