from datetime import datetime
import pandas as pd
import numpy as np

def filter_period(df, period):
    period_0 = str(datetime(period[0][0], period[0][1], period[0][2], 0, 0).date())
    period_1 = str(datetime(period[1][0], period[1][1], period[1][2], 0, 0).date())
    df = df.loc[(df.time >= period_0)&(df.btime <= period_1)]
    return df

def filter_extent(df, extent):
    min_x, max_x = extent[0], extent[1]
    min_y, max_y = extent[2], extent[3]
    
    if max_x<=180:
        df = df.loc[(df.bgn_x > min_x) & (df.bgn_x < max_x)]
        df = df.loc[(df.bgn_y > min_y) & (df.bgn_y < max_y)]
    elif max_x>180: 
        max_x_ = max_x -360
        df = df.loc[((df.bgn_x > min_x) & (df.bgn_x < 180)) | ((df.bgn_x > -180)&(df.bgn_x < max_x_))]
        df = df.loc[(df.bgn_y > min_y) & (df.bgn_y < max_y)]
        df = df.loc[((df.bgn_x <0 ) & (df.end_x >0)) | ((df.bgn_x >0) & (df.end_x <0 ))]
    return df

def filter_presure(df, pres):
    df = df.loc[(df.pres > pres[0])&(df.pres < pres[1])]
    return df

def filter_wmos(df, wmos):
    wmo = wmos[0]
    df = df.loc[(df.wmo==wmo)]
    return df



if __name__ == '__main__':

    extent = [-78, -60, 32, 47]  # lonmin, lonmax, latmin, latmax  
    period = [(2017, 1, 1), (2017, 12, 31)] 
    pressure = [950, 1050]
    wmos = [1901192]

    df = pd.read_csv('ANDRO.csv')
    df = filter_period(df, period) 
    df = filter_extent(df, extent)
    df = filter_presure(df, pressure) 
    df = filter_wmos(df, wmos)

    df.to_csv('data_out', index=False)
