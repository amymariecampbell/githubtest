# -*- coding: utf-8 -*-
"""
Created on Wed May  5 13:22:00 2021

@author: AC20
"""
import pandas as pd
import numpy as np

#Read MLST output csv
mlst=pd.read_csv("/Users/AC20/Documents/PhD/Year 1/South American and El Nino/Peru_Chile_Data/Genomes/mlst.csv")

#Read TSV output from NCBI isolates browser as a CSV
genomes=pd.read_csv("/Users/AC20/Documents/PhD/Year 1/South American and El Nino/Peru_Chile_Data/Genomes/isolates.tsv", sep="\t") 

#Merge them together based on the Assembly column 
data=pd.merge(genomes, mlst, on='Assembly')

#Save as CSV
data.to_csv('/Users/AC20/Documents/PhD/Year 1/South American and El Nino/Peru_Chile_Data/Genomes/data_beforeedit.csv')



#After openrefine, reopen data that has been editted
data=pd.read_csv("/Users/AC20/Documents/PhD/Year 1/South American and El Nino/Peru_Chile_Data/Genomes/data.csv")

ST=data.copy()
ST.set_index("Year", inplace=True)
ST=ST.drop(columns=['Column', '#Organism Group', 'Strain', 'Isolate identifiers', 'Serovar',
       'Isolate', 'Create date', 'Month', 'Day', 'Location', 'Lat/Lon',
       'Isolation Source', 'Isolation type', 'SNP cluster', 'Min-same',
       'Min-diff', 'BioSample', 'Assembly', 'AMR genotypes', 'Computed types',
       'Filename', 'Scheme', 'dnaE', 'gyrB', 'recA', 'dtdS', 'pntA',
       'pyrC', 'tnaA'])

#Extract the main 3 as their own columns- 3, 120, 36
ST['ST3']=np.where(ST['ST']== '3.0', '1', '0').astype(float)
ST['ST120']=np.where(ST['ST']== '120.0', '1', '0').astype(float)
ST['ST36']=np.where(ST['ST']== '36.0', '1', '0').astype(float)

ST=ST.drop(columns=['ST'])

ST=ST.groupby(["Year"]).sum()

ST.plot()

ST.reset_index()

#Climate data
import xarray
nc_sss=xarray.open_dataset("/Users/AC20/Documents/Data/Climate/sss_lowres.nc")
#Subset the netcdf files to the region and time period necessary 
sss1=nc_sss.loc[dict(lat=slice(-65, 0), lon=slice(-90, -50), time=slice('2010-01-01', '2018-12-31'))]
sss2=sss1['sss'].resample(time='1MS').mean(dim="time")
#As it is a dask array, need to load at the end 
sss=sss2.load()
sss.to_netcdf(path="/Users/AC20/Documents/Data/Climate/sss_Peru_Chile.nc")


sss_annual=sss.groupby("time.year").mean()
sss_annual=sss_annual.mean(dim=["lat", "lon"]).to_dataframe()

sss_annual=sss_annual.reset_index()
sss_annual.columns=["Year", "sss"]


sss_annual.plot()

#Make climate anomaly layer

climatology_mean = sss.groupby("time.year").mean("time")
climatology_std = sss.groupby("time.year").std("time")
stand_anomalies = xarray.apply_ufunc(
    lambda x, m, s: (x - m) / s,
    sss.groupby("time.year"),
    climatology_mean,
    climatology_std,
)
anomaly=stand_anomalies.mean(dim=['lat', 'lon']).to_dataframe()['sss']
annualanomaly=anomaly.groupby(anomaly.index.year).mean()
annualanomaly=anomaly.groupby(anomaly.index.year).mean()
annualanomaly.plot(figsize=(15,5))

annualanomaly=annualanomaly.reset_index()
annualanomaly.columns=["Year", "sss_anomaly"]




testdataset=ST.merge(sss_annual, on='Year')
testdataset=testdataset.merge(annualanomaly, on='Year')
testdataset.set_index('Year')

testdataset.plot()


import seaborn as sns
import matplotlib.pyplot as plt
pairplot=sns.pairplot(testdataset)
#sns.pairplot(df, hue='state')
#Save figure
#pairplot.savefig("/Users/amycampbell/Documents/YGT_Project/Results/Correlation_Analysis/coastal_ocean/pairplot")
plt.show()


fig, axs = plt.subplots(ncols=1,nrows=3, figsize=(5,20))
sns.regplot(x=testdataset ['ST3'], y=testdataset['sss_anomaly'], data=testdataset.loc[testdataset.index == "III"], scatter_kws={"s": 80},robust=True, ci=68, ax=axs[0])
sns.regplot(x=testdataset ['ST120'], y=testdataset['sss_anomaly'], data=testdataset.loc[testdataset.index == "III"], scatter_kws={"s": 80},robust=True, ci=68, ax=axs[1])
sns.regplot(x=testdataset ['ST36'], y=testdataset['sss_anomaly'], data=testdataset.loc[testdataset.index == "III"], scatter_kws={"s": 80},robust=True, ci=68, ax=axs[2])
#plt.savefig("/Users/amycampbell/Documents/YGT_Project/Results/Correlation_Analysis/Coastal_ocean/logx_regression")



fig, ax = plt.subplots() # Create the figure and axes object
# Plot the first x and y axes:
testdataset.plot(x = 'Year', y = 'sss_anomaly', ax = ax) 
# Plot the second x and y axes. By secondary_y = True a second y-axis is requested:
# (see https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html for details)
testdataset.plot(x = 'Year', y = 'ST36', ax = ax, secondary_y = True) 

