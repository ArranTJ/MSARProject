# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 11:09:46 2024

@author: Arran
"""

import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import time_support

from sunpy import timeseries as ts
from sunpy.net import Fido
from sunpy.net import attrs as a
import pandas as pd
import matplotlib.dates as md


import csv
import glob




tstart = "2017-09-06 07:00"
tend = "2017-09-06 15:00"
result = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"))
print(result)


result_goes16 = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"), a.goes.SatelliteNumber(16), a.Resolution("flx1s"))
print(result_goes16)


file_goes16 = Fido.fetch(result_goes16)

goes_16 = ts.TimeSeries(file_goes16)
# goes_16.peek()

fig, ax = plt.subplots(figsize=[8,5])
plt.title("GOES-16 X-ray Flux")
goes_16.plot()
# goes_16.plot(axes=ax, columns=["xrsb"])
# ax.set_yscale('log')

plt.savefig("GOES16Flux.pdf", dpi = 300)


rows = []
for filename in glob.glob("Kappa1\goodkappavals\*.csv"):
   with open(filename, newline="") as f:
      reader = csv.reader(f)
      rows.extend(reader)

dframe = pd.DataFrame(rows)


dframe = dframe.set_index(dframe[10])


fig, axs = plt.subplots(4, 2, figsize=[12,8])
df = goes_16.to_dataframe()
df = df[(df["xrsa_quality"] == 0) & (df["xrsb_quality"] == 0)]
goes_16 = ts.TimeSeries(df, goes_16.meta, goes_16.units)
# plt.tight_layout()
fig.subplots_adjust(hspace=0.05, wspace=0)
plt.suptitle("GOES-16 Comparison of Kappa Fits, \n for Kappa Function described in 'Jeffrey. N. et al (2016)'")
# axs[0] = plt.subplot(411)
axs[0][0] = goes_16.plot(axes=axs[0][0], columns=["xrsb"], color = 'r')
axs[0][0].set_ylabel('GOES-16, 1-8 $\AA$ \n W $m^{-2}$')
axs[0][0].set_title('Fe XVI GOES-16 Comparison')
axs[0][0].set_yscale('log')
axs[0][0].set_xlim(pd.Timestamp('2017-09-06 8:00:00'), pd.Timestamp('2017-09-06 14:50:00'))
axs[0][0].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
####Annotating regions
axs[0][0].plot()
axs[0][0].annotate('T1', xy = (pd.Timestamp('2017-09-06 8:30:00'), 1e-3))
axs[0][0].annotate('T2', xy = (pd.Timestamp('2017-09-06 9:05:00'), 1e-3))
axs[0][0].annotate('T3', xy = (pd.Timestamp('2017-09-06 10:00:00'), 1e-3))
axs[0][0].annotate('T4', xy = (pd.Timestamp('2017-09-06 11:20:00'), 1e-3))
axs[0][0].annotate('T5', xy = (pd.Timestamp('2017-09-06 12:10:00'), 1e-3))
axs[0][0].annotate('T6', xy = (pd.Timestamp('2017-09-06 13:00:00'), 1e-3))
plt.legend()

####Plotting Vertical Lines
axs[0][0].axvline(pd.Timestamp('2017-09-06 8:58:00'), color='k', ls=':')
axs[0][0].axvline(pd.Timestamp('2017-09-06 9:25:00'), color='k', ls=':')
axs[0][0].axvline(x=pd.Timestamp('2017-09-06 11:00:00'), color='k', ls=':')
axs[0][0].axvline(x=pd.Timestamp('2017-09-06 11:55:00'), color='k', ls=':')
axs[0][0].axvline(x=pd.Timestamp('2017-09-06 12:25:00'), color='k', ls=':')

plt.setp(axs[0][0].get_xticklabels(), visible=False)

axs[0][1] = goes_16.plot(axes=axs[0][1], columns=["xrsb"], color = 'r')
axs[0][1].label_outer()
axs[0][1].set_title('Fe XXIII GOES-16 Comparison')
axs[0][1].set_yscale('log')
axs[0][1].set_xlim(pd.Timestamp('2017-09-06 8:00:00'), pd.Timestamp('2017-09-06 14:50:00'))
axs[0][1].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
####Plotting Vertical Lines
axs[0][1].axvline(pd.Timestamp('2017-09-06 8:58:00'), color='k', ls=':')
axs[0][1].axvline(pd.Timestamp('2017-09-06 9:25:00'), color='k', ls=':')
axs[0][1].axvline(x=pd.Timestamp('2017-09-06 11:00:00'), color='k', ls=':')
axs[0][1].axvline(x=pd.Timestamp('2017-09-06 11:55:00'), color='k', ls=':')
axs[0][1].axvline(x=pd.Timestamp('2017-09-06 12:25:00'), color='k', ls=':')

####Annotating regions
axs[0][1].plot()
axs[0][1].annotate('T1', xy = (pd.Timestamp('2017-09-06 8:30:00'), 1e-3))
axs[0][1].annotate('T2', xy = (pd.Timestamp('2017-09-06 9:05:00'), 1e-3))
axs[0][1].annotate('T3', xy = (pd.Timestamp('2017-09-06 10:00:00'), 1e-3))
axs[0][1].annotate('T4', xy = (pd.Timestamp('2017-09-06 11:20:00'), 1e-3))
axs[0][1].annotate('T5', xy = (pd.Timestamp('2017-09-06 12:10:00'), 1e-3))
axs[0][1].annotate('T6', xy = (pd.Timestamp('2017-09-06 13:00:00'), 1e-3))
plt.legend()



###NOTE, all values in the dataframe imported from .csv exist as strings, 
###therefore, they need to be converted to floats, or 'time series' types when necessary

# axs[1] = plt.subplot(412)
axs[1][0].plot(pd.to_datetime(dframe.index), np.array(dframe[0]).astype(float), color='b')
axs[1][0].set_ylabel('Number of \n Good Kappa Fits')
axs[1][0].set_xlim(pd.Timestamp('8:00:00'), pd.Timestamp('14:50:00'))
axs[1][0].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))

plt.setp(axs[1][0].get_xticklabels(), visible=False)
# plt.setp(ax3.get_xticklabels(), visible=False)
####Plotting Vertical Lines
axs[1][0].axvline(pd.Timestamp('8:58:00'), color='k', ls=':')
axs[1][0].axvline(pd.Timestamp('9:25:00'), color='k', ls=':')
axs[1][0].axvline(x=pd.Timestamp('11:00:00'), color='k', ls=':')
axs[1][0].axvline(x=pd.Timestamp('11:55:00'), color='k', ls=':')
axs[1][0].axvline(x=pd.Timestamp('12:25:00'), color='k', ls=':')



axs[1][1].sharey(axs[1][0])
axs[1][1].label_outer()
axs[1][1].plot(pd.to_datetime(dframe.index), np.array(dframe[1]).astype(float), color='m')
# axs[1][1].set_ylabel('Number of \n Good Kappa Fits')
axs[1][1].set_xlim(pd.Timestamp('8:00:00'), pd.Timestamp('14:50:00'))
axs[1][1].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))

####Plotting Vertical Lines
axs[1][1].axvline(pd.Timestamp('8:58:00'), color='k', ls=':')
axs[1][1].axvline(pd.Timestamp('9:25:00'), color='k', ls=':')
axs[1][1].axvline(x=pd.Timestamp('11:00:00'), color='k', ls=':')
axs[1][1].axvline(x=pd.Timestamp('11:55:00'), color='k', ls=':')
axs[1][1].axvline(x=pd.Timestamp('12:25:00'), color='k', ls=':')


axs[2][0].sharex(axs[1][0])
axs[2][0].errorbar(pd.to_datetime(dframe.index), np.array(dframe[2]).astype(float), yerr =np.array(dframe[3]).astype(float), color='b', fmt='.')
# ax2.errorbar(pd.to_datetime(df.index), np.array(df[4]).astype(float), yerr =np.array(df[5]).astype(float), color='m',fmt='.')
axs[2][0].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
axs[2][0].set_ylabel('Mean $\kappa$')
axs[2][0].set_ylim(1,12.5)
plt.setp(axs[2][0].get_xticklabels(), visible=False)
####Plotting Vertical Lines
axs[2][0].axvline(pd.Timestamp('8:58:00'), color='k', ls=':')
axs[2][0].axvline(pd.Timestamp('9:25:00'), color='k', ls=':')
axs[2][0].axvline(x=pd.Timestamp('11:00:00'), color='k', ls=':')
axs[2][0].axvline(x=pd.Timestamp('11:55:00'), color='k', ls=':')
axs[2][0].axvline(x=pd.Timestamp('12:25:00'), color='k', ls=':')



axs[2][1].sharex(axs[1][1])
axs[2][1].sharey(axs[2][0])
axs[2][1].errorbar(pd.to_datetime(dframe.index), np.array(dframe[4]).astype(float), yerr=np.array(dframe[5]).astype(float), color='m',fmt='.')
axs[2][1].label_outer()
# ax2.errorbar(pd.to_datetime(df.index), np.array(df[4]).astype(float), yerr =np.array(df[5]).astype(float), color='m',fmt='.')
axs[2][1].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
# axs[2][1].set_ylabel('Mean $\kappa$')
axs[2][1].set_ylim(1,12.5)
####Plotting Vertical Lines
axs[2][1].axvline(pd.Timestamp('8:58:00'), color='k', ls=':')
axs[2][1].axvline(pd.Timestamp('9:25:00'), color='k', ls=':')
axs[2][1].axvline(x=pd.Timestamp('11:00:00'), color='k', ls=':')
axs[2][1].axvline(x=pd.Timestamp('11:55:00'), color='k', ls=':')
axs[2][1].axvline(x=pd.Timestamp('12:25:00'), color='k', ls=':')


axs[3][0].sharex(axs[1][0])
axs[3][0].errorbar(pd.to_datetime(dframe.index), np.array(dframe[6]).astype(float), yerr =np.array(dframe[7]).astype(float), color='b', fmt='.')
# ax2.errorbar(pd.to_datetime(df.index), np.array(df[4]).astype(float), yerr =np.array(df[5]).astype(float), color='m',fmt='.')
axs[3][0].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
axs[3][0].set_ylabel('Mean $\chi^{2}$')
axs[3][0].set_ylim(0.05)
# plt.setp(axs[3][0].get_xticklabels(), visible=False)
####Plotting Vertical Lines
axs[3][0].axvline(pd.Timestamp('8:58:00'), color='k', ls=':')
axs[3][0].axvline(pd.Timestamp('9:25:00'), color='k', ls=':')
axs[3][0].axvline(x=pd.Timestamp('11:00:00'), color='k', ls=':')
axs[3][0].axvline(x=pd.Timestamp('11:55:00'), color='k', ls=':')
axs[3][0].axvline(x=pd.Timestamp('12:25:00'), color='k', ls=':')


axs[3][1].sharex(axs[1][1])
axs[3][1].sharey(axs[3][0])
axs[3][1].errorbar(pd.to_datetime(dframe.index), np.array(dframe[8]).astype(float), yerr =np.array(dframe[9]).astype(float), color='m',fmt='.')
axs[3][1].label_outer()
# ax2.errorbar(pd.to_datetime(df.index), np.array(df[4]).astype(float), yerr =np.array(df[5]).astype(float), color='m',fmt='.')
axs[3][1].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
# axs[2][1].set_ylabel('Mean $\kappa$')
axs[3][1].set_ylim(0.05)
####Plotting Vertical Lines
axs[3][1].axvline(pd.Timestamp('8:58:00'), color='k', ls=':')
axs[3][1].axvline(pd.Timestamp('9:25:00'), color='k', ls=':')
axs[3][1].axvline(x=pd.Timestamp('11:00:00'), color='k', ls=':')
axs[3][1].axvline(x=pd.Timestamp('11:55:00'), color='k', ls=':')
axs[3][1].axvline(x=pd.Timestamp('12:25:00'), color='k', ls=':')



plt.savefig("GOESKappa1DataComparison.pdf", bbox_inches='tight')

#Function to get stats from data

def datastats(df):
    #Eliminating NaNs:
    df1 = np.array(df[0]).astype(float)[~np.isnan(np.array(df[0]).astype(float))]
    df2 = np.array(df[1]).astype(float)[~np.isnan(np.array(df[1]).astype(float))]
    df3 = np.array(df[2]).astype(float)[~np.isnan(np.array(df[2]).astype(float))]
    df4 = np.array(df[3]).astype(float)[~np.isnan(np.array(df[3]).astype(float))]
    df5 = np.array(df[4]).astype(float)[~np.isnan(np.array(df[4]).astype(float))]
    df6 = np.array(df[5]).astype(float)[~np.isnan(np.array(df[5]).astype(float))]
    df7 = np.array(df[6]).astype(float)[~np.isnan(np.array(df[6]).astype(float))]
    df8 = np.array(df[7]).astype(float)[~np.isnan(np.array(df[7]).astype(float))]
    df9 = np.array(df[8]).astype(float)[~np.isnan(np.array(df[8]).astype(float))]
    df10 = np.array(df[9]).astype(float)[~np.isnan(np.array(df[9]).astype(float))]
    
    #Eliminating Zeros (Exceptions are df1 and df2 where 0 counts for kappa makes sense):
    df3 = df3[df3 != 0 ]
    df4 = df4[df4 != 0 ]
    df5 = df5[df5 != 0 ]
    df6 = df6[df6 != 0 ]
    df7 = df7[df7 != 0 ]
    df8 = df8[df8 != 0 ]
    df9 = df9[df9 != 0 ]
    df10 = df10[df10!= 0 ]
    #Total Counts for 16
    Total16 = np.sum(df1)
    #Mean Count for 16
    MeanCt16 = np.mean(df1)
    #Std Dev for MeanCt16
    StdDevCt16 = np.std(df1)
    #Total Counts for 23
    Total23 = np.sum(df2)
    #Mean Count for 23
    MeanCt23 = np.mean(df2)
    #Std Dev for MeanCt23
    StdDevCt23 = np.std(df2)
    
    #MeanKappa16
    MeanK16 = np.mean(df3)
    #Error Propagation for Kappa16
    ErrK16 = np.sqrt(np.sum((df4)**2))/len(df4)
    
    #MeanKappa23
    MeanK23 = np.mean(df5)
    #Error Propagation for Kappa23
    ErrK23 = np.sqrt(np.sum((df6)**2))/len(df6)
    
    #MeanChi16
    MeanChi16 = np.mean(df7)
    #Error Propagation for Chi16
    ErrChi16 = np.sqrt(np.sum((df8)**2))/len(df8)
    #MeanChi23
    MeanChi23 = np.mean(df9)
    #Error Propagation for Chi23
    ErrChi23 = np.sqrt(np.sum((df10)**2))/len(df10)
    
    return Total16, MeanCt16, StdDevCt16, Total23, MeanCt23, StdDevCt23, MeanK16, ErrK16, MeanK23, ErrK23, MeanChi16, ErrChi16, MeanChi23, ErrChi23



#Define intervals of data for each phase of flare activity
df1 = dframe.loc['08:00:40':'08:57:35']
df2 = dframe.loc['09:00:25':'09:23:11']
df3 = dframe.loc['09:26:02':'10:59:57']
df4 = dframe.loc['11:02:48':'11:54:01']
df5 = dframe.loc['11:56:51':'12:22:28']
df6 = dframe.loc['12:25:19':'14:50:27']
#Obtain statistics from each interval
statT1 = datastats(df1)
statT2 = datastats(df2)
statT3 = datastats(df3)
statT4 = datastats(df4)
statT5 = datastats(df5)
statT6 = datastats(df6)


###Plotting Standalone Figure for X-ray curve
fig, ax = plt.subplots(figsize=[10,5])

ax = goes_16.plot(axes=ax, columns=["xrsb"], color = 'r')
ax.set_yscale('log')
ax.set_title('GOES-16, 1-8 $\AA$ X-ray Flux')
ax.set_ylabel('W $m^{-2}$')
ax.set_xlim(pd.Timestamp('2017-09-06 8:00:00'), pd.Timestamp('2017-09-06 14:50:00'))
plt.axvline(x=pd.Timestamp('2017-09-06 8:58:00'), color='k', ls=':')
plt.axvline(x=pd.Timestamp('2017-09-06 9:25:00'), color='k', ls=':')
plt.axvline(x=pd.Timestamp('2017-09-06 11:00:00'), color='k', ls=':')
plt.axvline(x=pd.Timestamp('2017-09-06 11:55:00'), color='k', ls=':')
plt.axvline(x=pd.Timestamp('2017-09-06 12:25:00'), color='k', ls=':')

####Annotating regions
ax.plot()
ax.annotate('T1', xy = (pd.Timestamp('2017-09-06 8:30:00'), 1e-3))
ax.annotate('T2', xy = (pd.Timestamp('2017-09-06 9:05:00'), 1e-3))
ax.annotate('T3', xy = (pd.Timestamp('2017-09-06 10:00:00'), 1e-3))
ax.annotate('T4', xy = (pd.Timestamp('2017-09-06 11:20:00'), 1e-3))
ax.annotate('T5', xy = (pd.Timestamp('2017-09-06 12:10:00'), 1e-3))
ax.annotate('T6', xy = (pd.Timestamp('2017-09-06 13:00:00'), 1e-3))
plt.legend()

plt.savefig("GOES16FluxLabelBounds.pdf", bbox_inches='tight')


plt.show()