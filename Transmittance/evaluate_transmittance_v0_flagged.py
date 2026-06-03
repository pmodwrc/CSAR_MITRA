### evaluate_transmittance_v0_flagged.py: Evaluates transmittance of MITRA window from the flagged data file 
### Author: Natalia Engler 
### Copyright: 05.2026, MeteoSwiss, PMOD/WRC Davos

import numpy as np
import time
import matplotlib.pyplot as plt
import math
import matplotlib.dates as md
import datetime
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# Read a flagged data file
filename = '20260407_103318_MITRA_Data_flagged.txt'
PATH = 'MITRA Data/'
init_tstamp = time.mktime(datetime.datetime.strptime(filename[0:15], "%Y%m%d_%H%M%S").timetuple())
MITRA_df = pd.read_table(PATH+filename,  parse_dates=['Time'])

# Calculate detector resistances 
# R_cav = Sensing detector resistance 
# R_ref = Reference detector resistance
# R_dark = Dark detector resistance
length_data = len(MITRA_df['Time'])
R_cav = np.zeros(length_data-2)
R_ref = np.zeros(length_data-2)
R_dark = np.zeros(length_data-2)

for i in range(1, length_data-1):
    R_cav[i-1] = 0.5*(0.5*(abs(MITRA_df['V_cav'][i-1]/MITRA_df['Current'][i-1]) + abs(MITRA_df['V_cav'][i]/MITRA_df['Current'][i]) ) + 
                      0.5*(abs(MITRA_df['V_cav'][i+1]/MITRA_df['Current'][i+1]) + abs(MITRA_df['V_cav'][i]/MITRA_df['Current'][i]) ))
    R_ref[i-1] = 0.5*(0.5*(abs(MITRA_df['V_ref'][i-1]/MITRA_df['Current'][i-1]) + abs(MITRA_df['V_ref'][i]/MITRA_df['Current'][i]) ) + 
                      0.5*(abs(MITRA_df['V_ref'][i+1]/MITRA_df['Current'][i+1]) + abs(MITRA_df['V_ref'][i]/MITRA_df['Current'][i]) ))
    R_dark[i-1] = 0.5*(0.5*(abs(MITRA_df['V_dark'][i-1]/MITRA_df['Current'][i-1]) + abs(MITRA_df['V_dark'][i]/MITRA_df['Current'][i]) ) + 
                       0.5*(abs(MITRA_df['V_dark'][i+1]/MITRA_df['Current'][i+1]) + abs(MITRA_df['V_dark'][i]/MITRA_df['Current'][i]) ))

# Plot resistances:
fig, (ax1) = plt.subplots(1, figsize=(13,7))
ax1.plot((MITRA_df['Time'][1:length_data-1]), R_cav, c='g', label='Sensing detector')
ax1.plot((MITRA_df['Time'][1:length_data-1]), R_ref, c='r', label='Reference detector')
ax1.plot((MITRA_df['Time'][1:length_data-1]), R_dark, c='k', label='Dark detector') 
ax1.set_xlabel('Time (CET)', fontsize=18)
ax1.set_ylabel('Resistance (Ohm)', fontsize=18)
plt.title(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4] + ': Resistances', fontsize=18)
plt.legend(loc="best", fancybox=True, prop={'size':15})
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.gcf().autofmt_xdate()
myFmt = md.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
fig.autofmt_xdate(bottom=0.2, rotation=0, ha='center')
#plt.savefig(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4]+'_resistances.jpeg', bbox_inches="tight")
plt.show()


### Coefficients from Wolfgang 07.04.2026 to calculate detector temperature:
#Sensor 1: R0=314.11197, TK=4.7901941e-3
#Sensor 2: R0=316.09389, TK=4.9585408e-3
#Sensor 3: R0=310.32896, TK=4.8034457e-3
T0 = 9.9
TK_ref, R0_ref  =  316.09389*0.0049585408, 316.09389   # Sensor ID 21-003  
TK_cav, R0_cav  =  314.11197*0.0047901941, 314.11197   # Sensor ID 21-022  
TK_dark, R0_dark = 310.32896*0.0048034457, 310.32896   # Dark


# Calculate detector temperatures and temperature differences
# T_cav = Sensing detector temperature 
# T_ref = Reference detector temperature
# T_dark = Dark detector temperature
T_cav = (R_cav/R0_cav - 1)*R0_cav/TK_cav + T0 
T_ref = (R_ref/R0_ref - 1)*R0_ref/TK_ref + T0       
T_dark = (R_dark/R0_dark - 1)*R0_dark/TK_dark + T0

# Plot temperatures:
fig, (ax1) = plt.subplots(1, figsize=(13,7))
ax1.plot((MITRA_df['Time'][1:length_data-1]), T_cav, c='g', label='Sensing detector')
ax1.plot((MITRA_df['Time'][1:length_data-1]), T_ref, c='r', label='Reference detector')
ax1.plot((MITRA_df['Time'][1:length_data-1]), T_dark, c='k', label='Dark detector') 
ax1.set_xlabel('Time (CET)', fontsize=18)
ax1.set_ylabel('Temperature (°C)', fontsize=18)
plt.title(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4]+ ': Temperatures', fontsize=18 )
plt.legend(loc="best", fancybox=True, prop={'size':15})
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.gcf().autofmt_xdate()
myFmt = md.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
fig.autofmt_xdate(bottom=0.2, rotation=0, ha='center')
#plt.savefig(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4]+'_temperatures.jpeg', bbox_inches="tight")
plt.show()

# Plot temperature differences:
fig, (ax1) = plt.subplots(1, figsize=(13,7))
ax1.plot((MITRA_df['Time'][1:length_data-1]), T_cav-T_ref, c='k', label='T_cav - T_ref')
ax1.plot((MITRA_df['Time'][1:length_data-1]), T_cav-T_dark, c='g', label='T_cav - T_dark') 
ax1.plot((MITRA_df['Time'][1:length_data-1]), T_ref-T_dark, c='r', label='T_ref - T_dark') 
ax1.set_xlabel('Time (CET)', fontsize=18)
ax1.set_ylabel('Temperature (°C)', fontsize=18)
plt.title(f'{filename[6:8]}.{filename[4:6]}.{filename[0:4]}: Temperature differences', fontsize=18)
plt.legend(loc="best", fancybox=True, prop={'size':15})
plt.gcf().autofmt_xdate()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
myFmt = md.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
plt.grid(True)
fig.autofmt_xdate(bottom=0.2, rotation=0, ha='center')
#plt.savefig(f'{filename[6:8]}.{filename[4:6]}.{filename[0:4]}_temperature_diff.jpeg', bbox_inches="tight")
plt.show()

# Plot detector ratio:
fig, (ax1) = plt.subplots(1, figsize=(13,7))
ax1.plot((MITRA_df['Time'][1:length_data-1]), (T_cav - T_dark)/(T_ref - T_dark), c='b', label='(T_cav - T_dark)/(T_ref - T_dark)')
ax1.set_xlabel('Time (CET)', fontsize=18)
ax1.set_ylabel('Ratio', fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4] + ': Detector ratio', fontsize=18)
plt.legend(loc="best", fancybox=True, prop={'size':15})
plt.gcf().autofmt_xdate()
myFmt = md.DateFormatter('%H:%M')
#plt.ylim([0.85, 1.03])
plt.gca().xaxis.set_major_formatter(myFmt)
fig.autofmt_xdate(bottom=0.2, rotation=0, ha='center')
#plt.savefig(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4] +'_Detector_ratio.jpeg', bbox_inches="tight")
plt.show()


# Define time interval for data selection:
mask_time1 =  filename[0:4]+ '-'+ filename[4:6]+ '-'+ filename[6:8]+ ' 08:00:00'
mask_time2 = filename[0:4]+ '-'+ filename[4:6]+ '-'+ filename[6:8]+' 17:50:00'

# Create a new dataframe with data within the selected time interval:
new_df = MITRA_df[['Time', 'Win_position']].iloc[:-1]
new_MITRA_df = new_df.iloc[1:]
new_MITRA_df['ts'] = new_MITRA_df.Time.astype('int64') // 10**9 
new_MITRA_df['Ratio'] = ((T_cav - T_dark)/(T_ref - T_dark)).tolist()
new_MITRA_df['R_ref'] = R_ref.tolist()
new_MITRA_df['R_cav'] = R_cav.tolist()
new_MITRA_df['R_dark'] = R_dark.tolist()
mask = (new_MITRA_df['Time'] >= mask_time1) & (new_MITRA_df['Time'] <= mask_time2)
selection_df = new_MITRA_df.loc[mask] 

### Plot selected resistances:
fig, (ax1) = plt.subplots(1, figsize=(13,7))
ax1.plot(selection_df['Time'], selection_df['R_cav'], c='g', label='Sensing detector')
ax1.plot(selection_df['Time'], selection_df['R_ref'], c='r', label='Reference detector')
ax1.plot(selection_df['Time'], selection_df['R_dark'], c='k', label='Dark detector')
ax1.set_xlabel('Time (CET)', fontsize=18)
ax1.set_ylabel('R', fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title('Fig.5: Resistances, '+ filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4])
plt.legend(loc="best", fancybox=True, prop={'size':15})
plt.gcf().autofmt_xdate()
myFmt = md.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
fig.autofmt_xdate(bottom=0.2, rotation=0, ha='center')
plt.show()


#####################################################################################################################################################
#               Calculate MITRA transmittance in intervals; Fit the ratio in open and closed detector states
#######################################################################################################################################################
# Define data recorded in open state:
max_ratio = 0.96
min_ratio = 0.925

# Split the data into categories
start_state = selection_df['Win_position'].iloc[0]
if start_state == 'Right Cavity Closed':
    end_state = 'Both Open'
else:
    end_state = 'Right Cavity Closed'
    
# cycle start
start_of_cycle = selection_df['Win_position'].eq(start_state) & selection_df['Win_position'].shift().eq(end_state)

# if the dataframe itself starts with start_state, that is cycle 1
start_of_cycle.iloc[0] = selection_df['Win_position'].iloc[0] == start_state

# cumulative count of cycle starts
selection_df['cycle'] = start_of_cycle.cumsum()
mask_open = (selection_df['Time'] >= mask_time1) & (selection_df['Time'] <= mask_time2) & (selection_df['Ratio'] >= max_ratio) \
                & (selection_df['Ratio'] < 1)  & (selection_df['Win_position'] == 'Both Open')

# Define data recorded in closed state:
mask_closed = (selection_df['Time'] >= mask_time1) & (selection_df['Time'] <= mask_time2) & (selection_df['Ratio'] <= min_ratio) \
                & (selection_df['Win_position'] != 'Both Open') 

# Print mean value and standard deviation of datapoints in open and closed states
ratio_open_mean = selection_df['Ratio'].loc[mask_open].mean()
ratio_closed_mean = selection_df['Ratio'].loc[mask_closed].mean()
ratio_open_std = selection_df['Ratio'].loc[mask_open].std()
ratio_closed_std = selection_df['Ratio'].loc[mask_closed].std()
sigma_open = np.full((len(selection_df['Ratio'].loc[mask_open])), ratio_open_std)
sigma_closed = np.full((len(selection_df['Ratio'].loc[mask_closed])), ratio_closed_std)

# Perform data fitting using a polynomial function:
# Degree of polynom fit
degree = 1
fit_coef_mask_open, Cov_open = np.polyfit(selection_df['ts'].loc[mask_open], selection_df['Ratio'].loc[mask_open], 
                                          w=1/sigma_open, deg=degree, cov=True) 
fit_coef_mask_closed, Cov_closed = np.polyfit(selection_df['ts'].loc[mask_closed], selection_df['Ratio'].loc[mask_closed], 
                                     w=1/sigma_closed, deg=degree, cov=True)

# Fit 2 degree polynomial
fit_coef_mask_open_2 = np.polyfit(selection_df['ts'].loc[mask_open], selection_df['Ratio'].loc[mask_open], 
                                          w=1/sigma_open, deg=degree+1) 
fit_coef_mask_closed_2 = np.polyfit(selection_df['ts'].loc[mask_closed], selection_df['Ratio'].loc[mask_closed], 
                                     w=1/sigma_closed, deg=degree+1)

# Calculate sigma from covariance matrix:
t = selection_df['ts'].loc[mask] 

# Matrix with rows 1, t, t**2, ...:
TT = np.vstack([t**(degree-i) for i in range(degree+1)]).T

# Matrix multiplication calculates the polynomial values
y_open = np.dot(TT, fit_coef_mask_open)  
y_closed = np.dot(TT, fit_coef_mask_closed) 
C_y_open = np.dot(TT, np.dot(Cov_open, TT.T))
C_y_closed = np.dot(TT, np.dot(Cov_closed, TT.T))

# Standard deviations are sqrt of diagonal
sigma_y_open = np.sqrt(abs(np.diag(C_y_open)))  
sigma_y_closed = np.sqrt(abs(np.diag(C_y_closed)))
rel_uncertainty_closed = sigma_y_closed/y_closed
rel_uncertainty_open = sigma_y_open/y_open   
rel_uncertainty_ratio =  np.sqrt(rel_uncertainty_closed**2 + rel_uncertainty_open**2)    
uncertainty_ratio_fit_all = np.poly1d(fit_coef_mask_closed)(selection_df['ts'].loc[mask]) \
         /np.poly1d(fit_coef_mask_open)(selection_df['ts'].loc[mask])*rel_uncertainty_ratio

# Check data:
fig,(ax1) = plt.subplots(1, figsize=(13,7))
ax1.scatter(selection_df['Time'].loc[mask_open], selection_df['Ratio'].loc[mask_open], c='b', s=30, label='Open state')  
ax1.scatter(selection_df['Time'].loc[mask_closed], selection_df['Ratio'].loc[mask_closed], c='g', s=30, label='Closed state')  
ax1.set_xlabel('Time (CET)', fontsize=18)
ax1.set_ylabel('Ratio', fontsize=18)
plt.title(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4]+ ': Detector ratio for selected data', fontsize=18)
plt.legend(loc="best", fancybox=True, prop={'size':15})
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.gcf().autofmt_xdate()
myFmt = md.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
fig.autofmt_xdate(bottom=0.2, rotation=0, ha='center')
#plt.savefig(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4] + '_dots.jpeg', bbox_inches="tight")
plt.show()

# Plot one cycle
cycle_number= 1
fig,(ax1) = plt.subplots(1, figsize=(13,7))
ax1.scatter(selection_df['Time'].loc[(selection_df['cycle']==cycle_number) & (selection_df['Win_position'] == 'Both Open')], 
            selection_df['Ratio'].loc[(selection_df['cycle']==cycle_number) & (selection_df['Win_position'] == 'Both Open')], c='b', s=30, label='Open state')  
ax1.scatter(selection_df['Time'].loc[(selection_df['cycle']==cycle_number) & (selection_df['Win_position'] == 'Right Cavity Closed')], 
            selection_df['Ratio'].loc[(selection_df['cycle']==cycle_number) & (selection_df['Win_position'] == 'Right Cavity Closed')], c='g', s=30, label='Closed state')  
ax1.set_xlabel('Time (CET)', fontsize=18)
ax1.set_ylabel('Ratio', fontsize=18)
plt.title(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4]+ ': Detector ratio for cycle '+ str(cycle_number), fontsize=18)
plt.legend(loc="best", fancybox=True, prop={'size':15})
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.gcf().autofmt_xdate()
myFmt = md.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
fig.autofmt_xdate(bottom=0.2, rotation=0, ha='center')
#plt.savefig(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4] + '_cycle8.jpeg', bbox_inches="tight")
plt.show()

# Plot mean of cycles:
num_cycles = selection_df['cycle'].max()

fig,(ax1) = plt.subplots(1, figsize=(13,7))
for i in range(1,num_cycles):
    ax1.errorbar(selection_df[(selection_df['cycle']==i) & (selection_df['Win_position'] == 'Both Open')]['Time'].iloc[0], 
            selection_df['Ratio'].loc[(selection_df['cycle']==i) & (selection_df['Win_position'] == 'Both Open')].mean(), 
            fmt='o', c='b', ms=10, yerr= selection_df['Ratio'].loc[(selection_df['cycle']==i) & (selection_df['Win_position'] == 'Both Open')].std() )#, label='Open state')  
    ax1.errorbar(selection_df[(selection_df['cycle']==i) & (selection_df['Win_position'] == 'Right Cavity Closed')]['Time'].iloc[0], 
            selection_df['Ratio'].loc[(selection_df['cycle']==i) & (selection_df['Win_position'] == 'Right Cavity Closed')].mean(),
            fmt='o', c='g', ms=10, yerr= selection_df['Ratio'].loc[(selection_df['cycle']==i) & (selection_df['Win_position'] == 'Right Cavity Closed')].std()) #, label='Closed state')
ax1.set_xlabel('Time (CET)', fontsize=18)
ax1.set_ylabel('Ratio', fontsize=18)
plt.title(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4]+ ': Detector ratio', fontsize=18)
plt.legend(loc="best", fancybox=True, prop={'size':15})
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.gcf().autofmt_xdate()
myFmt = md.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
fig.autofmt_xdate(bottom=0.2, rotation=0, ha='center')
#plt.savefig(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4] + '_mean_values_per_cycle_fits.jpeg', bbox_inches="tight")
plt.show()

# Calculate MITRA transmittance per cycle
# remove first 20 rows and last row within each (cycle, state) group
selection_df_filtered = (selection_df.groupby(["cycle", 'Win_position'], group_keys=False).apply(lambda g: g.iloc[20:-1] if len(g) > 21 else g.iloc[0:0])
      .reset_index(drop=True))
result = (selection_df_filtered.groupby(["cycle", 'Win_position']).agg(mean=('Ratio', "mean"), std=('Ratio', "std"),
          t_start=('Time', "min"), t_end=('Time', "max") ).reset_index())

result["std/mean"] = result["std"] / result["mean"]
result["t_mid"] = result["t_start"] + (result["t_end"] - result["t_start"]) / 2

# pivot to get open/closed stats on one row per cycle
wide = result.pivot(index="cycle", columns = 'Win_position', values=["mean", "std"])

# build a flat dataframe for merging
per_cycle = pd.DataFrame({"cycle": wide.index, "closed/open": wide[("mean", 'Right Cavity Closed')] / wide[("mean", 'Both Open')],
    "err_closed/open": np.sqrt((wide[("std", 'Right Cavity Closed')] / wide[("mean", 'Right Cavity Closed')])**2 + (wide[("std", 'Both Open')] 
                                                                            / wide[("mean", 'Both Open')])**2)}).reset_index(drop=True)

# merge back to both rows of the same cycle
result = result.merge(per_cycle, on="cycle", how="left")
result = result[["cycle", 'Win_position', "t_mid", "mean", "std", "std/mean", "closed/open", "err_closed/open"]]
result['unc_closed/open']= result['err_closed/open']*result['closed/open']

plot_df = (result.groupby("cycle").agg(ratio=("closed/open", "first"), ratio_err=("unc_closed/open", "first"), t_mean=("t_mid", "mean")).reset_index())
plot_df.to_csv('Mecano Data/Ratio_per_cycle_'+filename[0:4]+ '-'+ filename[4:6]+ '-'+ filename[6:8]+ '.csv', sep=';', index=False, header=True)
mean_value = np.nanmean(plot_df['ratio'])
std_value = np.nanstd(plot_df['ratio'])
sem = plot_df['ratio'].sem()

fig,(ax1) = plt.subplots(1, figsize=(13,7))
ax1.errorbar(pd.to_datetime(plot_df['t_mean'], unit='s'), plot_df['ratio'], fmt='o', c='dodgerblue', ms=10, 
             yerr = plot_df['ratio_err'], label='Transmittance per cycle')
ax1.plot([pd.to_datetime(selection_df['ts'].iloc[0], unit='s'), pd.to_datetime(selection_df['ts'].iloc[-1], unit='s')], 
                  [mean_value, mean_value], '--', c='dodgerblue', label='Mean value')
ax1.fill_between((pd.to_datetime(selection_df['ts'].iloc[0], unit='s'), pd.to_datetime(selection_df['ts'].iloc[-1], unit='s')), 
                 (mean_value - 2*std_value, mean_value - 2*std_value), (mean_value + 2*std_value, mean_value + 2*std_value), 
                 color='dodgerblue', alpha=0.2, label=r'2 $\sigma$')
ax1.set_xlabel('Time (CET)', fontsize=18)
ax1.set_ylabel('Transmittance', fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax1.set_ylim([0.92, 0.94])
plt.title(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4]+': Mean +/- sem ' + f'{mean_value:.6f}' + ' +/- ' + f'{sem:.6f}', fontsize=18)
plt.legend(loc="best", fancybox=True, prop={'size':15})
plt.gcf().autofmt_xdate()
myFmt = md.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
fig.autofmt_xdate(bottom=0.2, rotation=0, ha='center')
#plt.savefig(filename[6:8]+'.'+ filename[4:6]+'.'+ filename[0:4] + '_per_cycle.jpeg', bbox_inches="tight")
plt.show()