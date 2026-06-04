### environment.py: Investigate the impact of environmental parameters on the transmittance uncertainty
### Author: Natalia Engler 
### Copyright: 05.2026, MeteoSwiss, PMOD/WRC Davos

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as md
import pvlib

#################################################################################
# Calculate and plot sun azimut and elevation, ZA
#################################################################################
# Settings for Davos
latitude = 46.8021
longitude = 9.8372
altitude = 1560          # meters
timezone = "Europe/Zurich"

day = "2021-10-10"       # example day
# Create timestamps for the whole day (every minute)
times = pd.date_range(start=f"{day} 07:00:00", end=f"{day} 15:00:00", freq="10min", tz=timezone)

# Calculate solar position
solpos = pvlib.solarposition.get_solarposition(time=times, latitude=latitude, longitude=longitude, altitude=altitude)  #this return values in local time
result = pd.DataFrame({"timestamp": times, "azimuth_deg": solpos["azimuth"], "elevation_deg": solpos["elevation"], "zenith_deg": solpos["zenith"]})

plt.figure(figsize=(10, 5))
plt.plot(result["timestamp"], result["azimuth_deg"], label="Azimuth [deg]")
#plt.plot(result["timestamp"], result["elevation_deg"], label="Elevation [deg]")
plt.plot(result["timestamp"], result["zenith_deg"], label="Zenith [deg]")
plt.xlabel("Time")
plt.ylabel("Angle [deg]")
plt.title(f"Sun azimuth and elevation in Davos on {day}")
plt.legend()
plt.grid(True)
myFmt = md.DateFormatter('%H:%M', tz='Europe/Zurich')
#myFmt = md.DateFormatter('%H:%M')   # UTC
plt.gca().xaxis.set_major_formatter(myFmt)
plt.tight_layout()


#################################################################################
#  Open data files, time in UTC
#################################################################################
path = 'Data/'
day = '2025-10-08'

df_data = pd.read_csv(path + 'Ratio_per_cycle_' + day + '.csv', sep=';', engine='python') 
df_DAV = pd.read_csv(path+ 'DAV_N01_'+ day[0:4] + day[5:7] + day[8:10] + '.txt', sep=r'\s+', engine="python", skiprows=22, usecols= [0,1,4,6,7,8], 
                     names= ['hours', 'airmass', 'aod_500', 'alpha','beta', 'flag'])

#  Convert "hours" to timestamp
start_time = pd.Timestamp(day + ' 00:00:00')
df_DAV["time"] = start_time + pd.to_timedelta(df_DAV['hours'], unit="h") + pd.to_timedelta(2*3600, unit='s')

# Add ZA column, sort by time
df_DAV["ZA"]= np.degrees(np.arccos(1 /df_DAV['airmass']))
df_DAV = df_DAV.sort_values("time").reset_index(drop=True)

t_start = pd.Timestamp(day + ' 09:00:00')
t_end   = pd.Timestamp(day + ' 17:00:00')
df_DAV_sel = df_DAV[(df_DAV["time"] >= t_start) & (df_DAV["time"] <= t_end) & (df_DAV["flag"] != 13)].copy()   

# Open table with IWV data:
df_IWV = pd.read_csv(path + 'IWV_IPC-XIV.csv', sep=r",\s+", engine='python', names=['time', 'iwv'], header=None)
df_IWV['time'] = pd.to_datetime(df_IWV['time'], format='%d.%m.%y %H:%M', errors='coerce') + pd.to_timedelta(2*3600, unit='s')
df_IWV_sel = df_IWV[(df_IWV['time'] >= t_start) & (df_IWV['time'] <= t_end)].copy()

#  Calculate rate of change for parameters
columns_to_plot = ['ZA', 'aod_500', "alpha", 'beta']
dt_seconds = df_DAV_sel["time"].diff().dt.total_seconds()
myFmt = md.DateFormatter('%H:%M')
for col in columns_to_plot:
    df_DAV_sel[f"{col}_rate"] = df_DAV_sel[col].diff() / dt_seconds


# Parse t_mean in df_data
def parse_t_mean(x, day):
    """
    Handles t_mean stored as:
    - full datetime string
    - HH:MM:SS
    - HH:MM
    - numeric hours from midnight
    """
    if pd.isna(x):
        return pd.NaT
    # numeric: assume hours from midnight
    if isinstance(x, (int, float, np.integer, np.floating)):
        return pd.Timestamp(day + " 00:00:00") + pd.to_timedelta(x, unit="h")
    x_str = str(x).strip()
    # try full datetime first
    t = pd.to_datetime(x_str, errors="coerce")
    if pd.notna(t):
        # If pandas interpreted it without the correct date, fix it below
        if t.date() != pd.Timestamp(day).date() and len(x_str) <= 8:
            return pd.Timestamp(day + " " + x_str)
        return t
    # fallback: assume time of day
    return pd.Timestamp(day + " " + x_str)

df_data = df_data.copy()
df_data["t_mean_dt"] = df_data["t_mean"].apply(lambda x: parse_t_mean(x, day))
df_data = df_data.sort_values("t_mean_dt").reset_index(drop=True)

# Build intervals centered on each t_mean. Boundaries are halfway between neighboring t_mean values
t = df_data["t_mean_dt"]
left_edges = []
right_edges = []

for i in range(len(df_data)):
    if i == 0:
        # first interval: use half-distance to next t_mean
        half_width = (t.iloc[i + 1] - t.iloc[i]) / 2
        left = t.iloc[i] - half_width
    else:
        left = t.iloc[i] - (t.iloc[i] - t.iloc[i - 1]) / 2

    if i == len(df_data) - 1:
        # last interval: use half-distance to previous t_mean
        half_width = (t.iloc[i] - t.iloc[i - 1]) / 2
        right = t.iloc[i] + half_width
    else:
        right = t.iloc[i] + (t.iloc[i + 1] - t.iloc[i]) / 2
    left_edges.append(left)
    right_edges.append(right)

df_data["t_left"] = left_edges
df_data["t_right"] = right_edges

# Restrict intervals to selected time window
df_data["t_left"] = df_data["t_left"].clip(lower=t_start)
df_data["t_right"] = df_data["t_right"].clip(upper=t_end)

# Calculate IWV rate
df_IWV_sel = df_IWV_sel.sort_values("time").reset_index(drop=True)

# Rate between current IWV value and next hourly IWV value
df_IWV_sel["iwv_next"] = df_IWV_sel["iwv"].shift(-1)
df_IWV_sel["time_next"] = df_IWV_sel["time"].shift(-1)
df_IWV_sel["dt_next_seconds"] = (df_IWV_sel["time_next"] - df_IWV_sel["time"]).dt.total_seconds()
df_IWV_sel["iwv_rate"] = (df_IWV_sel["iwv_next"] - df_IWV_sel["iwv"]) / df_IWV_sel["dt_next_seconds"]

# Each IWV rate is assumed constant from time to time_next
df_IWV_rate_segments = df_IWV_sel.dropna(subset=["iwv_rate", "time_next"]).copy()

# Mean rates inside each t_mean-centered interval
mean_ZA_values = []
mean_aod_500_values = []
mean_alpha_values = []
mean_beta_values = []
mean_iwv_values = []
mean_ZA_rates = []
mean_aod_500_rates = []
mean_alpha_rates = []
mean_beta_rates = []
mean_iwv_rates = []
n_ZA = []
n_aod_500 = []
n_alpha = []
n_beta = []
n_iwv_segments = []

for _, row in df_data.iterrows():
    mask_dav = ((df_DAV_sel["time"] >= row["t_left"]) & (df_DAV_sel["time"] < row["t_right"]))
    dav_interval = df_DAV_sel.loc[mask_dav]
    mean_ZA_values.append(dav_interval["ZA"].mean())
    mean_aod_500_values.append(dav_interval["aod_500"].mean())
    mean_alpha_values.append(dav_interval["alpha"].mean())
    mean_beta_values.append(dav_interval["beta"].mean())
    mean_ZA_rates.append(dav_interval["ZA_rate"].abs().mean())
    mean_aod_500_rates.append(dav_interval["aod_500_rate"].abs().mean())
    mean_alpha_rates.append(dav_interval["alpha_rate"].abs().mean())
    mean_beta_rates.append(dav_interval["beta_rate"].abs().mean())
    n_ZA.append(dav_interval["ZA_rate"].notna().sum())
    n_aod_500.append(dav_interval["aod_500_rate"].notna().sum())
    n_alpha.append(dav_interval["alpha_rate"].notna().sum())
    n_beta.append(dav_interval["beta_rate"].notna().sum())

    # IWV rate averaged over overlap
    interval_start = row["t_left"]
    interval_end = row["t_right"]
    weighted_rate_sum = 0.0
    weighted_iwv_sum = 0.0
    total_overlap_seconds = 0.0
    used_segments = 0
    for _, seg in df_IWV_rate_segments.iterrows():
        seg_start = seg["time"]
        seg_end = seg["time_next"]
        overlap_start = max(interval_start, seg_start)
        overlap_end = min(interval_end, seg_end)
        overlap_seconds = (overlap_end - overlap_start).total_seconds()

        if overlap_seconds > 0:
            # Rate is assumed constant over the segment
            weighted_rate_sum += seg["iwv_rate"] * overlap_seconds
            # IWV value is assumed constant from seg_start to seg_end
            weighted_iwv_sum += seg["iwv"] * overlap_seconds
            total_overlap_seconds += overlap_seconds
            used_segments += 1
    if total_overlap_seconds > 0:
        mean_iwv_rate = weighted_rate_sum / total_overlap_seconds
        mean_iwv_value = weighted_iwv_sum / total_overlap_seconds
    else:
        mean_iwv_rate = np.nan
        mean_iwv_value = np.nan
    mean_iwv_rates.append(mean_iwv_rate)
    mean_iwv_values.append(mean_iwv_value)
    n_iwv_segments.append(used_segments)

df_data["mean_ZA"] = mean_ZA_values
df_data["mean_aod_500"] = mean_aod_500_values
df_data["mean_alpha"] = mean_alpha_values
df_data["mean_beta"] = mean_beta_values
df_data["mean_iwv"] = mean_iwv_values
df_data["mean_ZA_rate"] = mean_ZA_rates
df_data["mean_aod_500_rate"] = mean_aod_500_rates
df_data["mean_alpha_rate"] = mean_alpha_rates
df_data["mean_beta_rate"] = mean_beta_rates
df_data["mean_iwv_rate"] = mean_iwv_rates
df_data["n_ZA_rate"] = n_ZA
df_data["n_aod_500_rate"] = n_aod_500
df_data["n_alpha_rate"] = n_alpha
df_data["n_beta_rate"] = n_beta
df_data["n_iwv_segments"] = n_iwv_segments

df_data.to_csv('Data/Table_' + day + '.csv', sep=';', index=False, header=True)

#################################################################################
# Open saved tables for 6 days. All tables have local time
#################################################################################
days = ["2025-09-29", "2025-10-02", "2025-10-08", "2026-04-07", "2026-04-08", "2026-04-09"]
tables = []
for day in days:
    filename = path + f"Table_{day}.csv"
    df_tmp = pd.read_csv(filename, sep=';')
    # Add day label for plotting
    df_tmp["day"] = day
    tables.append(df_tmp)

df_all = pd.concat(tables, ignore_index=True)

rate_cols = ["mean_ZA_rate", "mean_aod_500_rate", "mean_alpha_rate", "mean_beta_rate", "mean_iwv_rate" ]  #, "mean_iwv_rate"
mean_cols = ["mean_ZA", "mean_aod_500", "mean_alpha", "mean_beta", "mean_iwv"]   #, "mean_iwv"

# Plot ratio_err versus absolute mean rates for both days
rate_columns = [ ("mean_ZA_rate", r"$|\overline{dZA/dt}|$  (deg/s)", "Zenith angle"), ("mean_aod_500_rate", r"$|\overline{dAOD_{500}/dt}|$  (1/s)", r"$AOD_{500}$"),
    ("mean_alpha_rate", r"$|\overline{d\alpha/dt|}$  (1/s)", "Alpha"), ("mean_beta_rate", r"$|\overline{d\beta/dt}|$   (1/s)", "Beta"),
    ("mean_iwv_rate", r"$|\overline{dIWV/dt}|$  (mm/s)", "IWV")]   # 

fig, axes = plt.subplots(5, 1, figsize=(9, 14.5), sharex=False) #15
for ax, (rate_col, xlabel, title) in zip(axes, rate_columns):
    for day in days:
        df_day = df_all[df_all["day"] == day].copy()
        valid = df_day[["ratio", "ratio_err", rate_col, "cycle"]].dropna().copy()
        valid[f"abs_{rate_col}"] = valid[rate_col].abs()
        ax.scatter(valid[f"abs_{rate_col}"], valid["ratio_err"]/valid["ratio"]*1e6, label=day, alpha=0.8)

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Rel. uncertainty (ppm)", fontsize=14)
    #ax.set_title(f"ratio_err versus absolute {title} rate")
    ax.grid(True, alpha=0.3)
axes[1].set_xlim(-0.00001, 0.0005)
plt.legend(bbox_to_anchor=(0.21, 0), loc="lower left",  bbox_transform=fig.transFigure, ncol=3, fontsize=13)
plt.tight_layout()
#plt.savefig('Parameter_rates_6days.jpeg', bbox_inches="tight")
plt.show()

# Plot ratio_err versus means for both days
mean_columns = [("mean_ZA", r"$\overline{ZA}$  (deg)", "Zenith angle"), ("mean_aod_500", r"$\overline{AOD_{500}}$", r"$AOD_{500}$"),
    ("mean_alpha", r"$\overline{\alpha}$", "Alpha"), ("mean_beta", r"$\overline{\beta}$", "Beta"), ("mean_iwv", r"$\overline{IWV}$  (mm)", "IWV")] # ]

fig, axes = plt.subplots(5, 1, figsize=(9, 14.5), sharex=False)
for ax, (mean_col, xlabel, title) in zip(axes, mean_columns):
    for day in days:
        df_day = df_all[df_all["day"] == day].copy()
        valid = df_day[['ratio',"ratio_err", mean_col, "cycle"]].dropna().copy()
        ax.scatter(valid[f"{mean_col}"], valid["ratio_err"]/valid["ratio"]*1e6, label=day, alpha=0.8)

    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Rel. uncertainty (ppm)", fontsize=14)
    #ax.set_title(f"ratio_err versus {title}")
    ax.grid(True, alpha=0.3)
    #ax.legend()
plt.legend(bbox_to_anchor=(0.21, 0), loc="lower left",  bbox_transform=fig.transFigure, ncol=3, fontsize=13)
plt.tight_layout()
#plt.savefig('Parameter_6days.jpeg', bbox_inches="tight")
plt.show()

# Correlations using all days together
print("Correlations using all days together:\n")
for rate_col, _, title in rate_columns:
    valid = df_all[["ratio_err", rate_col]].dropna().copy()
    valid[f"abs_{rate_col}"] = valid[rate_col].abs()
    if len(valid) > 1:
        pearson = valid["ratio_err"].corr(valid[f"abs_{rate_col}"], method="pearson")
        spearman = valid["ratio_err"].corr(
            valid[f"abs_{rate_col}"],
            method="spearman")
        print(f"{title}:")
        print(f"  Pearson  = {pearson:.3f}")
        print(f"  Spearman = {spearman:.3f}")
        print(f"  N points = {len(valid)}\n")
    else:
        print(f"{title}: not enough valid data\n")

for mean_col, _, title in mean_columns:
    valid = df_all[["ratio_err", mean_col]].dropna().copy()
    valid[f"abs_{mean_col}"] = valid[mean_col].abs()
    if len(valid) > 1:
        pearson = valid["ratio_err"].corr(valid[f"abs_{mean_col}"], method="pearson")
        spearman = valid["ratio_err"].corr(
            valid[f"abs_{mean_col}"],
            method="spearman")
        print(f"{title}:")
        print(f"  Pearson  = {pearson:.3f}")
        print(f"  Spearman = {spearman:.3f}")
        print(f"  N points = {len(valid)}\n")
    else:
        print(f"{title}: not enough valid data\n")


#################################################################################
#  Correlation with the wind speed and direction
#################################################################################

# Die Daten mit der Aufzeichnung der Wind Richtung und speed sind in csv-Format unter:
# "//ad.pmodwrc.ch/Institute/Departments/WRC/SRS/Mecano/PC208W/data/"
basefile = './Mecano Data/'
filename = 'py260408.csv'  

W_df = pd.read_table(basefile + filename, header=None, names = ['Jahr', 'Tag', 'Zeit', 'Wspeed', 'Wdir'], skiprows=0, usecols = [1,2,3,23,25], delimiter=',')   
# Make sure the columns are strings with leading zeros where needed
W_df["Jahr"] = W_df["Jahr"].astype(int).astype(str)
W_df["Tag"] = W_df["Tag"].astype(int).astype(str).str.zfill(3)   # day of year
W_df["Zeit"] = W_df["Zeit"].astype(int).astype(str).str.zfill(4) # HHMM

# Remove rows with invalid Zeit = 2400
W_df = W_df[W_df["Zeit"] != '2400'].copy()

# Create timestamp from year + day-of-year + HHMM
W_df["Time"] = pd.to_datetime(W_df["Jahr"] + W_df["Tag"] + W_df["Zeit"], format="%Y%j%H%M")

# Add 2*3600 seconds = 2 hours (UTC -> CEST Sommerzeit)
W_df["Time"] = W_df["Time"] + pd.to_timedelta(2*3600, unit="s")
#print(W_df)

# Visualize wind direction
# Convert degrees to radians
theta = np.deg2rad(W_df["Wdir"])
# Meteorological convention: direction wind comes FROM, clockwise from north.
W_df["u"] = -W_df["Wspeed"] * np.sin(theta)  # u is east-west component
W_df["v"] = -W_df["Wspeed"] * np.cos(theta)  # v is north-south component

# define time interval
start_t = pd.to_datetime("08:00:00").time()
end_t   = pd.to_datetime("16:00:00").time()
W_df_interval = W_df[ (W_df["Time"].dt.time >= start_t) & (W_df["Time"].dt.time <= end_t)].copy()
vmax = np.nanmax(np.abs(W_df_interval["v"]))

fig, ax = plt.subplots(figsize=(12, 6))
# Put all arrows on y=0 line
y = np.zeros(len(W_df_interval))
ax.quiver(W_df_interval["Time"], y, W_df_interval["u"], W_df_interval["v"], angles="uv", scale_units="y", scale=1)
#ax.plot(W_df_interval["Time"],W_df_interval["Wdir"], marker=".", ls="None")
ax.set_ylim(-1.1 * vmax, 1.1 * vmax)
ax.set_xlabel("Time", fontsize=18)
ax.set_title("08.04.2026: Wind magnitude and direction mesured by WSG sensor", fontsize=18)
ax.set_ylabel("N-S component of wind speed (m/s)", fontsize=18)
ax.xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax.grid(True)
plt.tight_layout()
#plt.savefig('2026-04-08_wind_speed_WSG.jpeg', bbox_inches="tight")
plt.show()

df_MS = pd.read_csv(basefile +"ogd-smn_dav_t_historical_2020-2029.csv", sep=r";", engine="python", names = ['Time', 'Wspeed', 'Wdir'], skiprows=0, usecols = [1,16,17])

# Convert columns to proper types
df_MS["Time"] = pd.to_datetime(df_MS["Time"], format="%d.%m.%Y %H:%M", errors="coerce") + pd.to_timedelta(1*3600, unit="s")
df_MS["Wspeed"] = pd.to_numeric(df_MS["Wspeed"], errors="coerce")
df_MS["Wdir"] = pd.to_numeric(df_MS["Wdir"], errors="coerce")

# Day and time must correspond day and time in first code to calculate sol_positions:
day = "2021-10-10"

# Filter interval
df_MS_interval = df_MS[(df_MS["Time"] >= pd.Timestamp(day + ' 08:00:00')) & (df_MS["Time"] <= pd.Timestamp(day + ' 16:00:00'))].copy()
print(df_MS_interval["Time"])
# Wind direction in radians
theta_MS_interval = np.deg2rad(df_MS_interval["Wdir"])
# Meteorological convention: direction wind comes FROM, clockwise from north.
df_MS_interval["u"] = -df_MS_interval["Wspeed"] * np.sin(theta_MS_interval)  # u is east-west component
df_MS_interval["v"] = -df_MS_interval["Wspeed"] * np.cos(theta_MS_interval)  # v is north-south component
print(df_MS_interval["v"])
ymax_abs = np.nanmax(np.abs(df_MS_interval["v"]))

fig, ax = plt.subplots(figsize=(12, 6))
# Put all arrows on y=0 line
y0 = np.zeros(len(df_MS_interval))
ax.quiver(df_MS_interval["Time"], y0, df_MS_interval["u"], df_MS_interval["v"], angles="uv", scale_units="y", scale=1)
#ax.quiver(W_df_interval["Time"], y, W_df_interval["u"], W_df_interval["v"], angles="uv", color='r', scale_units="xy", scale=1)
#ax.plot(df_MS_interval["Time"], df_MS_interval["Wdir"], marker=".", ls="None")
ax.set_ylim(-1.5 * ymax_abs, 1)
ax.set_ylabel("N-S component of wind speed (m/s)", fontsize=18)
ax.set_xlabel("Time", fontsize=18)
ax.set_title("14.10.2021: Wind magnitude and direction for Davos (averaged over 10 min)", fontsize=18)
ax.xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
ax.grid(True)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
#plt.savefig(f'{day}_wind_speed_Davos.jpeg', bbox_inches="tight")
plt.show()

# angle difference in degrees
azimuth = result["azimuth_deg"].to_numpy()
wdir = df_MS_interval["Wdir"].to_numpy()
wspeed = df_MS_interval["Wspeed"].to_numpy()

delta_angle_deg = (azimuth - wdir)
delta_angle_rad = np.deg2rad(delta_angle_deg)
df_MS_interval["v"] = -wspeed * np.cos(delta_angle_rad)

# projection of wind speed onto perpendicular sun azimuth direction
df_MS_interval["u"] = wspeed * np.sin(delta_angle_rad)

fig, ax = plt.subplots(figsize=(12, 6))
# Put all arrows on y=0 line
y0 = np.zeros(len(df_MS_interval))
ax.quiver(df_MS_interval["Time"], y0, df_MS_interval["u"], df_MS_interval["v"], angles="uv", scale_units="y", scale=1)
#ax.quiver(W_df_interval["Time"], y, W_df_interval["u"], W_df_interval["v"], angles="uv", color='r', scale_units="xy", scale=1)
#ax.plot(df_MS_interval["Time"], df_MS_interval["Wdir"], marker=".", ls="None")
ax.set_ylim(-2, 4)
ax.set_ylabel("Sun direction component of wind speed (m/s)", fontsize=18)
ax.set_xlabel("Time", fontsize=18)
ax.set_title("14.10.2021:: Projection into sun azimuth", fontsize=18)
ax.xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
ax.grid(True)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
#plt.savefig(f'{day}_wind_speed_Sun_projection.jpeg', bbox_inches="tight")
plt.show()