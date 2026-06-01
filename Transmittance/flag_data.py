### flag_data.py: interactive data flagging tool to inspect the data in a time series plot and assign flags to individual data points or groups of data points.
### Author: Natalia Engler Copyright: 2026, MeteoSwiss, PMOD/WRC Davos

import tkinter as tk
from tkinter import filedialog, messagebox
from pathlib import Path
from typing import Any
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.widgets import Button, RadioButtons, CheckButtons
from typing import Optional
import numpy as np

# Configuration
plt.rcParams.update({"figure.figsize": (16, 9), "figure.dpi": 120, "font.size": 20, "axes.titlesize": 22, "axes.labelsize": 20, "xtick.labelsize": 18,
    "ytick.labelsize": 18, "legend.fontsize": 18,})
dtm = 'dtm'
flags = '_flag_'
colors = '_color_'
flag_col_prefix = 'flag_'

# State / key mapping
keys: dict[str, dict[str, Any]] = {
    'escape': {'flag': None, 'color': 'magenta', 'meaning': 'unflagged'},
    '0': {'flag': 0, 'color': 'blue', 'meaning': 'valid'},
    '1': {'flag': 1, 'color': 'red', 'meaning': 'invalid'}}

# Columns that should not appear as selectable variables
exclude_variables = {dtm, 'sec_str', 'V_cav', 'V_ref', 'V_dark', 'Current', 'Current_flag'}

# Functions
def _order_key_items(keys: dict[str, dict[str, Any]]) -> list[tuple[str, dict[str, Any]]]:
    items: list[tuple[str, dict[str, Any]]] = []
    if 'escape' in keys:
        items.append(('escape', keys['escape']))
    items += [(k, keys[k]) for _, k in sorted((int(k), k) for k in keys if k.isdigit())]
    items += [(k, v) for k, v in keys.items() if k not in {'escape'} and not k.isdigit()]
    return items


def add_legend_below_axes(ax, keys: dict[str, dict[str, Any]], ncol: Optional[int] = None, bottom_pad: float = 0.22):

    fig = ax.figure

    # remove previous legends
    for lg in list(fig.legends):
        try:
            lg.remove()
        except Exception:
            pass

    handles = []
    for _, info in _order_key_items(keys):
        flag = info.get('flag', '')
        color = info.get('color', '#888')
        meaning = info.get('meaning', '')
        handles.append(Line2D([0], [0], marker='o', linestyle='', mfc=color, markeredgecolor='#333', ms=16, label=f'{flag}: {meaning}'))

    if ncol is None:
        ncol = min(6, len(handles)) if handles else 1
    fig.subplots_adjust(bottom=bottom_pad, left=0.22)
    fig.legend(handles=handles, loc='lower center', ncol=ncol, bbox_to_anchor=(0.5, 0.01), frameon=True)


def try_read_txt(path: Path) -> pd.DataFrame:
    '''
    Read txt/csv files. Add column 'Ratio' with calculated 
    '''
    try:
        MITRA_df = pd.read_table(path, header=None, names = ['Time', 'sec_str', 'V_cav', 'V_ref', 'V_dark', 'Current', 'Current_flag', 
                          'Win_position'], skiprows=1, sep=None, engine='python')
        length_data = len(MITRA_df['sec_str'])
        R_cav = np.zeros(length_data-2)
        R_ref = np.zeros(length_data-2)
        R_dark = np.zeros(length_data-2)
        for i in range(1, length_data-1):
            R_cav[i-1] = 0.5*(0.5*(abs(MITRA_df['V_cav'].iloc[i-1]/MITRA_df['Current'].iloc[i-1]) + abs(MITRA_df['V_cav'].iloc[i]/MITRA_df['Current'].iloc[i]) ) +0.5*(abs(MITRA_df['V_cav'].iloc[i+1]/MITRA_df['Current'].iloc[i+1]) + abs(MITRA_df['V_cav'].iloc[i]/MITRA_df['Current'].iloc[i]) ))
            R_ref[i-1] = 0.5*(0.5*(abs(MITRA_df['V_ref'][i-1]/MITRA_df['Current'][i-1]) + abs(MITRA_df['V_ref'][i]/MITRA_df['Current'][i]) ) + 0.5*(abs(MITRA_df['V_ref'][i+1]/MITRA_df['Current'][i+1]) + abs(MITRA_df['V_ref'][i]/MITRA_df['Current'][i]) ))
            R_dark[i-1] = 0.5*(0.5*(abs(MITRA_df['V_dark'][i-1]/MITRA_df['Current'][i-1]) + abs(MITRA_df['V_dark'][i]/MITRA_df['Current'][i]) ) + 0.5*(abs(MITRA_df['V_dark'][i+1]/MITRA_df['Current'][i+1]) + abs(MITRA_df['V_dark'][i]/MITRA_df['Current'][i]) ))

        # Old TKs and R0 from temperature fits (R0 are at 0°C), used for 20220921 and IPC-13:
        T0 = 0
        TK_ref, R0_ref =  1.57222846, 300.78440976   # Sensor ID 21-003
        TK_cav, R0_cav =  1.50948573, 299.65496045  # Sensor ID 21-022
        TK_dark, R0_dark = 1.49504373, 296.01205954  # Sensor ID 21-025

        """
        # Coefficients 09.03.2026:
        TK_ref, R0_ref  =  331.5758*0.0047281, 331.5758   # Sensor ID 21-003  0.00506617
        TK_cav, R0_cav  =  328.9755*0.0045752, 328.9755   # Sensor ID 21-022   0.0112037728
        TK_dark, R0_dark = 325.0509*0.0045872, 325.0509   # Dark
        """
        
        T_cav = (R_cav/R0_cav - 1)*R0_cav/TK_cav +19
        T_ref = (R_ref/R0_ref - 1)*R0_ref/TK_ref +19
        T_dark = (R_dark/R0_dark - 1)*R0_dark/TK_dark +19
        ratio = (T_cav-T_dark)/(T_ref-T_dark)
        MITRA_df['Ratio'] = np.r_[1, ratio, 1]

        return MITRA_df

    except Exception:
        pass
        raise ValueError(f'Could not read TXT file: {path}')


def find_dtm_column(df: pd.DataFrame) -> str:
    '''
    Prefer exact 'dtm', otherwise try common datetime column names.
    '''
    candidates = ['dtm', 'datetime', 'DateTime', 'timestamp', 'Timestamp', 'time', 'Time', 'date', 'Date']
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError('No datetime column found. Expected for example "dtm", "datetime", or "timestamp".')


def is_numeric_series(s: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(s)


# Main app
class FlagData:
    def __init__(self, width: int = 18, height: int = 9):
        self.df: pd.DataFrame | None = None
        self.selected_file: Path | None = None
        self.variable: str | None = None
        self.dtm_col: str | None = None
        self.scatter = None

        # Tkinter root only for dialogs
        self.root = tk.Tk()
        self.root.withdraw()
        self.fig, self.ax = plt.subplots(figsize=(width, height))
        try:
            self.fig.canvas.manager.set_window_title('Flag data')
        except Exception:
            pass

        # top buttons
        ax_open = plt.axes([0.01, 0.9, 0.15, 0.06])
        ax_save = plt.axes([0.01, 0.8, 0.15, 0.06])
        ax_check = plt.axes([0.01, 0.68, 0.15, 0.1])
        self.btn_open = Button(ax_open, 'Open')
        self.btn_save = Button(ax_save, 'Save')
        self.check_remove_invalid = CheckButtons(ax_check, ["Remove invalid\nrows on save"], [False])
        self.ax_radio = plt.axes([0.01, 0.32, 0.15, 0.3])
        self.radio = RadioButtons(self.ax_radio, ['-'], active=0)
        self.status_text = self.fig.text(0.22, 0.95, 'Open a TXT file to start', fontsize=20)
        self.btn_open.on_clicked(self.on_open_file)
        self.btn_save.on_clicked(self.on_save_file)
        self.radio.on_clicked(self.on_variable_selected)
        self.fig.canvas.mpl_connect('pick_event', self.on_picked_flag_point)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_pressed_flag_points)
        add_legend_below_axes(self.ax, keys, bottom_pad=0.15)

    def set_status(self, text: str):
        self.status_text.set_text(text)
        self.fig.canvas.draw_idle()

    def refresh_variable_selector(self):
        options = self.get_variable_options()
        if not options:
            options = ['-']
        self.ax_radio.clear()
        self.radio = RadioButtons(self.ax_radio, options, active=0)
        self.radio.on_clicked(self.on_variable_selected)
        self.ax_radio.set_title('Variable', fontsize=16)
        self.fig.canvas.draw_idle()
        if options and options[0] != '-':
            self.on_variable_selected(options[0])

    def remove_invalid_on_save_enabled(self) -> bool:
        return self.check_remove_invalid.get_status()[0]

    def get_variable_options(self) -> list[str]:
        if self.df is None:
            return []
        options = []
        for c in self.df.columns:
            if c in exclude_variables:
                continue
            if str(c).startswith(flag_col_prefix):
                continue
            if c == self.dtm_col:
                continue
            if is_numeric_series(self.df[c]):
                options.append(c)
        return options

    def ensure_plotting_columns_for_variable(self, variable: str):
        '''
        Create or refresh:
        - _flag_   : current plotting flag column
        - _color_  : current plotting color column
        '''
        if self.df is None:
            return

        f_variable = f'{flag_col_prefix}{variable}'
        if f_variable not in self.df.columns:
            self.df[f_variable] = pd.Series([pd.NA] * len(self.df), dtype='Int64')

        self.df[flags] = self.df[f_variable]
        self.df[colors] = keys['escape']['color']
        for _, info in _order_key_items(keys):
            flag = info['flag']
            color = info['color']
            if flag is not None:
                self.df.loc[self.df[flags] == flag, colors] = color

    # File actions
    def on_open_file(self, event=None):
        file_path = filedialog.askopenfilename(
            title='Select TXT file',
            filetypes=[('TXT files', '*.txt'), ('All files', '*.*')])
        if not file_path:
            return
        path = Path(file_path)

        try:
            df = try_read_txt(path)
            dtm_col = find_dtm_column(df)
            df[dtm_col] = pd.to_datetime(df[dtm_col], errors='coerce')
            df = df.dropna(subset=[dtm_col]).sort_values(dtm_col).reset_index(drop=True)
            self.df = df
            self.selected_file = path
            self.dtm_col = dtm_col
            self.variable = None
            self.refresh_variable_selector()
            self.set_status(f'Opened: {path}')

        except Exception as e:
            messagebox.showerror('Open failed', str(e))
            self.set_status(f'Open failed: {e}')

    def on_save_file(self, event=None):
        if self.df is None:
            self.set_status('No data loaded.')
            return

        # Ensure selected variable's plot flags are already stored in f_<var>
        out = self.df.copy()

        # Optional: remove rows flagged as invalid for current variable
        if self.remove_invalid_on_save_enabled() and self.variable is not None:
            f_variable = f"{flag_col_prefix}{self.variable}"
            print(f_variable)
            if f_variable in out.columns:
                flag_values = pd.to_numeric(out[f_variable], errors="coerce")
                print(flag_values.value_counts(dropna=False))
                mask = flag_values.ne(1) | flag_values.isna()
                out = out[mask].copy()
                print(out)     

        # Remove helper columns used only for plotting
        out = out.drop(columns=[flags, colors], errors='ignore')
        save_path = filedialog.asksaveasfilename(
            title='Save flagged TXT file',
            defaultextension='.txt',
            filetypes=[('TXT files', '*.txt'), ('All files', '*.*')],
            initialfile=(self.selected_file.stem + '_flagged.txt') if self.selected_file else 'flagged.txt')
        if not save_path:
            return

        save_path = Path(save_path)
        try:
            out.to_csv(save_path, sep='\t', index=False)
            self.set_status(f'Saved: {save_path}')
        except Exception as e:
            messagebox.showerror('Save failed', str(e))
            self.set_status(f'Save failed: {e}')

    # Plots
    def on_variable_selected(self, label):
        if self.df is None or label == '-':
            return
        self.variable = label
        self.plot_current_variable()

    def plot_current_variable(self):
        if self.df is None or self.variable is None or self.dtm_col is None:
            return
        if self.variable not in self.df.columns:
            self.set_status(f'Variable not found: {self.variable}')
            return
        if self.df[self.variable].dropna().empty:
            self.set_status(f'No data available for {self.variable}')
            return
        self.ensure_plotting_columns_for_variable(self.variable)

        self.ax.clear()
        """
        title = 'Data flagging'
        if self.selected_file is not None:
            title += f'\n{self.selected_file}'
        self.ax.set_title(title)
        """
        self.scatter = self.ax.scatter(self.df[self.dtm_col], self.df[self.variable], c=self.df[colors].tolist(), alpha=0.7, s=15, picker=5)
        self.ax.set_xlabel(self.dtm_col)
        self.ax.set_ylabel(self.variable)
        self.ax.grid(True, alpha=0.3)
        self.ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))  #('%Y-%m-%d\n%H:%M')
        add_legend_below_axes(self.ax, keys, bottom_pad=0.15)
        self.fig.canvas.draw_idle()
        self.set_status(f'Selected variable: {self.variable}')

    # Interaction
    def on_picked_flag_point(self, event):
        '''
        Flag single picked point(s), but only if zoom/pan mode is OFF.
        '''
        if self.df is None or self.variable is None or self.scatter is None:
            return

        if self.ax.get_navigate_mode() is not None:
            self.set_status('Disable zoom/pan before clicking single points.')
            return
        
        key = event.mouseevent.key
        if key not in keys:
            self.set_status(f'Point selected, but key "{key}" is not assigned.')
            return

        flag = keys[key]['flag']
        color = keys[key]['color']
        f_variable = f'{flag_col_prefix}{self.variable}'
        idx = event.ind
        self.df.loc[idx, flags] = flag
        self.df.loc[idx, f_variable] = flag
        if flag is None:
            self.df.loc[idx, colors] = keys['escape']['color']
        else:
            self.df.loc[idx, colors] = color

        self.scatter.set_color(self.df[colors].tolist())
        self.fig.canvas.draw_idle()
        self.set_status(f'Point(s) {list(idx)} flagged as {flag}')

    def on_key_pressed_flag_points(self, event):
        '''
        If zoom mode is ON, apply current key-flag to all visible points.
        If zoom mode is OFF, only show status.
        '''
        if self.df is None or self.variable is None or self.scatter is None or self.dtm_col is None:
            return

        if event.key not in keys:
            return

        flag = keys[event.key]['flag']
        color = keys[event.key]['color']
        meaning = keys[event.key]['meaning']
        f_variable = f'{flag_col_prefix}{self.variable}'

        # only bulk-flag when zoom/pan is active
        if self.ax.get_navigate_mode() is None:
            self.set_status(f'Key "{event.key}" = flag {flag} ({meaning})')
            return

        zoom_xlim = self.ax.get_xlim()
        zoom_ylim = self.ax.get_ylim()
        x0 = matplotlib.dates.num2date(zoom_xlim[0]).replace(tzinfo=None)
        x1 = matplotlib.dates.num2date(zoom_xlim[1]).replace(tzinfo=None)
        y0, y1 = zoom_ylim
        condition = ((self.df[self.dtm_col] > x0) & (self.df[self.dtm_col] < x1) & (self.df[self.variable] > y0) & (self.df[self.variable] < y1))
        self.df.loc[condition, flags] = flag
        self.df.loc[condition, f_variable] = flag
        if flag is None:
            self.df.loc[condition, colors] = keys['escape']['color']
        else:
            self.df.loc[condition, colors] = color
        self.scatter.set_color(self.df[colors].tolist())
        self.fig.canvas.draw_idle()
        n = int(condition.sum())
        self.set_status(f'{n} visible point(s) flagged as {flag} ({meaning})')


    def run(self):
        plt.show()

if __name__ == '__main__':
    app = FlagData(width=20, height=10)
    app.run()