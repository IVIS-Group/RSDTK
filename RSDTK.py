"""
Remote Sensing Data Toolkit (RSDTK) v1.2.1
A comprehensive GUI application for batch processing and converting
satellite and airborne remote sensing data to analysis-ready formats.
 
Supported Sensors:
  Spaceborne Multispectral:
    - ASTER (V004 GeoTIFF + V003 HDF products)
    - Landsat 8/9 (Collection 2 L1 and L2SP)
    - Sentinel-2 (L2A Surface Reflectance)
    - ECOSTRESS (L1CT Radiance, L2 LSTE)
    - VIIRS (VNP02IMG, VNP02MOD, VNP21)
 
  Spaceborne Hyperspectral:
    - EMIT (L2A Surface Reflectance)
    - AIRS (L1B Radiance)
 
  Airborne:
    - MASTER (L1B Radiance, L2 Emissivity/LST)
    - HyTES (L1 Radiance, L2 Emissivity/LST)
 
Dependencies: customtkinter, rasterio, numpy, scipy, h5py, pyhdf
"""
 
import os
import sys
import threading
import queue
from datetime import datetime
 
try:
    import customtkinter as ctk
except ImportError:
    print("ERROR: customtkinter is required. Install with: pip install customtkinter")
    sys.exit(1)
 
from tkinter import filedialog, messagebox
 
 
# ============================================================================
# Application configuration
# ============================================================================
 
APP_NAME    = "Remote Sensing Data Toolkit"
APP_VERSION = "1.2.1"
APP_SIZE    = "1200x900"
APP_MIN_SIZE = (1000, 800)
 
COLORS = {
    'sidebar_bg':     '#1a1a2e',
    'sidebar_hover':  '#16213e',
    'sidebar_active': '#0f3460',
    'accent':         '#e94560',
    'accent_hover':   '#c73e54',
    'success':        '#2ecc71',
    'warning':        '#f39c12',
    'error':          '#e74c3c',
    'text_primary':   ('#1a1a2e', '#ffffff'),   # (light_mode, dark_mode)
    'text_secondary': ('#555555', '#b0b0b0'),
    'card_bg':        ('#e8e8ee', '#2a2a3e'),
    'input_bg':       ('#d8d8e0', '#1e1e30'),
}
 
SENSOR_GROUPS = {
    'Spaceborne Multispectral': [
        {'id': 'aster', 'name': 'ASTER', 'icon': 'AS',
         'desc': 'Terra ASTER multispectral VNIR/SWIR/TIR with auto-scaling',
         'products': {
             'V004 GeoTIFF': {'sensor_id': 'aster_v004',
                              'desc': 'AST_05, AST_07, AST_08, AST_09T, AST_L1T'},
             'V003 HDF':     {'sensor_id': 'aster_v003',
                              'desc': 'AST_05, AST_07, AST_07M, AST_08, AST_09T'},
         }},
        {'id': 'landsat', 'name': 'Landsat', 'icon': 'LS',
         'desc': 'Collection 2 — sensor auto-detected from MTL metadata',
         'products': {
             'Level-1 (TOA)':     {'sensor_id': 'landsat_l1',
                                   'desc': 'LC08/09, LE07, LT04/05, LM01-05 (L1TP/L1GT)'},
             'Level-2 (Surface)': {'sensor_id': 'landsat_l2',
                                   'desc': 'LC08/09, LE07, LT04/05 only (L2SP — no MSS)'},
         }},
        {'id': 'sentinel2', 'name': 'Sentinel-2', 'icon': 'S2',
         'desc': 'Surface Reflectance (10m + 20m bands)',
         'products': {
             'L2A/L2B/L2C': {'sensor_id': 'sentinel2',
                             'desc': 'R10m + R20m Reflectance + Scene Classification'},
         }},
        {'id': 'sentinel3', 'name': 'Sentinel-3', 'icon': 'S3',
         'desc': 'SLSTR/OLCI/Synergy — auto-detected from .SEN3 folder',
         'products': {
             'SLSTR L1B RBT':    {'sensor_id': 's3_slstr_l1',
                                   'desc': 'S1-S6 Radiance (500m) + S7-S9/F1-F2 BT (1km)'},
             'SLSTR L2 LST':     {'sensor_id': 's3_slstr_l2',
                                   'desc': 'Land Surface Temperature (1km)'},
             'OLCI L1B EFR':     {'sensor_id': 's3_olci_l1',
                                   'desc': 'TOA Radiance, 21 bands (300m)'},
             'Synergy L2 SYN':   {'sensor_id': 's3_syn_l2',
                                   'desc': 'Surface Reflectance — OLCI + SLSTR (300m)'},
         }},
        {'id': 'ecostress', 'name': 'ECOSTRESS', 'icon': 'EC',
         'desc': 'ISS thermal infrared (5 TIR bands)',
         'products': {
             'L1CT Radiance': {'sensor_id': 'ecostress_l1ct',
                               'desc': 'Terrain-corrected at-sensor radiance (7 bands)'},
             'L2 LSTE':       {'sensor_id': 'ecostress_l2',
                               'desc': 'Land Surface Temperature and Emissivity'},
         }},
        {'id': 'viirs', 'name': 'VIIRS', 'icon': 'VI',
         'desc': 'Suomi NPP / NOAA-20 multispectral imager',
         'products': {
             'VNP02IMG / VNP02MOD / VNP21': {'sensor_id': 'viirs',
                                              'desc': 'Calibrated Radiance + LST/Emissivity'},
         }},
        {'id': 'modis', 'name': 'MODIS', 'icon': 'MO',
         'desc': 'Terra/Aqua gridded products (sinusoidal projection)',
         'products': {
             'Surface Reflectance': {'sensor_id': 'modis_refl',
                                     'desc': 'MOD09GA / MYD09GA — 7 bands, 500m'},
             'LST / Emissivity':    {'sensor_id': 'modis_lst',
                                     'desc': 'MOD11A1 / MYD11A1 — LST + Emis, 1km'},
         }},
    ],
    'Spaceborne Hyperspectral': [
        {'id': 'emit', 'name': 'EMIT', 'icon': 'EM',
         'desc': 'ISS imaging spectrometer (285 bands, 380–2500 nm)',
         'products': {
             'L2A Reflectance': {'sensor_id': 'emit',
                                 'desc': 'Surface Reflectance + Mask + RMSE Uncertainty'},
         }},
        {'id': 'airs', 'name': 'AIRS', 'icon': 'AI',
         'desc': 'Aqua atmospheric sounder (2378 channels, 3.7–15.4 µm)',
         'products': {
             'L1B Radiance': {'sensor_id': 'airs',
                              'desc': 'Calibrated Radiance Cube'},
         }},
    ],
    'Airborne': [
        {'id': 'master', 'name': 'MASTER', 'icon': 'MA',
         'desc': 'MODIS/ASTER Airborne Simulator (50 bands, VSWIR + TIR)',
         'products': {
             'L1B Radiance':      {'sensor_id': 'master_l1b',
                                   'desc': 'Calibrated Radiance (VSWIR + TIR)'},
             'L2 Emissivity/LST': {'sensor_id': 'master_l2',
                                   'desc': 'Emissivity and Land Surface Temperature'},
         }},
        {'id': 'hytes', 'name': 'HyTES', 'icon': 'HT',
         'desc': 'Hyperspectral Thermal Emission Spectrometer (256 TIR bands)',
         'products': {
             'L1 Radiance':       {'sensor_id': 'hytes_l1',
                                   'desc': 'TIR Radiance (7.5–12 µm)'},
             'L2 Emissivity/LST': {'sensor_id': 'hytes_l2',
                                   'desc': 'Emissivity + LST + PC Emissivity'},
         }},
        {'id': 'aviris', 'name': 'AVIRIS', 'icon': 'AV',
         'desc': 'Airborne VSWIR Imaging Spectrometer family',
         'products': {
             'AVIRIS-3/5 L1B':  {'sensor_id': 'aviris35',
                                 'desc': 'AVIRIS-3 or AVIRIS-5 L1B Radiance (NetCDF)'},
             'AVIRIS-3/5 L2A':  {'sensor_id': 'aviris35_l2a',
                                 'desc': 'AVIRIS-3/5 L2A Surface Reflectance (NetCDF ORT)'},
             'NG L1B Radiance': {'sensor_id': 'aviris_ng_l1b',
                                 'desc': 'AVIRIS-NG orthorectified radiance (ENVI)'},
             'NG L2 Reflectance': {'sensor_id': 'aviris_ng_l2',
                                   'desc': 'AVIRIS-NG corrected reflectance (ENVI)'},
             'Classic L1B Radiance': {'sensor_id': 'aviris_classic_l1b',
                                      'desc': 'Classic AVIRIS orthorectified radiance (ENVI)'},
             'Classic L2 Reflectance': {'sensor_id': 'aviris_classic_l2',
                                        'desc': 'Classic AVIRIS surface reflectance (ENVI)'},
         }},
    ],
}

# Sensors that support spatial subsetting (bbox)
BBOX_SENSORS = {
    'ecostress_l1ct', 'ecostress_l2', 'emit', 'aviris35', 'aviris35_l2a',
    'aviris_ng_l1b', 'aviris_ng_l2', 'aviris_classic_l1b', 'aviris_classic_l2',
    'master_l1b', 'master_l2', 'hytes_l1', 'hytes_l2', 'airs',
    's3_slstr_l1', 's3_slstr_l2', 's3_olci_l1', 's3_syn_l2',
}

# Sensors that support spectral subsetting
SPECTRAL_SENSORS = {
    'emit', 'aviris35', 'aviris35_l2a',
    'aviris_ng_l1b', 'aviris_ng_l2', 'aviris_classic_l1b', 'aviris_classic_l2',
    'master_l1b', 'hytes_l1', 'hytes_l2', 'airs',
}

 
# ============================================================================
# Logging system (thread-safe)
# ============================================================================
 
class LogQueue:
    def __init__(self):
        self.queue = queue.Queue()
 
    def log(self, message, level='info'):
        timestamp = datetime.now().strftime('%H:%M:%S')
        self.queue.put((timestamp, level, message))
 
    def get_all(self):
        messages = []
        while not self.queue.empty():
            try:
                messages.append(self.queue.get_nowait())
            except queue.Empty:
                break
        return messages
 
 
# ============================================================================
# Sidebar Navigation
# ============================================================================
 
class SidebarButton(ctk.CTkButton):
    def __init__(self, master, text, icon='', command=None, **kwargs):
        super().__init__(
            master,
            text=f"  {icon}  {text}" if icon else f"  {text}",
            command=command,
            anchor='w',
            height=45,
            corner_radius=8,
            font=ctk.CTkFont(size=14),
            fg_color='transparent',
            text_color='#b0b0b0',
            hover_color=COLORS['sidebar_hover'],
            **kwargs
        )
 
    def set_active(self, active=True):
        if active:
            self.configure(fg_color=COLORS['sidebar_active'],
                           text_color='#ffffff')
        else:
            self.configure(fg_color='transparent',
                           text_color='#b0b0b0')
 
 
# ============================================================================
# Sensor Card Widget
# ============================================================================
 
class SensorCard(ctk.CTkFrame):
    def __init__(self, master, sensor_info, on_select, **kwargs):
        super().__init__(master, corner_radius=12, fg_color=COLORS['card_bg'],
                         cursor='hand2', **kwargs)
 
        self.sensor_info = sensor_info
        self.on_select   = on_select
        self.is_selected = False
 
        self.grid_columnconfigure(0, minsize=50, weight=0)
        self.grid_columnconfigure(1, weight=1)
 
        self.icon_label = ctk.CTkLabel(self, text=sensor_info['icon'],
                     font=ctk.CTkFont(size=14, weight='bold'),
                     width=42, height=42, corner_radius=8,
                     fg_color=COLORS['accent'], text_color='#ffffff')
        self.icon_label.grid(row=0, column=0, rowspan=2, padx=(12, 4), pady=10)
 
        self.name_label = ctk.CTkLabel(self, text=sensor_info['name'],
                     font=ctk.CTkFont(size=15, weight='bold'), anchor='w')
        self.name_label.grid(row=0, column=1, padx=10, pady=(10, 0), sticky='w')
 
        self.desc_label = ctk.CTkLabel(self, text=sensor_info['desc'],
                     font=ctk.CTkFont(size=12),
                     text_color=COLORS['text_secondary'], anchor='w')
        self.desc_label.grid(row=1, column=1, padx=10, pady=(0, 10), sticky='w')
 
        self.bind('<Button-1>', lambda e: self._on_click())
        for child in self.winfo_children():
            child.bind('<Button-1>', lambda e: self._on_click())
 
    def _on_click(self):
        self.on_select(self.sensor_info)
 
    def set_selected(self, selected):
        self.is_selected = selected
        if selected:
            self.configure(fg_color=COLORS['accent'],
                           border_color=COLORS['accent'], border_width=2)
            self.name_label.configure(text_color='#ffffff')
            self.desc_label.configure(text_color='#ffffff')
        else:
            self.configure(fg_color=COLORS['card_bg'], border_width=0)
            self.name_label.configure(text_color=COLORS['text_primary'])
            self.desc_label.configure(text_color=COLORS['text_secondary'])
 
 
# ============================================================================
# Processing Panel
# ============================================================================
 
class ProcessingPanel(ctk.CTkFrame):
    def __init__(self, master, log_queue, **kwargs):
        super().__init__(master, fg_color='transparent', **kwargs)
 
        self.log_queue      = log_queue
        self.current_sensor = None
        self.processing     = False
 
        self.grid_columnconfigure(0, weight=1)
        # Top config area gets more space, log is compact but scrollable
        self.grid_rowconfigure(0, weight=3)   # scrollable config area
        self.grid_rowconfigure(1, weight=1)   # log section — compact

        # ================================================================
        # TOP: Scrollable config area (header + I/O + options + subsetting)
        # ================================================================
        self.config_scroll = ctk.CTkScrollableFrame(self, fg_color='transparent')
        self.config_scroll.grid(row=0, column=0, sticky='nsew', padx=0, pady=0)
        self.config_scroll.grid_columnconfigure(0, weight=1)

        # --- Sensor header ---
        self.header_frame = ctk.CTkFrame(self.config_scroll, corner_radius=12,
                                          fg_color=COLORS['card_bg'])
        self.header_frame.grid(row=0, column=0, sticky='ew', padx=10, pady=(10, 5))
        self.header_frame.grid_columnconfigure(0, weight=1)
 
        self.sensor_title = ctk.CTkLabel(
            self.header_frame, text="Select a sensor to begin",
            font=ctk.CTkFont(size=20, weight='bold'), anchor='w')
        self.sensor_title.grid(row=0, column=0, padx=20, pady=(15, 5), sticky='w')
 
        self.sensor_desc = ctk.CTkLabel(
            self.header_frame, text="Choose a sensor from the left panel",
            font=ctk.CTkFont(size=13),
            text_color=COLORS['text_secondary'], anchor='w')
        self.sensor_desc.grid(row=1, column=0, padx=20, pady=(0, 15), sticky='w')
 
        # --- I/O section ---
        self.io_frame = ctk.CTkFrame(self.config_scroll, corner_radius=12,
                                      fg_color=COLORS['card_bg'])
        self.io_frame.grid(row=1, column=0, sticky='ew', padx=10, pady=5)
        self.io_frame.grid_columnconfigure(1, weight=1)
 
        ctk.CTkLabel(self.io_frame, text="Input:",
                     font=ctk.CTkFont(size=13, weight='bold')
                     ).grid(row=0, column=0, padx=(20, 10), pady=(15, 5), sticky='w')
        self.input_entry = ctk.CTkEntry(
            self.io_frame, placeholder_text="Select input folder...",
            font=ctk.CTkFont(size=13), height=36)
        self.input_entry.grid(row=0, column=1, padx=5, pady=(15, 5), sticky='ew')
        ctk.CTkButton(self.io_frame, text="Browse", width=80, height=36,
                      command=self._browse_input, corner_radius=8
                      ).grid(row=0, column=2, padx=(5, 20), pady=(15, 5))
 
        ctk.CTkLabel(self.io_frame, text="Output:",
                     font=ctk.CTkFont(size=13, weight='bold')
                     ).grid(row=1, column=0, padx=(20, 10), pady=(5, 15), sticky='w')
        self.output_entry = ctk.CTkEntry(
            self.io_frame,
            placeholder_text="Select output folder (default: input/converted)...",
            font=ctk.CTkFont(size=13), height=36)
        self.output_entry.grid(row=1, column=1, padx=5, pady=(5, 15), sticky='ew')
        ctk.CTkButton(self.io_frame, text="Browse", width=80, height=36,
                      command=self._browse_output, corner_radius=8
                      ).grid(row=1, column=2, padx=(5, 20), pady=(5, 15))
 
        # --- Options section (format + sensor-specific) ---
        self.options_frame = ctk.CTkFrame(self.config_scroll, corner_radius=12,
                                           fg_color=COLORS['card_bg'])
        self.options_frame.grid(row=2, column=0, sticky='ew', padx=10, pady=5)
        self.options_frame.grid_columnconfigure((0, 1, 2, 3), weight=1)
 
        ctk.CTkLabel(self.options_frame, text="Output Format:",
                     font=ctk.CTkFont(size=13, weight='bold')
                     ).grid(row=0, column=0, padx=20, pady=(15, 5), sticky='w')
        self.format_var  = ctk.StringVar(value='both')
        self.format_menu = ctk.CTkSegmentedButton(
            self.options_frame, values=['GeoTIFF', 'ENVI', 'Both'],
            variable=self.format_var, font=ctk.CTkFont(size=12))
        self.format_menu.grid(row=0, column=1, columnspan=2,
                               padx=10, pady=(15, 5), sticky='ew')
        self.format_menu.set('Both')

        # Overwrite checkbox
        self.overwrite_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(self.options_frame, text="Overwrite existing files",
                        variable=self.overwrite_var,
                        font=ctk.CTkFont(size=12)
                        ).grid(row=0, column=3, padx=(10, 20), pady=(15, 5), sticky='w')
 
        self.sensor_options_frame = ctk.CTkFrame(self.options_frame,
                                                   fg_color='transparent')
        self.sensor_options_frame.grid(row=1, column=0, columnspan=4,
                                        padx=20, pady=(5, 15), sticky='ew')

        # --- Subsetting / Spectral section (shared, shown conditionally) ---
        self.subset_spectral_frame = ctk.CTkFrame(self.config_scroll, corner_radius=12,
                                                    fg_color=COLORS['card_bg'])
        # Row 3 in config_scroll — shown/hidden depending on sensor
        self.subset_spectral_frame.grid_columnconfigure((0, 1, 2, 3), weight=1)
        self._build_subset_spectral_panel()

        # ================================================================
        # BOTTOM: Fixed log + buttons (always visible, never squeezed)
        # ================================================================

        # --- Log section ---
        self.log_frame = ctk.CTkFrame(self, corner_radius=12,
                                       fg_color=COLORS['card_bg'])
        self.log_frame.grid(row=1, column=0, sticky='nsew', padx=10, pady=5)
        self.log_frame.grid_columnconfigure(0, weight=1)
        self.log_frame.grid_rowconfigure(1, weight=1)
 
        log_header = ctk.CTkFrame(self.log_frame, fg_color='transparent')
        log_header.grid(row=0, column=0, sticky='ew', padx=15, pady=(10, 0))
        log_header.grid_columnconfigure(0, weight=1)
 
        ctk.CTkLabel(log_header, text="Processing Log",
                     font=ctk.CTkFont(size=14, weight='bold')
                     ).grid(row=0, column=0, sticky='w')
        ctk.CTkButton(log_header, text="Clear", width=60, height=28,
                      font=ctk.CTkFont(size=11), corner_radius=6,
                      fg_color='transparent', border_width=1,
                      command=self._clear_log
                      ).grid(row=0, column=1)
 
        self.log_text = ctk.CTkTextbox(
            self.log_frame, font=ctk.CTkFont(family='Consolas', size=12),
            fg_color=COLORS['input_bg'], corner_radius=8, wrap='word',
            height=100)
        self.log_text.grid(row=1, column=0, sticky='nsew', padx=15, pady=(5, 10))
 
        # --- Process + Cancel buttons + progress ---
        self.btn_frame = ctk.CTkFrame(self, fg_color='transparent')
        self.btn_frame.grid(row=2, column=0, sticky='ew', padx=10, pady=10)
        self.btn_frame.grid_columnconfigure(0, weight=1)

        self.process_btn = ctk.CTkButton(
            self.btn_frame, text="▶  Process Data", height=48,
            font=ctk.CTkFont(size=16, weight='bold'),
            corner_radius=10, fg_color=COLORS['accent'],
            hover_color=COLORS['accent_hover'],
            command=self._start_processing)
        self.process_btn.grid(row=0, column=0, sticky='ew', padx=(0, 5))

        self.cancel_btn = ctk.CTkButton(
            self.btn_frame, text="■  Cancel", height=48, width=120,
            font=ctk.CTkFont(size=14, weight='bold'),
            corner_radius=10, fg_color=COLORS['error'],
            hover_color='#c0392b', state='disabled',
            command=self._cancel_processing)
        self.cancel_btn.grid(row=0, column=1, padx=(5, 5))

        self.open_folder_btn = ctk.CTkButton(
            self.btn_frame, text="📂  Open File Location", height=48, width=180,
            font=ctk.CTkFont(size=14, weight='bold'),
            corner_radius=10, fg_color=COLORS['sidebar_active'],
            hover_color=COLORS['sidebar_hover'], state='disabled',
            command=self._open_output_folder)
        self.open_folder_btn.grid(row=0, column=2, padx=(5, 0))

        self._cancel_flag = False
 
        self.progress = ctk.CTkProgressBar(self, height=4, corner_radius=2)
        self.progress.grid(row=3, column=0, sticky='ew', padx=10, pady=(0, 10))
        self.progress.set(0)
 
        self._poll_log()

    # ------------------------------------------------------------------
    # Subsetting / Spectral panel (shared across sensors)
    # ------------------------------------------------------------------

    def _build_subset_spectral_panel(self):
        """Build the shared subsetting and spectral options panel."""
        frame = self.subset_spectral_frame

        # --- Spatial subsetting row ---
        self.bbox_enable_var = ctk.BooleanVar(value=False)
        self.bbox_check = ctk.CTkCheckBox(
            frame, text="Spatial Subset (Bounding Box)",
            variable=self.bbox_enable_var,
            font=ctk.CTkFont(size=13, weight='bold'),
            command=self._toggle_bbox)
        self.bbox_check.grid(row=0, column=0, columnspan=2,
                             padx=20, pady=(15, 5), sticky='w')

        self.bbox_fields_frame = ctk.CTkFrame(frame, fg_color='transparent')
        self.bbox_fields_frame.grid(row=1, column=0, columnspan=4,
                                     padx=30, pady=(0, 5), sticky='ew')
        self.bbox_fields_frame.grid_columnconfigure((1, 3, 5, 7), weight=1)

        labels = ['Min Lat:', 'Max Lat:', 'Min Lon:', 'Max Lon:']
        self.bbox_entries = {}
        for i, lbl in enumerate(labels):
            key = lbl.replace(':', '').replace(' ', '_').lower()
            ctk.CTkLabel(self.bbox_fields_frame, text=lbl,
                         font=ctk.CTkFont(size=11)
                         ).grid(row=0, column=i*2, padx=(5, 2), pady=5, sticky='e')
            entry = ctk.CTkEntry(self.bbox_fields_frame, width=90, height=28,
                                  font=ctk.CTkFont(size=11),
                                  placeholder_text='0.0', state='disabled')
            entry.grid(row=0, column=i*2+1, padx=(2, 5), pady=5, sticky='w')
            self.bbox_entries[key] = entry

        # Subset file browse
        self.subset_file_var = ctk.StringVar(value='')
        ctk.CTkLabel(self.bbox_fields_frame, text="Or file:",
                     font=ctk.CTkFont(size=11)
                     ).grid(row=1, column=0, padx=(5, 2), pady=5, sticky='e')
        self.subset_file_entry = ctk.CTkEntry(
            self.bbox_fields_frame, textvariable=self.subset_file_var,
            placeholder_text="Shapefile, KML, GeoJSON...",
            font=ctk.CTkFont(size=11), height=28, state='disabled')
        self.subset_file_entry.grid(row=1, column=1, columnspan=5,
                                     padx=2, pady=5, sticky='ew')
        self.subset_file_btn = ctk.CTkButton(
            self.bbox_fields_frame, text="Browse", width=60, height=28,
            font=ctk.CTkFont(size=11), state='disabled',
            command=self._browse_subset_file)
        self.subset_file_btn.grid(row=1, column=6, columnspan=2,
                                   padx=(2, 5), pady=5, sticky='w')

        # --- Spectral options row ---
        sep = ctk.CTkFrame(frame, height=1, fg_color=COLORS['text_secondary'])
        sep.grid(row=2, column=0, columnspan=4, sticky='ew', padx=20, pady=5)

        self.spectral_label = ctk.CTkLabel(
            frame, text="Spectral Options",
            font=ctk.CTkFont(size=13, weight='bold'))
        self.spectral_label.grid(row=3, column=0, padx=20, pady=(5, 5), sticky='w')

        self.spectral_inner = ctk.CTkFrame(frame, fg_color='transparent')
        self.spectral_inner.grid(row=4, column=0, columnspan=4,
                                  padx=30, pady=(0, 15), sticky='ew')
        self.spectral_inner.grid_columnconfigure((1, 3), weight=1)

        ctk.CTkLabel(self.spectral_inner, text="Mode:",
                     font=ctk.CTkFont(size=11)
                     ).grid(row=0, column=0, padx=(0, 5), pady=5, sticky='e')
        self.spectral_mode_var = ctk.StringVar(value='Default')
        self.spectral_mode_menu = ctk.CTkOptionMenu(
            self.spectral_inner,
            values=['Default', 'Region', 'Match Sensor', 'Every Nth', 'Custom Range'],
            variable=self.spectral_mode_var,
            font=ctk.CTkFont(size=11), width=150,
            command=self._on_spectral_mode_change)
        self.spectral_mode_menu.grid(row=0, column=1, padx=5, pady=5, sticky='w')

        # Spectral value field (context-dependent)
        self.spectral_value_label = ctk.CTkLabel(
            self.spectral_inner, text="",
            font=ctk.CTkFont(size=11))
        self.spectral_value_label.grid(row=0, column=2, padx=(15, 5), pady=5, sticky='e')

        self.spectral_value_var = ctk.StringVar(value='')
        self.spectral_value_widget = ctk.CTkOptionMenu(
            self.spectral_inner, values=['VNIR', 'SWIR', 'VIS', 'NIR', 'Full'],
            variable=self.spectral_value_var,
            font=ctk.CTkFont(size=11), width=150)
        self.spectral_value_widget.grid(row=0, column=3, padx=5, pady=5, sticky='w')
        self.spectral_value_widget.grid_remove()  # hidden by default

        # Entry version for custom range / every-nth
        self.spectral_value_entry = ctk.CTkEntry(
            self.spectral_inner, font=ctk.CTkFont(size=11), width=150, height=28,
            placeholder_text='')
        self.spectral_value_entry.grid(row=0, column=3, padx=5, pady=5, sticky='w')
        self.spectral_value_entry.grid_remove()

    def _toggle_bbox(self):
        """Enable/disable bbox entry fields."""
        state = 'normal' if self.bbox_enable_var.get() else 'disabled'
        for entry in self.bbox_entries.values():
            entry.configure(state=state)
        self.subset_file_entry.configure(state=state)
        self.subset_file_btn.configure(state=state)

    def _browse_subset_file(self):
        path = filedialog.askopenfilename(
            title="Select Subset File",
            filetypes=[
                ("Vector/Text files", "*.shp *.kml *.kmz *.geojson *.json *.txt *.csv"),
                ("All files", "*.*"),
            ])
        if path:
            self.subset_file_var.set(path)

    def _on_spectral_mode_change(self, mode):
        """Show/hide the appropriate spectral value widget."""
        self.spectral_value_widget.grid_remove()
        self.spectral_value_entry.grid_remove()
        self.spectral_value_label.configure(text='')

        if mode == 'Region':
            self.spectral_value_label.configure(text='Region:')
            self.spectral_value_widget.configure(
                values=['VNIR', 'SWIR', 'VIS', 'NIR', 'Full'])
            self.spectral_value_var.set('VNIR')
            self.spectral_value_widget.grid()
        elif mode == 'Match Sensor':
            self.spectral_value_label.configure(text='Sensor:')
            self.spectral_value_widget.configure(
                values=['Landsat 8/9', 'Sentinel-2', 'ASTER', 'MODIS'])
            self.spectral_value_var.set('Landsat 8/9')
            self.spectral_value_widget.grid()
        elif mode == 'Every Nth':
            self.spectral_value_label.configure(text='N:')
            self.spectral_value_entry.configure(placeholder_text='5')
            self.spectral_value_entry.grid()
        elif mode == 'Custom Range':
            self.spectral_value_label.configure(text='Ranges (nm):')
            self.spectral_value_entry.configure(placeholder_text='450-900,2000-2400')
            self.spectral_value_entry.grid()

    def _update_subset_spectral_visibility(self, sensor_id):
        """Show or hide the subsetting/spectral panel based on sensor."""
        show_bbox = sensor_id in BBOX_SENSORS
        show_spectral = sensor_id in SPECTRAL_SENSORS

        if show_bbox or show_spectral:
            self.subset_spectral_frame.grid(row=3, column=0, sticky='ew',
                                             padx=10, pady=5)
        else:
            self.subset_spectral_frame.grid_forget()

        # Show/hide bbox row
        if show_bbox:
            self.bbox_check.grid()
            self.bbox_fields_frame.grid()
        else:
            self.bbox_check.grid_remove()
            self.bbox_fields_frame.grid_remove()
            self.bbox_enable_var.set(False)

        # Show/hide spectral row
        if show_spectral:
            self.spectral_label.grid()
            self.spectral_inner.grid()
        else:
            self.spectral_label.grid_remove()
            self.spectral_inner.grid_remove()
            self.spectral_mode_var.set('Default')

    # ------------------------------------------------------------------
    # Helper: build bbox/spectral kwargs from GUI state
    # ------------------------------------------------------------------

    def _get_bbox_tuple(self):
        """Return (min_lat, max_lat, min_lon, max_lon) or None."""
        if not self.bbox_enable_var.get():
            return None

        # Check for subset file first
        sfile = self.subset_file_var.get().strip()
        if sfile and os.path.isfile(sfile):
            # Return as subset_opts dict for ECOSTRESS, or convert for others
            return ('file', sfile)

        try:
            vals = {k: float(e.get()) for k, e in self.bbox_entries.items()
                    if e.get().strip()}
            if len(vals) == 4:
                return (vals['min_lat'], vals['max_lat'],
                        vals['min_lon'], vals['max_lon'])
        except ValueError:
            pass
        return None

    def _get_spectral_opts(self):
        """Build spectral_opts dict from GUI state, or None."""
        mode = self.spectral_mode_var.get()
        if mode == 'Default':
            return None

        sensor_map = {
            'Landsat 8/9': 'landsat89',
            'Sentinel-2':  'sentinel2',
            'ASTER':       'aster',
            'MODIS':       'modis_land',
        }

        if mode == 'Region':
            return {'mode': 'region', 'region': self.spectral_value_var.get().lower(),
                    'exclude_water': True}
        elif mode == 'Match Sensor':
            sname = self.spectral_value_var.get()
            return {'mode': 'sensor', 'sensor': sensor_map.get(sname, 'landsat89'),
                    'nearest': False, 'rsr_dir': None, 'exclude_water': True}
        elif mode == 'Every Nth':
            try:
                n = int(self.spectral_value_entry.get().strip())
            except ValueError:
                n = 5
            return {'mode': 'nth', 'every_nth': n, 'exclude_water': True}
        elif mode == 'Custom Range':
            ranges_str = self.spectral_value_entry.get().strip()
            return {'mode': 'bands', 'ranges': ranges_str, 'exclude_water': True}
        return None

    # ------------------------------------------------------------------
    # Sensor options
    # ------------------------------------------------------------------
 
    def set_sensor(self, sensor_info):
        self.current_sensor = sensor_info
        self.sensor_title.configure(
            text=f"{sensor_info['icon']}  {sensor_info['name']}")
        self.sensor_desc.configure(text=sensor_info['desc'])

        # Clear existing sensor options
        for widget in self.sensor_options_frame.winfo_children():
            widget.destroy()

        # Build product dropdown if sensor has multiple products
        products = sensor_info.get('products', {})
        product_names = list(products.keys())

        if len(product_names) > 1:
            # Multi-product sensor — show dropdown
            prod_frame = ctk.CTkFrame(self.sensor_options_frame, fg_color='transparent')
            prod_frame.grid(row=0, column=0, columnspan=4, sticky='ew')
            prod_frame.grid_columnconfigure(1, weight=1)

            ctk.CTkLabel(prod_frame, text="Product:",
                         font=ctk.CTkFont(size=13, weight='bold')
                         ).grid(row=0, column=0, padx=(0, 10), pady=5, sticky='e')

            self.product_var = ctk.StringVar(value=product_names[0])
            self.product_menu = ctk.CTkOptionMenu(
                prod_frame, values=product_names,
                variable=self.product_var,
                font=ctk.CTkFont(size=12), width=280,
                command=self._on_product_change)
            self.product_menu.grid(row=0, column=1, padx=5, pady=5, sticky='w')

            # Product description label
            self.product_desc_label = ctk.CTkLabel(
                prod_frame, text=products[product_names[0]]['desc'],
                font=ctk.CTkFont(size=11),
                text_color=COLORS['text_secondary'])
            self.product_desc_label.grid(row=0, column=2, padx=(15, 0),
                                          pady=5, sticky='w')
        else:
            self.product_var = ctk.StringVar(
                value=product_names[0] if product_names else '')
            self.product_desc_label = None

        # Get the active product sensor_id and build its options
        self._active_sensor_id = self._get_active_sensor_id()
        self._build_sensor_options(self._active_sensor_id)
        self._update_subset_spectral_visibility(self._active_sensor_id)

    def _on_product_change(self, product_name):
        """Handle product dropdown change — rebuild sensor-specific options."""
        products = self.current_sensor.get('products', {})
        prod_info = products.get(product_name, {})

        # Update description label
        if self.product_desc_label is not None:
            self.product_desc_label.configure(text=prod_info.get('desc', ''))

        # Clear and rebuild sensor-specific options (keep product dropdown)
        # Find and destroy everything except the product frame (row 0)
        children = self.sensor_options_frame.winfo_children()
        for widget in children:
            grid_info = widget.grid_info()
            if grid_info.get('row', 0) > 0:
                widget.destroy()

        self._active_sensor_id = self._get_active_sensor_id()
        self._build_sensor_options(self._active_sensor_id)
        self._update_subset_spectral_visibility(self._active_sensor_id)

    def _get_active_sensor_id(self):
        """Get the sensor_id for the currently selected product."""
        if self.current_sensor is None:
            return None
        products = self.current_sensor.get('products', {})
        product_name = self.product_var.get()
        prod_info = products.get(product_name, {})
        return prod_info.get('sensor_id', self.current_sensor['id'])
 
    def _build_sensor_options(self, sensor_id):
        frame = self.sensor_options_frame
        frame.grid_columnconfigure((0, 1, 2, 3), weight=1)

        # Determine start row (1 if product dropdown exists in row 0, else 0)
        products = self.current_sensor.get('products', {})
        start_row = 1 if len(products) > 1 else 0

        if sensor_id in ('aster_v004', 'aster_v003'):
            ctk.CTkLabel(frame, text="ASTER Product:",
                         font=ctk.CTkFont(size=12)
                         ).grid(row=start_row, column=0, padx=(0, 5), pady=5, sticky='e')
            self.aster_product_var = ctk.StringVar(value='AST_05')
            if sensor_id == 'aster_v003':
                product_values = ['AST_05', 'AST_07', 'AST_07M', 'AST_08', 'AST_09T']
            else:
                product_values = ['AST_05', 'AST_07', 'AST_08', 'AST_09T', 'AST_L1T']
            self.aster_product_menu = ctk.CTkOptionMenu(
                frame, values=product_values,
                variable=self.aster_product_var,
                font=ctk.CTkFont(size=12), width=140)
            self.aster_product_menu.grid(row=start_row, column=1, padx=5, pady=5, sticky='w')

            self.aster_batch_var = ctk.BooleanVar(value=True)
            ctk.CTkCheckBox(frame, text="Batch mode (subfolder per granule)",
                            variable=self.aster_batch_var,
                            font=ctk.CTkFont(size=12)
                            ).grid(row=start_row, column=2, columnspan=2,
                                   padx=10, pady=5, sticky='w')

        elif sensor_id in ('aviris35', 'aviris_ng_l1b', 'aviris_ng_l2',
                           'aviris_classic_l1b', 'aviris_classic_l2'):
            self.av3_water_var = ctk.BooleanVar(value=True)
            ctk.CTkCheckBox(frame, text="Exclude water vapor bands",
                            variable=self.av3_water_var,
                            font=ctk.CTkFont(size=12)
                            ).grid(row=start_row, column=0, columnspan=2,
                                   padx=10, pady=5, sticky='w')
            if sensor_id == 'aviris35':
                ctk.CTkLabel(frame, text="Resolution (m):",
                             font=ctk.CTkFont(size=12)
                             ).grid(row=start_row, column=2, padx=(10, 5), pady=5, sticky='e')
                self.resolution_entry = ctk.CTkEntry(
                    frame, placeholder_text="auto (from GLT)", width=120, height=30,
                    font=ctk.CTkFont(size=12))
                self.resolution_entry.grid(row=start_row, column=3, padx=5, pady=5, sticky='w')
            if sensor_id in ('aviris_ng_l1b', 'aviris_ng_l2'):
                self.av_bbl_var = ctk.BooleanVar(value=True)
                ctk.CTkCheckBox(frame, text="Use bad bands list (BBL)",
                                variable=self.av_bbl_var,
                                font=ctk.CTkFont(size=12)
                                ).grid(row=start_row+1, column=0, columnspan=2,
                                       padx=10, pady=5, sticky='w')

        elif sensor_id == 'landsat_l1':
            ctk.CTkLabel(
                frame,
                text="Bands 1-7,9 → TOA Reflectance  |  Bands 10-11 → Radiance  |  Band 8 → Pan Reflectance",
                font=ctk.CTkFont(size=12), text_color=COLORS['text_secondary']
            ).grid(row=start_row, column=0, columnspan=4, pady=5)

        elif sensor_id == 'landsat_l2':
            ctk.CTkLabel(
                frame,
                text="SR Bands 1-7 → Surface Reflectance  |  ST Band 10 → Temperature (K)  |  ST Ancillary layers",
                font=ctk.CTkFont(size=12), text_color=COLORS['text_secondary']
            ).grid(row=start_row, column=0, columnspan=4, pady=5)

        elif sensor_id == 'sentinel2':
            ctk.CTkLabel(
                frame,
                text="R10m: B02,B03,B04,B08  |  R20m: B01-B07,B8A,B11,B12  |  Scene Classification",
                font=ctk.CTkFont(size=12), text_color=COLORS['text_secondary']
            ).grid(row=start_row, column=0, columnspan=4, pady=5)

        elif sensor_id == 'modis_refl':
            ctk.CTkLabel(
                frame,
                text="Bands 1-7 → Surface Reflectance (500m, sinusoidal projection)",
                font=ctk.CTkFont(size=12), text_color=COLORS['text_secondary']
            ).grid(row=start_row, column=0, columnspan=4, pady=5)

        elif sensor_id == 'modis_lst':
            ctk.CTkLabel(
                frame,
                text="LST Day/Night (K) + Emis 31/32 (1km, sinusoidal projection)",
                font=ctk.CTkFont(size=12), text_color=COLORS['text_secondary']
            ).grid(row=start_row, column=0, columnspan=4, pady=5)

        elif sensor_id in ('ecostress_l2', 'master_l1b', 'hytes_l1', 'hytes_l2', 'airs'):
            ctk.CTkLabel(frame, text="Projection:",
                         font=ctk.CTkFont(size=12)
                         ).grid(row=start_row, column=0, padx=(0, 5), pady=5, sticky='e')
            self.projection_var = ctk.StringVar(value='geographic')
            proj_btn = ctk.CTkSegmentedButton(
                frame, values=['Geographic', 'UTM'],
                variable=self.projection_var, font=ctk.CTkFont(size=11))
            proj_btn.grid(row=start_row, column=1, padx=5, pady=5, sticky='w')
            proj_btn.set('Geographic')

        if sensor_id == 'viirs':
            ctk.CTkLabel(frame, text="Bands:",
                         font=ctk.CTkFont(size=12)
                         ).grid(row=start_row, column=0, padx=(0, 5), pady=5, sticky='e')
            self.viirs_bands_var = ctk.StringVar(value='All bands')
            ctk.CTkOptionMenu(
                frame,
                values=['All bands', 'TIR only (M14-M16)', 'Split-window (M15, M16)',
                        'MIR + TIR (M12-M16)', 'VSWIR only (M01-M11)',
                        'I-band TIR (I04, I05)', 'VNP21 core (LST + Emissivity)', 'Custom'],
                variable=self.viirs_bands_var,
                font=ctk.CTkFont(size=12), width=220,
                command=self._on_viirs_bands_change
            ).grid(row=start_row, column=1, padx=5, pady=5, sticky='w')
            self.viirs_custom_label = ctk.CTkLabel(frame, text="Bands:", font=ctk.CTkFont(size=11))
            self.viirs_custom_label.grid(row=start_row, column=2, padx=(10, 5), pady=5, sticky='e')
            self.viirs_custom_entry = ctk.CTkEntry(
                frame, placeholder_text="M14, M15, M16", font=ctk.CTkFont(size=11), width=180, height=28)
            self.viirs_custom_entry.grid(row=start_row, column=3, padx=5, pady=5, sticky='w')
            self.viirs_custom_label.grid_remove()
            self.viirs_custom_entry.grid_remove()

        if sensor_id == 'master_l1b':
            self.reflectance_var = ctk.BooleanVar(value=False)
            ctk.CTkCheckBox(frame, text="Convert VSWIR to TOA Reflectance",
                            variable=self.reflectance_var, font=ctk.CTkFont(size=12)
                            ).grid(row=start_row, column=2, columnspan=2, padx=10, pady=5, sticky='w')

        if sensor_id in ('master_l1b', 'master_l2', 'hytes_l1', 'hytes_l2'):
            ctk.CTkLabel(frame, text="Resolution (m):", font=ctk.CTkFont(size=12)
                         ).grid(row=start_row+1, column=0, padx=(0, 5), pady=5, sticky='e')
            self.resolution_entry = ctk.CTkEntry(
                frame, placeholder_text="auto", width=80, height=30, font=ctk.CTkFont(size=12))
            self.resolution_entry.grid(row=start_row+1, column=1, padx=5, pady=5, sticky='w')

        if sensor_id == 'emit':
            ctk.CTkLabel(frame, text="Bands:", font=ctk.CTkFont(size=12)
                         ).grid(row=start_row, column=0, padx=(0, 5), pady=5, sticky='e')
            self.emit_bands_var = ctk.StringVar(value='Good wavelengths')
            ctk.CTkOptionMenu(
                frame, values=['Good wavelengths', 'All bands', 'VNIR only', 'SWIR only',
                               'Match Landsat', 'Match Sentinel-2', 'Match ASTER', 'Match MODIS'],
                variable=self.emit_bands_var, font=ctk.CTkFont(size=12), width=180
            ).grid(row=start_row, column=1, padx=5, pady=5, sticky='w')
            self.emit_rmse_var = ctk.BooleanVar(value=True)
            ctk.CTkCheckBox(frame, text="RMSE Uncertainty", variable=self.emit_rmse_var,
                            font=ctk.CTkFont(size=12)
                            ).grid(row=start_row, column=2, padx=10, pady=5, sticky='w')

        if sensor_id == 'hytes_l2':
            self.pc_emis_var = ctk.BooleanVar(value=False)
            ctk.CTkCheckBox(frame, text="Include PC Emissivity (256 bands)",
                            variable=self.pc_emis_var, font=ctk.CTkFont(size=12)
                            ).grid(row=start_row+1, column=2, columnspan=2, padx=10, pady=5, sticky='w')

    # ------------------------------------------------------------------
    # VIIRS band selection handler
    # ------------------------------------------------------------------

    def _on_viirs_bands_change(self, selection):
        """Show/hide custom band entry based on VIIRS band selection."""
        if selection == 'Custom':
            self.viirs_custom_label.grid()
            self.viirs_custom_entry.grid()
        else:
            self.viirs_custom_label.grid_remove()
            self.viirs_custom_entry.grid_remove()

    def _get_viirs_selected_bands(self):
        """
        Resolve the VIIRS band selection dropdown into a list of band names,
        or None for 'All bands'.
        """
        if not hasattr(self, 'viirs_bands_var'):
            return None

        mode = self.viirs_bands_var.get()

        band_presets = {
            'All bands':               None,
            'TIR only (M14-M16)':      ['M14', 'M15', 'M16'],
            'Split-window (M15, M16)': ['M15', 'M16'],
            'MIR + TIR (M12-M16)':     ['M12', 'M13', 'M14', 'M15', 'M16'],
            'VSWIR only (M01-M11)':    [f'M{i:02d}' for i in range(1, 12)],
            'I-band TIR (I04, I05)':   ['I04', 'I05'],
            'VNP21 core (LST + Emissivity)': ['LST', 'Emis_14', 'Emis_15', 'Emis_16'],
        }

        if mode in band_presets:
            return band_presets[mode]
        elif mode == 'Custom':
            text = self.viirs_custom_entry.get().strip()
            if text:
                return [b.strip() for b in text.split(',') if b.strip()]
        return None

    # ------------------------------------------------------------------
    # I/O helpers
    # ------------------------------------------------------------------
 
    def _browse_input(self):
        path = filedialog.askdirectory(title="Select Input Folder")
        if path:
            self.input_entry.delete(0, 'end')
            self.input_entry.insert(0, path)
 
    def _browse_output(self):
        path = filedialog.askdirectory(title="Select Output Folder")
        if path:
            self.output_entry.delete(0, 'end')
            self.output_entry.insert(0, path)
 
    def _clear_log(self):
        self.log_text.configure(state='normal')
        self.log_text.delete('1.0', 'end')
 
    def _append_log(self, timestamp, level, message):
        self.log_text.configure(state='normal')
        self.log_text.insert('end', f"[{timestamp}] {message}\n")
        self.log_text.see('end')
 
    def _poll_log(self):
        for timestamp, level, message in self.log_queue.get_all():
            self._append_log(timestamp, level, message)
        self.after(100, self._poll_log)
 
    # ------------------------------------------------------------------
    # Processing
    # ------------------------------------------------------------------
 
    def _start_processing(self):
        if self.current_sensor is None:
            messagebox.showwarning("No Sensor Selected",
                                   "Please select a sensor from the left panel.")
            return
 
        input_path = self.input_entry.get().strip()
        if not input_path or not os.path.exists(input_path):
            messagebox.showwarning("Invalid Input",
                                   "Please select a valid input folder.")
            return
 
        output_path = self.output_entry.get().strip() or None
 
        fmt           = self.format_menu.get()
        write_geotiff = fmt in ('GeoTIFF', 'Both')
        write_envi    = fmt in ('ENVI', 'Both')
 
        self.processing = True
        self._cancel_flag = False
        self.process_btn.configure(text="⏳  Processing...", state='disabled')
        self.cancel_btn.configure(state='normal')
        self.progress.set(0)
        self.progress.start()
 
        kwargs = {
            'input_dir':    input_path,
            'output_dir':   output_path,
            'write_geotiff': write_geotiff,
            'write_envi':    write_envi,
            'sensor_id':     self._active_sensor_id,
            'overwrite':     self.overwrite_var.get(),
        }
 
        if hasattr(self, 'aster_product_var'):
            kwargs['aster_product'] = self.aster_product_var.get()
        if hasattr(self, 'aster_version_var'):
            kwargs['aster_version'] = self.aster_version_var.get()
        if hasattr(self, 'av3_water_var'):
            kwargs['exclude_water'] = self.av3_water_var.get()
        if hasattr(self, 'projection_var'):
            kwargs['projection'] = self.projection_var.get().lower()
        if hasattr(self, 'reflectance_var'):
            kwargs['vswir_reflectance'] = self.reflectance_var.get()
        if hasattr(self, 'resolution_entry'):
            res_text = self.resolution_entry.get().strip()
            kwargs['resolution_m'] = (float(res_text)
                                      if res_text and res_text != 'auto' else None)
        if hasattr(self, 'emit_bands_var'):
            kwargs['emit_bands'] = self.emit_bands_var.get()
            kwargs['emit_rmse']  = self.emit_rmse_var.get()
        if hasattr(self, 'pc_emis_var'):
            kwargs['include_pc'] = self.pc_emis_var.get()

        # VIIRS band selection
        viirs_bands = self._get_viirs_selected_bands()
        if viirs_bands is not None:
            kwargs['viirs_selected_bands'] = viirs_bands

        # --- Bbox from shared panel ---
        bbox_result = self._get_bbox_tuple()
        if bbox_result is not None:
            if isinstance(bbox_result, tuple) and bbox_result[0] == 'file':
                kwargs['subset_file'] = bbox_result[1]
            else:
                kwargs['bbox'] = bbox_result

        # --- Spectral opts from shared panel ---
        spectral = self._get_spectral_opts()
        if spectral is not None:
            kwargs['spectral_opts'] = spectral
 
        threading.Thread(target=self._run_processing,
                         args=(kwargs,), daemon=True).start()
 
    def _run_processing(self, kwargs):
        sensor_id = kwargs.pop('sensor_id')
        bbox = kwargs.pop('bbox', None)
        subset_file = kwargs.pop('subset_file', None)
        spectral_opts = kwargs.pop('spectral_opts', None)
        overwrite = kwargs.pop('overwrite', True)
 
        try:
            # --- Resolve subset file to bbox tuple for non-ECOSTRESS sensors ---
            if subset_file and bbox is None and sensor_id != 'ecostress_l2':
                try:
                    from subset_tool import bbox_from_file
                    resolved = bbox_from_file(subset_file)
                    bbox = (resolved.min_lat, resolved.max_lat,
                            resolved.min_lon, resolved.max_lon)
                    self.log_queue.log(
                        f"Resolved subset file to bbox: "
                        f"lat [{bbox[0]:.4f}, {bbox[1]:.4f}], "
                        f"lon [{bbox[2]:.4f}, {bbox[3]:.4f}]")
                except ImportError:
                    self.log_queue.log(
                        "WARNING: subset_tool.py not found. "
                        "Cannot resolve file to bounding box.", 'warning')
                except Exception as e:
                    self.log_queue.log(
                        f"WARNING: Could not read subset file: {e}", 'warning')

            self.log_queue.log(
                f"Starting {self.current_sensor['name']} processing...")
            self.log_queue.log(f"Input: {kwargs['input_dir']}")
            if bbox:
                self.log_queue.log(
                    f"Bbox: lat [{bbox[0]}, {bbox[1]}], lon [{bbox[2]}, {bbox[3]}]")
            if subset_file and sensor_id == 'ecostress_l2':
                self.log_queue.log(f"Subset file: {subset_file}")
            if spectral_opts:
                self.log_queue.log(f"Spectral: {spectral_opts}")

            # Check cancel flag
            if self._cancel_flag:
                self.log_queue.log("Processing cancelled by user.")
                return
 
            # ------------------------------------------------------------
            # Converter dispatch
            # ------------------------------------------------------------
 
            if sensor_id == 'aster_v004':
                from ASTER_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    aster_product=kwargs.get('aster_product', 'AST_05'),
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )

            elif sensor_id == 'aster_v003':
                from ASTER_Converter import process_directory_v003
                process_directory_v003(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    aster_product=kwargs.get('aster_product', 'AST_05'),
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'landsat_l1':
                from Landsat_L1_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'landsat_l2':
                from Landsat_L2_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'sentinel2':
                from Sentinel2_L2A_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id in ('s3_slstr_l1', 's3_slstr_l2',
                               's3_olci_l1', 's3_syn_l2'):
                from Sentinel3_Converter import process_directory
                # Map sensor_id to product filter
                s3_filter = {
                    's3_slstr_l1': 'SL_1_RBT',
                    's3_slstr_l2': 'SL_2_LST',
                    's3_olci_l1':  'OL_1_EFR',
                    's3_syn_l2':   'SY_2_SYN',
                }.get(sensor_id)
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite,
                    product_filter=s3_filter
                )

            elif sensor_id == 'ecostress_l1ct':
                from ECOSTRESS_L1TC_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    output_geotiff=kwargs['write_geotiff'],
                    output_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'ecostress_l2':
                from ECOSTRESS_L2_Converter import process_directory
                # Build subset_opts for ECOSTRESS (which uses its own format)
                eco_subset = None
                if subset_file:
                    eco_subset = {'mode': 'subset', 'bbox': None,
                                  'tile_size': None, 'subset_file': subset_file}
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    projection=kwargs.get('projection', 'geographic'),
                    subset_opts=eco_subset,
                    bbox=bbox,
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'viirs':
                from VIIRS_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite,
                    selected_bands=kwargs.get('viirs_selected_bands')
                )
 
            elif sensor_id in ('modis_refl', 'modis_lst'):
                from MODIS_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )

            elif sensor_id == 'emit':
                from EMIT_L2A_Converter import process_directory
                spectral_emit = spectral_opts
                good_only     = True
 
                band_mode = kwargs.get('emit_bands', 'Good wavelengths')
                if band_mode == 'All bands':
                    good_only = False
                elif band_mode in ('Match Landsat', 'Match Sentinel-2',
                                   'Match ASTER', 'Match MODIS'):
                    sensor_map = {
                        'Match Landsat':    'landsat89',
                        'Match Sentinel-2': 'sentinel2',
                        'Match ASTER':      'aster',
                        'Match MODIS':      'modis_land',
                    }
                    spectral_emit = {
                        'mode':          'sensor',
                        'sensor':        sensor_map[band_mode],
                        'nearest':       False,
                        'rsr_dir':       None,
                        'exclude_water': True,
                    }
                elif band_mode == 'VNIR only':
                    spectral_emit = {'mode': 'region', 'region': 'vnir',
                                     'exclude_water': True}
                elif band_mode == 'SWIR only':
                    spectral_emit = {'mode': 'region', 'region': 'swir',
                                     'exclude_water': True}
 
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    good_only=good_only,
                    include_rmse=kwargs.get('emit_rmse', True),
                    spectral_opts=spectral_emit,
                    bbox=bbox,
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'airs':
                from AIRS_L1B_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    bbox=bbox,
                    spectral_opts=spectral_opts,
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'master_l1b':
                from MASTER_L1B_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    vswir_as_reflectance=kwargs.get('vswir_reflectance', False),
                    resolution_m=kwargs.get('resolution_m') or 50.0,
                    bbox=bbox,
                    spectral_opts=spectral_opts,
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'master_l2':
                from MASTER_L2_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    resolution_m=kwargs.get('resolution_m') or 50.0,
                    bbox=bbox,
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'hytes_l1':
                from HyTES_L1_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    resolution_m=kwargs.get('resolution_m'),
                    bbox=bbox,
                    spectral_opts=spectral_opts,
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'hytes_l2':
                from HyTES_L2_Converter import process_directory
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    resolution_m=kwargs.get('resolution_m'),
                    include_pc=kwargs.get('include_pc', False),
                    bbox=bbox,
                    spectral_opts=spectral_opts,
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )
 
            elif sensor_id == 'aviris35':
                self.log_queue.log("Importing AVIRIS-3/5 converter...")
                from AVIRIS3_L1B_Converter import process_directory
                self.log_queue.log("Starting AVIRIS-3/5 processing...")
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    exclude_water=kwargs.get('exclude_water', True),
                    spectral_opts=spectral_opts,
                    bbox=bbox,
                    resolution_m=kwargs.get('resolution_m'),
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )

            elif sensor_id == 'aviris35_l2a':
                self.log_queue.log("Importing AVIRIS-3/5 L2A converter...")
                from AVIRIS3_L2A_Converter import process_directory
                self.log_queue.log("Starting AVIRIS-3/5 L2A processing...")
                process_directory(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    exclude_water=kwargs.get('exclude_water', True),
                    spectral_opts=spectral_opts,
                    bbox=bbox,
                    write_geotiff=kwargs['write_geotiff'],
                    write_envi=kwargs['write_envi'],
                    overwrite=overwrite
                )

            elif sensor_id in ('aviris_ng_l1b', 'aviris_ng_l2',
                               'aviris_classic_l1b', 'aviris_classic_l2'):
                self.log_queue.log("Importing AVIRIS-NG/Classic converter...")
                from AVIRIS_Converter import process_file, discover_files, \
                    parse_envi_header, parse_float_list, get_band_mask, \
                    parse_map_info, build_transform, build_crs
                from AVIRIS_Converter import process_directory as aviris_process_dir
                self.log_queue.log("Starting AVIRIS-NG/Classic processing...")
                aviris_process_dir(
                    kwargs['input_dir'], kwargs.get('output_dir'),
                    geotiff=kwargs['write_geotiff'],
                    envi=kwargs['write_envi'],
                    exclude_water=kwargs.get('exclude_water', True),
                    use_bbl=getattr(self, 'av_bbl_var', ctk.BooleanVar(value=True)).get()
                        if hasattr(self, 'av_bbl_var') else True,
                    overwrite=overwrite
                )
 
            else:
                self.log_queue.log(
                    f"ERROR: No converter found for sensor '{sensor_id}'", 'error')
                return
 
            self.log_queue.log("Processing complete!", 'success')

            # Note: output file organization into GeoTIFF/ENVI subfolders
            # is handled by each converter internally.

        except Exception as e:
            self.log_queue.log(f"ERROR: {e}", 'error')
            import traceback
            self.log_queue.log(traceback.format_exc(), 'error')
 
        finally:
            self.after(0, self._processing_complete)
 
    def _processing_complete(self):
        self.processing = False
        self._cancel_flag = False
        self.process_btn.configure(text="▶  Process Data", state='normal')
        self.cancel_btn.configure(state='disabled', text="■  Cancel")
        self.open_folder_btn.configure(state='normal')
        self.progress.stop()
        self.progress.set(1)

    def _cancel_processing(self):
        """Set cancel flag — processing thread checks this and stops."""
        if self.processing:
            self._cancel_flag = True
            self.log_queue.log("Cancelling... (will stop after current file completes)", 'warning')
            self.cancel_btn.configure(state='disabled', text="■  Cancelling...")

    def _open_output_folder(self):
        """Open the output folder in the system file explorer."""
        import subprocess

        # Determine output path: use explicit output, or default to input/converted
        out_path = self.output_entry.get().strip()
        if not out_path:
            in_path = self.input_entry.get().strip()
            if in_path:
                out_path = os.path.join(in_path, 'converted')

        if out_path and os.path.isdir(out_path):
            # Windows: open in Explorer
            subprocess.Popen(['explorer', os.path.normpath(out_path)])
        elif out_path:
            # Try parent directory if exact path doesn't exist yet
            parent = os.path.dirname(out_path)
            if os.path.isdir(parent):
                subprocess.Popen(['explorer', os.path.normpath(parent)])
            else:
                from tkinter import messagebox
                messagebox.showinfo("Folder Not Found",
                    f"Output folder does not exist yet:\n{out_path}")
        else:
            from tkinter import messagebox
            messagebox.showinfo("No Output Path",
                "Please run processing first or set an output folder.")

    @staticmethod
    def _organize_output_files(out_dir):
        """
        Move GeoTIFF (.tif) and ENVI (.dat/.hdr) files into GeoTIFF/ and ENVI/
        subfolders. For AST_09T products, further separates SurfaceRadiance
        and SkyIrradiance into their own subfolders.
        Skips files already in subfolders.
        """
        import shutil

        geotiff_dir = os.path.join(out_dir, 'GeoTIFF')
        envi_dir    = os.path.join(out_dir, 'ENVI')

        moved = 0
        for fname in os.listdir(out_dir):
            fpath = os.path.join(out_dir, fname)
            if not os.path.isfile(fpath):
                continue

            ext = os.path.splitext(fname)[1].lower()
            fname_lower = fname.lower()

            if ext == '.tif':
                # AST_09T subfolder separation
                if 'surfaceradiance' in fname_lower:
                    dest = os.path.join(geotiff_dir, 'SurfaceRadiance')
                elif 'skyirradiance' in fname_lower:
                    dest = os.path.join(geotiff_dir, 'SkyIrradiance')
                else:
                    dest = geotiff_dir
                os.makedirs(dest, exist_ok=True)
                shutil.move(fpath, os.path.join(dest, fname))
                moved += 1
            elif ext in ('.dat', '.hdr'):
                # AST_09T subfolder separation
                if 'surfaceradiance' in fname_lower:
                    dest = os.path.join(envi_dir, 'SurfaceRadiance')
                elif 'skyirradiance' in fname_lower:
                    dest = os.path.join(envi_dir, 'SkyIrradiance')
                else:
                    dest = envi_dir
                os.makedirs(dest, exist_ok=True)
                shutil.move(fpath, os.path.join(dest, fname))
                moved += 1

        if moved > 0:
            print(f"  Organized {moved} file(s) into GeoTIFF/ and ENVI/ subfolders")
 
 
# ============================================================================
# Main Application
# ============================================================================
 
class RemoteSensingToolkit(ctk.CTk):
    def __init__(self):
        super().__init__()
 
        self.title(f"{APP_NAME} v{APP_VERSION}")
        self.minsize(*APP_MIN_SIZE)

        # Maximize window on startup (Windows)
        try:
            self.state('zoomed')
        except Exception:
            # Fallback for non-Windows: set geometry to screen size
            self.geometry(f"{self.winfo_screenwidth()}x{self.winfo_screenheight()}+0+0")
 
        ctk.set_appearance_mode('dark')
        ctk.set_default_color_theme('blue')
 
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
 
        self.log_queue = LogQueue()
        self._setup_stdout_redirect()
        self._build_sidebar()
        self._build_main_area()
 
        self.sensor_cards      = []
        self.current_sensor_id = None
 
        self._show_sensor_selection()
 
    def _setup_stdout_redirect(self):
        log_q = self.log_queue
 
        class _Redirector:
            def __init__(self):
                self.original = sys.stdout
            def write(self, text):
                t = text.strip()
                if t:
                    log_q.log(t)
                if self.original is not None:
                    self.original.write(text)
            def flush(self):
                if self.original is not None:
                    self.original.flush()
 
        sys.stdout = _Redirector()
        sys.stderr = _Redirector()  # Also redirect stderr to the log
 
    def _build_sidebar(self):
        self.sidebar = ctk.CTkFrame(self, width=260, corner_radius=0,
                                     fg_color=COLORS['sidebar_bg'])
        self.sidebar.grid(row=0, column=0, sticky='nsew')
        self.sidebar.grid_propagate(False)
        self.sidebar.grid_columnconfigure(0, weight=1)
 
        title_frame = ctk.CTkFrame(self.sidebar, fg_color='transparent')
        title_frame.grid(row=0, column=0, sticky='ew', padx=15, pady=(20, 10))
 
        ctk.CTkLabel(title_frame, text="🛰️ RSDTK",
                     font=ctk.CTkFont(size=24, weight='bold'),
                     text_color=COLORS['accent']).pack(anchor='w')
        ctk.CTkLabel(title_frame, text="Remote Sensing Data Toolkit",
                     font=ctk.CTkFont(size=11),
                     text_color='#b0b0b0').pack(anchor='w')
        ctk.CTkLabel(title_frame, text=f"v{APP_VERSION}",
                     font=ctk.CTkFont(size=10),
                     text_color='#b0b0b0').pack(anchor='w', pady=(0, 5))
 
        ctk.CTkFrame(self.sidebar, height=1,
                     fg_color=COLORS['sidebar_hover']
                     ).grid(row=1, column=0, sticky='ew', padx=15, pady=5)
 
        nav_frame = ctk.CTkFrame(self.sidebar, fg_color='transparent')
        nav_frame.grid(row=2, column=0, sticky='ew', padx=10)
        nav_frame.grid_columnconfigure(0, weight=1)
 
        self.nav_buttons = {}
 
        btn_select = SidebarButton(nav_frame, text="Sensor Selection", icon="📋",
                                    command=self._show_sensor_selection)
        btn_select.grid(row=0, column=0, sticky='ew', pady=2)
        self.nav_buttons['select'] = btn_select
 
        btn_process = SidebarButton(nav_frame, text="Processing", icon="⚙️",
                                     command=self._show_processing)
        btn_process.grid(row=1, column=0, sticky='ew', pady=2)
        self.nav_buttons['process'] = btn_process
 
        ctk.CTkFrame(self.sidebar, height=1,
                     fg_color=COLORS['sidebar_hover']
                     ).grid(row=3, column=0, sticky='ew', padx=15, pady=10)
 
        bottom = ctk.CTkFrame(self.sidebar, fg_color='transparent')
        bottom.grid(row=4, column=0, sticky='sew', padx=15, pady=15)
        bottom.grid_columnconfigure(0, weight=1)
        self.sidebar.grid_rowconfigure(4, weight=1)
 
        ctk.CTkLabel(bottom, text="Appearance:",
                     font=ctk.CTkFont(size=12),
                     text_color='#b0b0b0').pack(anchor='w', pady=(0, 5))
        self.appearance_menu = ctk.CTkOptionMenu(
            bottom, values=['Dark', 'Light', 'System'],
            command=self._change_appearance, width=140,
            font=ctk.CTkFont(size=12))
        self.appearance_menu.pack(anchor='w')
        self.appearance_menu.set('Dark')
 
    def _build_main_area(self):
        self.main_frame = ctk.CTkFrame(self, fg_color='transparent')
        self.main_frame.grid(row=0, column=1, sticky='nsew')
        self.main_frame.grid_columnconfigure(0, weight=1)
        self.main_frame.grid_rowconfigure(0, weight=1)
 
        self.sensor_scroll = ctk.CTkScrollableFrame(self.main_frame,
                                                      fg_color='transparent')
        self.sensor_scroll.grid_columnconfigure(0, weight=1)
 
        self.processing_panel = ProcessingPanel(self.main_frame, self.log_queue)
 
    def _show_sensor_selection(self):
        self.processing_panel.grid_forget()
        self.sensor_scroll.grid(row=0, column=0, sticky='nsew', padx=10, pady=10)
        self._set_active_nav('select')
        self._populate_sensors()
 
    def _show_processing(self):
        self.sensor_scroll.grid_forget()
        self.processing_panel.grid(row=0, column=0, sticky='nsew')
        self._set_active_nav('process')
 
    def _set_active_nav(self, key):
        for k, btn in self.nav_buttons.items():
            btn.set_active(k == key)
 
    def _populate_sensors(self):
        for widget in self.sensor_scroll.winfo_children():
            widget.destroy()
        self.sensor_cards = []
 
        for group_name, sensors in SENSOR_GROUPS.items():
            ctk.CTkLabel(self.sensor_scroll, text=group_name,
                         font=ctk.CTkFont(size=16, weight='bold'),
                         anchor='w').grid(sticky='w', padx=15, pady=(15, 8))
            for sensor in sensors:
                card = SensorCard(self.sensor_scroll, sensor,
                                   on_select=self._on_sensor_selected)
                card.grid(sticky='ew', padx=15, pady=3)
                self.sensor_cards.append(card)
 
    def _on_sensor_selected(self, sensor_info):
        for card in self.sensor_cards:
            card.set_selected(card.sensor_info['id'] == sensor_info['id'])
        self.current_sensor_id = sensor_info['id']
        self.processing_panel.set_sensor(sensor_info)
        self.after(300, self._show_processing)
 
    def _change_appearance(self, mode):
        ctk.set_appearance_mode(mode.lower())
 
 
# ============================================================================
# Entry point
# ============================================================================
 
if __name__ == '__main__':
    app = RemoteSensingToolkit()
    app.mainloop()