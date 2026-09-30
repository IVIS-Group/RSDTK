; ============================================================================
; Remote Sensing Data Toolkit (RSDTK) — Inno Setup Installer Script
; ============================================================================
; Save as: RSDTK_Installer.iss
; Build:   Open in Inno Setup Compiler → Build → Compile (Ctrl+F9)
; Output:  Installer_Output\RSDTK_v1.2_Setup.exe
; ============================================================================

[Setup]
AppName=Remote Sensing Data Toolkit
AppVersion=1.2
AppPublisher=University of Pittsburgh — Geology & Environmental Science
DefaultDirName={autopf}\RSDTK
DefaultGroupName=RSDTK
UninstallDisplayIcon={app}\RSDTK.exe
OutputDir=Installer_Output
OutputBaseFilename=RSDTK_v1.2_Setup
Compression=lzma2
SolidCompression=yes
SetupIconFile=icon.ico
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
WizardStyle=modern
DisableProgramGroupPage=yes
PrivilegesRequired=lowest

[Files]
; The compiled .exe from PyInstaller
Source: "dist\RSDTK.exe"; DestDir: "{app}"; Flags: ignoreversion
; Include the icon
Source: "icon.ico"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
; Start Menu shortcut
Name: "{userprograms}\RSDTK"; Filename: "{app}\RSDTK.exe"; IconFilename: "{app}\icon.ico"
Name: "{userprograms}\Uninstall RSDTK"; Filename: "{uninstallexe}"
; Desktop shortcut (optional, user chooses during install)
Name: "{userdesktop}\RSDTK"; Filename: "{app}\RSDTK.exe"; IconFilename: "{app}\icon.ico"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"

[Run]
; Option to launch after install
Filename: "{app}\RSDTK.exe"; Description: "Launch RSDTK"; Flags: nowait postinstall skipifsilent
