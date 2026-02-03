!include "MUI2.nsh"

!ifndef APP_NAME
  !define APP_NAME "MID Nano"
!endif
!ifndef APP_ID
  !define APP_ID "mid-nano"
!endif
!ifndef APP_VERSION
  !define APP_VERSION "0.0.0"
!endif
!ifndef SOURCE_DIR
  !define SOURCE_DIR "stage"
!endif
!ifndef OUT_FILE
  !define OUT_FILE "${APP_NAME}-${APP_VERSION}-win-x64.exe"
!endif

Name "${APP_NAME}"
OutFile "${OUT_FILE}"

InstallDir "$LOCALAPPDATA\\${APP_NAME}"
RequestExecutionLevel user

!define MUI_ABORTWARNING

!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES

!insertmacro MUI_LANGUAGE "English"

Section "Install"
  SetOutPath "$INSTDIR"
  File /r "${SOURCE_DIR}\\*"

  WriteUninstaller "$INSTDIR\\Uninstall.exe"

  CreateDirectory "$SMPROGRAMS\\${APP_NAME}"
  CreateShortcut "$SMPROGRAMS\\${APP_NAME}\\${APP_NAME}.lnk" "$INSTDIR\\${APP_NAME}.exe" "" "$INSTDIR\\${APP_NAME}.exe" 0
  CreateShortcut "$DESKTOP\\${APP_NAME}.lnk" "$INSTDIR\\${APP_NAME}.exe" "" "$INSTDIR\\${APP_NAME}.exe" 0
  CreateShortcut "$SMPROGRAMS\\${APP_NAME}\\Uninstall.lnk" "$INSTDIR\\Uninstall.exe"
SectionEnd

Section "Uninstall"
  Delete "$SMPROGRAMS\\${APP_NAME}\\${APP_NAME}.lnk"
  Delete "$DESKTOP\\${APP_NAME}.lnk"
  Delete "$SMPROGRAMS\\${APP_NAME}\\Uninstall.lnk"
  RMDir "$SMPROGRAMS\\${APP_NAME}"

  RMDir /r "$INSTDIR"
SectionEnd
