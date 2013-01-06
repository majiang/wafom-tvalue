set PATH=C:\Program Files\D\dmd2\windows\bin;C:\Program Files\Microsoft SDKs\Windows\v6.0A\\bin;%PATH%
dmd -g -debug -X -Xf"Debug\wafom_vs_t.json" -deps="Debug\wafom_vs_t.dep" -of"Debug\wafom_vs_t.exe_cv" -map "Debug\wafom_vs_t.map" -L/NOMAP -unittest digitalnet.d finitefield.d main.d sobol.d tvalue.d wafom.d
if errorlevel 1 goto reportError
if not exist "Debug\wafom_vs_t.exe_cv" (echo "Debug\wafom_vs_t.exe_cv" not created! && goto reportError)
echo Converting debug information...
"C:\Program Files\VisualD\cv2pdb\cv2pdb.exe" -D2.043 "Debug\wafom_vs_t.exe_cv" "Debug\wafom_vs_t.exe"
if errorlevel 1 goto reportError
if not exist "Debug\wafom_vs_t.exe" (echo "Debug\wafom_vs_t.exe" not created! && goto reportError)

goto noError

:reportError
echo Building Debug\wafom_vs_t.exe failed!

:noError
