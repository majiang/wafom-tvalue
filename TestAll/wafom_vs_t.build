set PATH=C:\Program Files\D\dmd2\windows\bin;C:\Program Files\Microsoft SDKs\Windows\v6.0A\\bin;%PATH%

echo asianoption.d >TestAll\wafom_vs_t.build.rsp
echo digitalnet.d >>TestAll\wafom_vs_t.build.rsp
echo finitefield.d >>TestAll\wafom_vs_t.build.rsp
echo graycode.d >>TestAll\wafom_vs_t.build.rsp
echo integral.d >>TestAll\wafom_vs_t.build.rsp
echo main.d >>TestAll\wafom_vs_t.build.rsp
echo pointset.d >>TestAll\wafom_vs_t.build.rsp
echo sobol.d >>TestAll\wafom_vs_t.build.rsp
echo testfunction.d >>TestAll\wafom_vs_t.build.rsp
echo tvalue.d >>TestAll\wafom_vs_t.build.rsp
echo wafom.d >>TestAll\wafom_vs_t.build.rsp

dmd -g -debug -X -Xf"TestAll\wafom_vs_t.json" -deps="TestAll\wafom_vs_t.dep" -of"TestAll\wafom_vs_t.exe" -map "TestAll\wafom_vs_t.map" -L/NOMAP -unittest -debug=verbose -debug=working @TestAll\wafom_vs_t.build.rsp
if errorlevel 1 goto reportError
if not exist "TestAll\wafom_vs_t.exe" (echo "TestAll\wafom_vs_t.exe" not created! && goto reportError)

goto noError

:reportError
echo Building TestAll\wafom_vs_t.exe failed!

:noError
