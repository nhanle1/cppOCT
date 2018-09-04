@echo off
set /A param1=0
echo 	Processing mode: %param1%

set /A temp=%param1% %% 8
if %temp% GTR 3 (echo 		File compression ON) else (
	echo 		File compression OFF)
set /A temp1=%temp% %% 4
set /A temp=%temp1% %% 2
if %temp% EQU 1 (echo 		Interpolation ON
	if %temp1% GTR 1 (echo 		- Cubic spline) else (
		echo 		- Linear spline)
	) else (
		echo 		Interpolation OFF)

echo.
echo 	Batch name found: 
for %%a in (*.oct) do (
	set octbatch=!octbatch! %%~nxa
	echo 		%%~nxa)
echo.
echo		Press any keys to start OCT processing...
pause>nul
for %%a in (*.oct) do (start UWgingiva.exe 00 %%~nxa)

:: for interpolation change parameter 1 
:: (after UWgingiva.exe) first digit to 1 (ON) and 0 (OFF)
:: second digit controls the output binary compression
:: eg: 10 interpolation ON, compression OFF



:: Previous .bat version manual naming
:: start  UWgingiva.exe ^
:: 0^
:: 2017.10.20_14.41.31_s03F_G1_8x8mm_r8.oct^
:: 2017.10.20_14.41.31_s03F_G1_8x8mm_r8-Copy.oct


