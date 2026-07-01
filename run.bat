@echo off

echo Compiling Java files...
javac ExpressionEvaluator.java EnergyPredictSBGE.java GEIndividual.java BNFGrammar.java

echo.
echo Select mode:
echo 1 - Structural-Based GE (SBGE)

set /p choice=Enter choice (1): 

if "%choice%"=="1" goto sbge
goto invalid

:sbge
echo Running Structural-Based GE (EnergyPredictSBGE)...
java EnergyPredictSBGE
goto end

:invalid
echo Invalid choice.

:end
pause