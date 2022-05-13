@REM much less dangerous version of git push, with (hopefully) no merge conflicts and well-organized, informative commit history.
@REM first checks that the working directory is up-to-date with remote, and has no uncommitted changes.
@REM Then soft reset, squashing commits, with backup branch.  
@REM squashed commit has same commit msg as last commit on build.
@REM limitation: not meant for multiple devs.  i.e. this works for my use-case - solo dev - and not much else.
	
@cd %~dp0
@set "ERRORLEVEL=0"

@for /f "delims=" %%a in ('git.exe rev-parse --abbrev-ref HEAD') do set "branchname=%%a"
@if %ERRORLEVEL% neq 0 goto :ERROR
@git.exe fetch
@if %ERRORLEVEL% neq 0 goto :ERROR
@for /f %%i in ('git.exe log HEAD..origin/master --oneline') do @if /i "%%~a" NEQ "" set "ERRORLEVEL=1"
@if %ERRORLEVEL% neq 0 goto :BEHIND
@for /f %%i in ('git.exe status --porcelain') do @if /i "%%~a" NEQ "" set "ERRORLEVEL=1"
@if %ERRORLEVEL% neq 0 goto :NOTCLEAN
git.exe checkout -b build-config
@if %ERRORLEVEL% neq 0 goto :ERROR
git.exe branch -u origin/master
@if %ERRORLEVEL% neq 0 goto :ERROR
git.exe reset --soft origin/master
@if %ERRORLEVEL% neq 0 goto :ERROR
git commit --reuse-message=%branchname%@{0}
@if %ERRORLEVEL% neq 0 goto :ERROR
git.exe push origin HEAD:master
@if %ERRORLEVEL% neq 0 goto :ERROR
git.exe branch -D lm
git.exe branch -m lm
@if %ERRORLEVEL% neq 0 goto :ERROR
git.exe checkout %branchname%
@if %ERRORLEVEL% neq 0 goto :ERROR
git.exe reset origin/master
@if %ERRORLEVEL% neq 0 goto :ERROR
git.exe fetch
@goto :EOF

:BEHIND
@echo currently not up-to-date with remote repo, script aborted.
@set "ERRORLEVEL=0"
@exit /b 1
:NOTCLEAN
@echo working directory not clean, script aborted.
@set "ERRORLEVEL=0"
@exit /b 1
:ERROR
@echo Failed with error #%errorlevel%.
@exit /b %errorlevel%