@REM git commit without the annoying parts.
@REM never commits outside of the correct branch (build).
@REM doesn't nag me for a commit message.
@REM git --allow-empty-message doesn't work in gitextensions, but this still allows commentless commits.

@cd %~dp0
@set "ERRORLEVEL=0"
git.exe checkout build
@if %ERRORLEVEL% neq 0 goto :EOF
if ["%~1"]==[""] (
git.exe commit -m "." --allow-empty) else (
git.exe commit -m "%*" --allow-empty)
@if %ERRORLEVEL% neq 0 goto :EOF
git.exe --no-pager log -n 5 --pretty=oneline
:EOF