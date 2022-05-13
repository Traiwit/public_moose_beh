@REM git status but more informative.
@cd %~dp0
@git.exe branch
@git.exe --no-pager log -n 10 --pretty=oneline
@git.exe status
