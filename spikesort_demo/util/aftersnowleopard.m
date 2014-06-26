function bool = aftersnowleopard()
if ~ismac(), bool = false; return; end

[~,systemrelease] = system('uname -r');
systemrelease = sscanf(systemrelease, '%d.%d.%d.');
bool = systemrelease(1) > 10;