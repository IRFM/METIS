function watchdog_mdsplus_callback(myTimerObj, thisEvent)

global watchdog_mdsplus

watchdog_mdsplus = myTimerObj;
if ~isempty(which('mdsipmex'))
    clear mdsipmex
else
    mdsdisconnect;
end
warning off
delete(watchdog_mdsplus);
warning on
disp('MDS connnection closed');
watchdog_mdsplus=[];