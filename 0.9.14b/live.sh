#clear; 
echo "Launch NGS LIVE";

NGS_SCRIPT=`dirname $0`
echo $NGS_SCRIPT

TERMINAL="xterm -hold" # gnome-terminal xterm -hold

$TERMINAL -e $NGS_SCRIPT/monitor.live.sh --noclose 2>/dev/null &
$TERMINAL -e $NGS_SCRIPT/analysis.time.live.sh --noclose 2>/dev/null &
$TERMINAL -e $NGS_SCRIPT/tail.live.sh --noclose 2>/dev/null &

exit 0;

