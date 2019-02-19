
ANALYSISTIME=$(./analysis.time.sh); while [ "$(echo "$ANALYSISTIME" | grep RUNNING -c)" == "1" ]; do clear; echo "$ANALYSISTIME"; sleep 10; ANALYSISTIME=$(./analysis.time.sh); done; clear; echo "$ANALYSISTIME";
