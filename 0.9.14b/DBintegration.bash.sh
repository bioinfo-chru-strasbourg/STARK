

# var
res=/media/IRCV2/V2/RES/ALL
config=config.db2.ini
log=DBintegration.bash.db2.log
err=DBintegration.bash.db2.err

# init
echo "" > $log
echo "" > $err

# run
for f in `ls $res`
do
	time ./DBintegration.sh $config $f 1>>$log 2>>$err
done
