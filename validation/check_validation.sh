
TMP=/tmp/$RANDOM

DIR=$1
PATTERN=$2

if [ "$DIR" == "" ] ; then DIR="."; fi

ERROR=0

for d in $DIR; do
	#for t in $d/VALIDATION_*; do
	for t in $d/$PATTERN*; do
		if [ -d "$t" ]; then
			MSG=" OK "
			if (($(cat $t/ana*log $t/run*log 2>/dev/null | grep "\*\*\*" -c))); then
				MSG=" ERROR "
				((ERROR++))
			fi;
			echo "CHECK '"$(basename $t)"': $MSG" >> $TMP;
		fi;
	done;
done;

if [ -s $TMP ] && (($(grep ^ -c $TMP))); then
	cat $TMP 2>/dev/null | column -t
	if (($ERROR)); then
		echo -e "VALIDATION ERROR"
	else
		echo -e "VALIDATION OK"
	fi;
	rm -f $TMP
fi;



