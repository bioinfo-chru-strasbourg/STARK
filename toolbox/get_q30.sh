#############################################
###
### 	> get Q30 <
###
#############################################

##########################
### USAGE

#./get_q30.sh <PATH to metrics.fastqc.txt>

##########################
### SCRIPT



FILE=$1
ALL=$(cat $FILE | awk '/>>Per sequence quality scores/,/>>END_MODULE/' | head -n -1 | tail -n+3 | awk '{s+=$2}END{print s}');
Q30=$(cat $FILE | awk '/>>Per sequence quality scores/,/>>END_MODULE/' | awk '/>>Per sequence quality scores/,/30\t/' | head -n -1 | tail -n+3 | awk '{s+=$2}END{print s}');
RES=$(bc <<< "scale = 4; (($ALL - $Q30) / $ALL * 100)");
#echo ALL=$ALL
#echo Q30=$Q30
#echo RES=$RES
echo Q30: $RES

