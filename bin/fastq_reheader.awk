# $1 is read
# $2 is index1
# $3 is index2 (if exist)
# READ is the read number (usually "1" or "2")
BEGIN {
		H=""
		R=""
		I1=""
		I2=""
		if (READ=="") {
			READ="1"
		}
}
{ 
	if(match($1, /^@/) )
	{ 
		H=$1
	} else if( H != "") {
		R=$1
		I1=$2
		I2=$3
		I=""
		if (I1!="") {
			I=I I1
		}
		if (I!="") {SEP="+"}
		if (I2!="") {
			I=I SEP I2
		}
		gsub(" [0-9]*:N:0:[^ $]*", "",H)
		H=H " " READ ":N:0:" I
		print H 
		print R
		H=""
  	} else {
		print $1
	}
}