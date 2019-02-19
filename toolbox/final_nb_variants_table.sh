grep ^# -cv */*reports/*.final.vcf | tr "/" " " | sed s/_env_/:/ | sed s/_sh/:/ | cut -d":" -f1,2,4 --output-delimiter=" " | column -t
