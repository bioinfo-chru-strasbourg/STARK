#!/usr/bin/env python
import json
import sys
from pprint import pprint

jdata = open(sys.argv[1])

data = json.load(jdata)

def keyvalue(obj,key=""):
	if type(obj) == dict:
		for k, v in obj.items():
			if hasattr(v, '__iter__'):
				keyvalue(v,k)
			else:
				print '--%s=%s' % (k, '"'+str(v).strip()+'"')
	elif type(obj) == list:
		list_concat=""
		for v in obj:
			sep=""
			if hasattr(v, '__iter__'):
				keyvalue(v)
			else:
				if list_concat != "":
					sep=","
				list_concat=list_concat+sep+'"'+str(v).strip()+'"'
		if list_concat != "":
			print "--"+key+"="+list_concat
	else:
		print obj

keyvalue(data)

jdata.close()
