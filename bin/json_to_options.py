#!/usr/bin/env python
import json
import sys
from pprint import pprint

jdata = open(sys.argv[1])

data = json.load(jdata)

def keyvalue(obj,key=""):
	#print(type(obj))
	if type(obj) == dict:
	#if type(obj) is dict:
	#if isinstance(obj, dict):
		for k, v in obj.items():
			if hasattr(v, '__iter__'):
				keyvalue(v,k)
			else:
				if v == None or type(v) == bool:
				#if v == None or type(v) is bool:
					print('--'+k)
				else:
					#print ('--%s=%s') % (k, '"'+str(v).strip()+'"')
					#print ('--1=2') #% (k, '"'+str(v).strip()+'"')
					print('--'+k+'="'+str(v).strip()+'"')
	elif type(obj) == list:
	#elif type(obj) is list:
	#elif isinstance(obj, list):
		#list_concat=None
		list_concat=""
		for v in obj:
			sep=""
			if hasattr(v, '__iter__'):
				keyvalue(v)
			else:
				if list_concat != "":
					sep=","
				if v != None:
					list_concat=list_concat+sep+str(v).strip()
		if list_concat == None:
			print("--"+key)
		else: #if list_concat != "":
			print("--"+key+'="'+list_concat+'"')
	else:
		print('obj'+obj)

keyvalue(data)

jdata.close()
