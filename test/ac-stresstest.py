#!/usr/bin/env python3
# vim: tabstop=4:ai:

import argparse
import requests
import sys
import numpy as np
import time


achost = None
timeout = 300

def oneid(unid,phrange,vlist):
	for v in vlist:
		for ph in np.arange(*phrange):
			before = time.time()
			res = oneshot(unid,ph,v)
			after = time.time()
			print(unid,ph,v,res,after-before)

def perr(*args,**kwargs):
	print(*args,file=sys.stderr,**kwargs)

def oneshot(unid,ph,v):
	try:
		r = requests.get(f"{achost}/calculate_charges/{unid}?ph={ph}&alphafold_prediction_version={v}",timeout=timeout) 
		if r.status_code != 200:
			perr(r.json())
			return f'HTTP_{r.status_code}'

		if r.json()['status'] != 'partial atomic charges successfully calculated':
			perr(r.json())
			return 'CALC_ERROR'

	except Exception as e:
		perr(e)
		return 'HTTP_EXCEPTION'

	try:
		rt = requests.get(f"{achost}/download_file/{r.json()['ID']}/txt")
		if rt.status_code != 200:
			perr(r.json())
			return f'HTTP_{r.status_code}'

		resp = rt.text.split('\n')
		if len(resp) != 3:
			return 'BAD_TXT'

		chg = np.array(list(map(float,resp[1].split())))
		if np.min(chg) < -3.0 or np.max(chg) > 3.0:
			return 'BAD_CHARGE'

	except Exception as e:
		perr(e)
		return 'HTTP_EXCEPTION'

	return 'OK'


if __name__ == '__main__':

	p = argparse.ArgumentParser(description='Stress test Alpha Charge server with a series of API calls and check results for sanity.')

	p.add_argument('ids',help='File with one Uniprot ID per line')
	p.add_argument('-v','--versions',help='Comma separated list of Alpha Fold prediction versions',default='3,4',required=False)
	p.add_argument('-p','--phrange',help='Range of pH to try as MIN,MAX,STEP',default='6.0,8.1,1',required=False)
	p.add_argument('-s','--server',help='URL prefix of the server',default='http://147.251.115.173/',required=False)

	a = p.parse_args()
	achost = a.server

	phrange = list(map(float,a.phrange.split(',')))
	vlist = tuple(map(int,a.versions.split(',')))
	with open(a.ids) as ids:
		for unid in ids:
			unid = unid.rstrip()
			oneid(unid,phrange,vlist)
	
