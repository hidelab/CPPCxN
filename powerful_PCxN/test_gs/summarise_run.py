
import pprint


filenames_1 = ["job_stats.txt", "job_stats2.txt"]

dct = {"network":0, "filename":0, "jobname":0, "part":0, "gs":0, "jobnumber":0, "taskid":0, "ru_wallclock":0, "ru_utime":0, "ru_stime":0, "cpu":0, "mem":0, "maxvmem":0, "start_time":0, "end_time":0, "exit_status":0}

# print("n0etwork\tfilename\tpart\tgs\tjobnumber\ttaskid\tru_wallclock\tru_utime\tru_stime\tcpu\tmem\tmaxvmem\tstart_time\tend_time")
print("\t".join("{}".format(k) for k, v in dct.items()))

for i, item in enumerate(filenames_1):
	for line in open(item, 'rU').readlines()[1:]:
		if line.startswith("=="):
			# pprint.pprint((dct.values()))
			# print(pprint.pformat((dct)))
			print("\t".join("{}".format(v) for k, v in dct.items()))
			dct = {x: 0 for x in dct}
			continue
		#assign filename 
		dct["filename"] = item

		f = line.split(" ", 1)
		prop = f[0]
		val = f[1].strip()
		if prop in dct.keys():
			# print(prop)
			dct[prop] = val
			# print(dct[prop])
		else:
			continue

		# print(dct)
		# print(dct["jobname"])
		jobname = dct['jobname']
		# assign what part it is
		if jobname.find('_1_') != -1:
			dct["part"] = 1
			# print(dct["part"])
		elif jobname.find('_2_') != -1:
			dct["part"] = 2
		elif jobname.find('_0_') != -1:
			dct["part"] = 0
		elif jobname.find('_3_') != -1:
			dct["part"] = 3

		#assign network  type
		if jobname.find('PDxN') != -1:
			dct["network"] = "PDxN"
		elif jobname.find('PCxN') != -1:
			dct["network"] = "PCxN"

		#assign number of gs
		if jobname.find('v800') != -1:
			dct["gs"] = 800
		elif jobname.find('v1200') != -1:
			dct["gs"] = 1200
		elif jobname.find('v1') != -1:
			dct["gs"] = 10
			# print(dct["part"])
		elif jobname.find('v2') != -1:
			dct["gs"] = 20
		elif jobname.find('v3') != -1:
			dct["gs"] = 50
		elif jobname.find('v4') != -1:
			dct["gs"] = 100
		elif jobname.find('v5') != -1:
			dct["gs"] = 200
		elif jobname.find('v6') != -1:
			dct["gs"] = 400
		elif jobname.find('v7') != -1:
			dct["gs"] = 1000
		elif jobname.find('v8') != -1:
			dct["gs"] = 1473



