#!/usr/bin/env python3

import os
import glob

current_dir=os.getcwd()


def main():
	replays = []
	with open("Configs/Runs.txt") as fp:
		line = fp.readline()
		while line:
			if len(line)==5:
				replays.append(line.strip())
			line = fp.readline()


	for replay in replays:
		for name in glob.glob(current_dir+'/Replay_'+replay+'/Data/*/*psa_*'):
			directory_name = os.path.dirname(name)

			file_name = os.path.basename(name)			

			if file_name.startswith('psa_0'):
				print(directory_name, ' name :', file_name)
				try:
					os.rename(name, directory_name+'/SRM_AGATA_'+file_name)
				except PermissionError:
					print("Permission erro for :", name)
			

if __name__ == '__main__':
	main()
