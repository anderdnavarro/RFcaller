#!/usr/bin/env python3

from sys import stdin

mutID = 1
with stdin as vcf:
	for line in vcf:
		if line.startswith("#"):
			print(line.strip())
			continue
		else:
			columns = line.strip().split("\t")
			print("\t".join(columns[0:2]), f"{columns[2]}_{mutID}", "\t".join(columns[3:]), sep = "\t")
			mutID += 1
			continue