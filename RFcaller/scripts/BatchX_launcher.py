#!/usr/bin/python3
import json
import os
import subprocess

# BatchX CPUs
threads = os.environ['BX_VCPUS']

# Parse json
with open("/batchx/input/input.json", "r") as inputFile:
	inputJson = inputFile.read()
	parsedJson = json.loads(inputJson)

# Make our input directory
os.mkdir("/input")

# Make BAM's directories
## Normal
os.mkdir("/input/normal")
normalBam = "/input/normal/" + parsedJson["normal"]["normalBam"].split("/")[-1]
os.symlink(parsedJson["normal"]["normalBam"], normalBam)
if "normalIndex" in parsedJson["normal"]:
	normalIndex = "/input/normal/" + parsedJson["normal"]["normalIndex"].split("/")[-1]
	os.symlink(parsedJson["normal"]["normalIndex"], normalIndex)
## Tumor
os.mkdir("/input/tumor")
tumorBam = "/input/tumor/" + parsedJson["tumor"]["tumorBam"].split("/")[-1]
os.symlink(parsedJson["tumor"]["tumorBam"], tumorBam)
if "tumorIndex" in parsedJson["tumor"]:
	tumorIndex = "/input/tumor/" + parsedJson["tumor"]["tumorIndex"].split("/")[-1]
	os.symlink(parsedJson["tumor"]["tumorIndex"], tumorIndex)

# Change genome files
os.mkdir("/input/genome")
genome = "/input/genome/" + parsedJson["genome"]["genomeFasta"].split("/")[-1]
os.symlink(parsedJson["genome"]["genomeFasta"], genome)
if "genomeIndexFai" in parsedJson["genome"]:
	genomeIndexFai = "/input/genome/" + parsedJson["genome"]["genomeIndexFai"].split("/")[-1]
	os.symlink(parsedJson["genome"]["genomeIndexFai"], genomeIndexFai)
if "genomeIndexGzi" in parsedJson["genome"]:
	genomeIndexGzi = "/input/genome/" + parsedJson["genome"]["genomeIndexGzi"].split("/")[-1]
	os.symlink(parsedJson["genome"]["genomeIndexGzi"], genomeIndexGzi)

# Check dbSNP
if parsedJson["dbSNP"] == "Other":
	dbSNP = "/home/databases/dbSNP/" + parsedJson["dbSNPuser"]["dbSNPfile"].split("/")[-1]
	os.symlink(parsedJson["dbSNPuser"]["dbSNPfile"], dbSNP)
	if "dbSNPindex" in parsedJson["dbSNPuser"]:
		dbSNPindex = "/home/databases/dbSNP/" + parsedJson["dbSNPuser"]["dbSNPindex"].split("/")[-1]
		os.symlink(parsedJson["dbSNPuser"]["dbSNPindex"], dbSNPindex)
elif parsedJson["dbSNP"] == "hg19":
	dbSNP = "hg19"
elif parsedJson["dbSNP"] == "hg38":
	dbSNP = "hg38"

# Check PoN
if "PoN" in parsedJson:
	if parsedJson["PoN"] == "Other":
		PoN = parsedJson["PoNuser"]["PoNfile"]
	else:
		PoN = parsedJson["PoN"]
else:
	PoN = "NULL"

# Check regions
if "regions" in parsedJson:
	regions = parsedJson["regions"]
else:
	regions = "NULL"

# Check Ploidy file
if ("assembly" in parsedJson) and ("assemblyTag" in parsedJson["assembly"]):
	ploidyFile = parsedJson["assembly"]["assemblyTag"]
	if ploidyFile == "Other":
		ploidyFile = parsedJson["assembly"]["assemblyFile"]
else:
	ploidyFile = "NULL"

# Check defaults
## Contamination
if "contamination" in parsedJson:
	contamination = parsedJson["contamination"]
else:
	contamination = 0.05

## Tumor coverage for SNVs
if "TD_cov_SNV" in parsedJson:
	TD_cov_SNV = parsedJson["TD_cov_SNV"]
else:
	TD_cov_SNV = 7
## Tumor coverage for INDELs
if "TD_cov_INDEL" in parsedJson:
	TD_cov_INDEL = parsedJson["TD_cov_INDEL"]
else:
	TD_cov_INDEL = 7
## Normal coverage for SNVs
if "ND_cov_SNV" in parsedJson:
	ND_cov_SNV = parsedJson["ND_cov_SNV"]
else:
	ND_cov_SNV = 7
## Normal coverage for INDELs
if "ND_cov_INDEL" in parsedJson:
	ND_cov_INDEL = parsedJson["ND_cov_INDEL"]
else:
	ND_cov_INDEL = 7
## Minimum number of mut reads in tumor for SNVs
if "TD_mut_SNV" in parsedJson:
	TD_mut_SNV = parsedJson["TD_mut_SNV"]
else:
	TD_mut_SNV = 3
## Minimum number of mut reads in tumor for INDELs
if "TD_mut_INDEL" in parsedJson:
	TD_mut_INDEL = parsedJson["TD_mut_INDEL"]
else:
	TD_mut_INDEL = 4
## Max number of mut reads in normal for SNVs
if "ND_mut_SNV" in parsedJson:
	ND_mut_SNV = parsedJson["ND_mut_SNV"]
else:
	ND_mut_SNV = 3
## Max number of mut reads in normal for INDELs
if "ND_mut_INDEL" in parsedJson:
	ND_mut_INDEL = parsedJson["ND_mut_INDEL"]
else:
	ND_mut_INDEL = 2
## Window size to analyze INDELs
if "ND_window" in parsedJson:
	ND_window = parsedJson["ND_window"]
else:
	ND_window = 10
## SNV_threshold
if "SNV_threshold" in parsedJson:
	SNV_threshold = parsedJson["SNV_threshold"]
else:
	SNV_threshold = 10.726
## INDEL_threshold
if "INDEL_threshold" in parsedJson:
	INDEL_threshold = parsedJson["INDEL_threshold"]
else:
	INDEL_threshold = 32.1418
## polyINDEL_threshold
if "polyINDEL_threshold" in parsedJson:
	polyINDEL_threshold = parsedJson["polyINDEL_threshold"]
else:
	polyINDEL_threshold = 0.7723

# Search for incompatibilities between versions
if ((parsedJson["dbSNP"] == "hg19" and ploidyFile == "GRCh38") or (parsedJson["dbSNP"] == "hg38" and ploidyFile == "GRCh37") or (parsedJson["dbSNP"] == "hg19" and ploidyFile == "GRCm38") or (parsedJson["dbSNP"] == "hg19" and ploidyFile == "GRCm39") or (parsedJson["dbSNP"] == "hg38" and ploidyFile == "GRCm38") or (parsedJson["dbSNP"] == "hg38" and ploidyFile == "GRCm39")):
	raise ValueError("dbSNP and assemblyTag fields must correspond to the same genome versions.")

if (parsedJson["dbSNP"] == "Other" and (ploidyFile == "None" or ploidyFile == "NULL")):
	raise ValueError("If you are using your own dbSNP, you must choose between: GRCh37, GRCh38, GRCm38, GRCm39 or provide your own ploidy_file.")

# Launch RFcaller
try:
	cmd = "RFcaller" \
		+ " -@ " + str(threads) \
		+ " --normalBam " + normalBam \
		+ " --normal " + parsedJson["normal"]["normalName"] \
		+ " --tumorBam " + tumorBam \
		+ " --tumor " + parsedJson["tumor"]["tumorName"] \
		+ " --output " + parsedJson["outputName"] \
		+ " --genome " + genome \
		+ " --dbSNP " + dbSNP \
		+ " --contamination " + str(contamination) \
		+ " --TD_cov_SNV " + str(TD_cov_SNV) \
		+ " --TD_cov_INDEL " + str(TD_cov_INDEL) \
		+ " --ND_cov_SNV " + str(ND_cov_SNV) \
		+ " --ND_cov_INDEL " + str(ND_cov_INDEL) \
		+ " --TD_mut_SNV " + str(TD_mut_SNV) \
		+ " --TD_mut_INDEL " + str(TD_mut_INDEL) \
		+ " --ND_mut_SNV " + str(ND_mut_SNV) \
		+ " --ND_mut_INDEL " + str(ND_mut_INDEL) \
		+ " --ND_window " + str(ND_window) \
		+ " --SNV_threshold " + str(SNV_threshold) \
		+ " --INDEL_threshold " + str(INDEL_threshold) \
		+ " --polyINDEL_threshold " + str(polyINDEL_threshold)
	if (PoN != "NULL"):
		cmd += " --PoN " + PoN
	if (regions != "NULL"):
		cmd += " --regions " + regions
	if (ploidyFile != "NULL"):
		cmd += " --ploidy_file " + ploidyFile
	print(cmd, flush=True)
	subprocess.check_call(cmd, shell=True)
except subprocess.CalledProcessError as e:
		print(e)
		exit(e.returncode)

# Output json file
outputJson = {}
outputJson["vcf"] = "/batchx/output/mutations_" + parsedJson["outputName"] + ".vcf"
outputJson["log"] = "/batchx/output/RFcaller.log"

with open('/batchx/output/output.json', 'w+') as json_file:
	json.dump(outputJson, json_file)
