#!/usr/bin/env python3

import pandas as pd
import pathlib
from pprint import pprint
import re
import sys


common_size = { 'Bos_taurus' : 2870000000,
				'Homo_sapiens' : 3100000000,
				'Zea_mays' : 2800000000}

samples = pd.read_table(config['samples'], comment='#').set_index('sample', drop=False)

##############################
# METHOD
##############################

def complexExpand(wildcards):
	pattern = "subsampling/{sample}/{sample}-covx.bam.stats"
	samplesPart_list = [ s.replace('cov','{cov}') for s in expand(pattern, sample=samples.index)]
	
	full = list()
	for i in range(len(samples)):
		cov_list = samples.loc[samples.index[i], 'cov'].split(',')
		full = full + expand(samplesPart_list[i], cov=cov_list)
		
	return full

def get_TotLen(stat_file):
	pattern1 = 'SN'
	pattern2 = 'total length.*?(\d+)'
	SN_list = list()
	SN_pass = False
	with open(stat_file,'rt') as fh:
		for line in fh:
			# ~ pprint('Line: '+line)
			if re.search(pattern1,line):
				SN_list.append(line,)
				SN_pass = True
			elif (SN_pass):
				break
			else:
				continue
	
	# ~ pprint(SN_list)
	for SN in SN_list:
		extract = re.search(pattern2,SN)
		if extract:
			return int(extract.group(1))

def get_Cigar(stat_file):
	pattern1 = 'SN'
	pattern2 = '\(cigar\).*?(\d+)'
	SN_list = list()
	SN_pass = False
	with open(stat_file,'rt') as fh:
		for line in fh:
			# ~ pprint('Line: '+line)
			if re.search(pattern1,line):
				SN_list.append(line,)
				SN_pass = True
			elif (SN_pass):
				break
			else:
				continue
	
	# ~ pprint(SN_list)
	for SN in SN_list:
		extract = re.search(pattern2,SN)
		if extract:
			return int(extract.group(1))

def get_Size(value):
	''' Retrieve genome size, either directly from file as numerical value or
		from built-in dict with species name keyword
	'''
	# pprint(value)
	if re.search('[^\d]+',value):
		if value in common_size.keys():
			return int(common_size[value])
		else:
			sys.exit('Error in get_Size ! '+value+' is not an accepted keyword.')
	elif re.search('\d+',value):
		return int(value)
	else:
		pprint('Value type not expected')

def get_Bam(wildcards):
	''' Get bam file
	'''
	return samples.loc[wildcards.sample, "path"]

def get_Ratio(wildcards, resources):
	''' Calculate coverage ratio to produce asked coverage from base bam
	'''
	size = get_Size(samples.loc[wildcards.sample, "size"])
	cov = float(wildcards.cov)
	cov_base = float(resources.Cigar/size)
	ratio = format(cov/cov_base, '.4f')
	return ratio

def get_threads(rule, default):
	cluster_config = snakemake.workflow.cluster_config
	if rule in cluster_config and "threads" in cluster_config[rule]:
		return cluster_config[rule]["threads"]
	elif "default" in cluster_config and "threads" in cluster_config["default"]:
		return cluster_config["default"]["threads"]
	else:
		return default

##############################
# RULES
##############################

rule all:
	''' Main rule to produce default output
	'''
	input:
		'summary/reports.tsv'

rule base_bam:
	''' Rename starting bam to organize output
	'''
	input:
		bam = get_Bam
	output:
		bam = "mapping/{sample}.bam"
	shell:
		"ln -s {input.bam} {output.bam}"

rule base_index:
	''' Produce bam index for starting bam
	'''
	input:
		bam = "mapping/{sample}.bam"
	output:
		bai = "mapping/{sample}.bam.bai"
	conda:
		'envs/samtools-env.yaml'
	threads:
		get_threads('base_index',10)
	shell:
		"samtools index -@ {threads} {input.bam}"

rule base_stats:
	''' Produce stats file from samtools with base bam
	'''
	input:
		bam = "mapping/{sample}.bam",
		bai = "mapping/{sample}.bam.bai"
	output:
		stats = "mapping/{sample}.bam.stats"
	conda:
		'envs/samtools-env.yaml'
	threads:
		get_threads('base_stats',4)
	shell:
		"samtools stats -@ {threads} {input.bam} > {output.stats}"

rule subsampling:
	''' Subsample base bam with various cov
	'''
	input:
		bam = "mapping/{sample}.bam",
		stats = "mapping/{sample}.bam.stats"
	output:
		bam = "subsampling/{sample}/{sample}-{cov}x.bam"
	conda:
		'envs/samtools-env.yaml'
	resources:
		TotLen = lambda wildcards, input: get_TotLen(input.stats),
		Cigar = lambda wildcards, input: get_Cigar(input.stats)
	params:
		ratio = lambda wildcards, resources: get_Ratio(wildcards, resources)
	threads:
		get_threads('subsampling',10)
	shell:
		"samtools view -@ {threads} -s {params.ratio} -b {input.bam} > {output.bam}"
		
rule index:
	''' Produce bam index for starting bam
	'''
	input:
		bam = "subsampling/{sample}/{sample}-{cov}x.bam"
	output:
		bai = "subsampling/{sample}/{sample}-{cov}x.bam.bai"
	conda:
		'envs/samtools-env.yaml'
	threads:
		get_threads('index',4)
	shell:
		"samtools index -@ {threads} {input.bam}"

rule stats:
	''' Produce stats file from samtools for subsample file
	'''
	input:
		bam = "subsampling/{sample}/{sample}-{cov}x.bam",
		bai = "subsampling/{sample}/{sample}-{cov}x.bam.bai"
	output:
		stats = "subsampling/{sample}/{sample}-{cov}x.bam.stats"
	conda:
		'envs/samtools-env.yaml'
	threads:
		get_threads('stats',4)
	shell:
		"samtools stats -@ {threads} {input.bam} > {output.stats}"

rule reports:
	''' Summary of subsampling manip
	'''
	input:
		complexExpand
	output:
		'summary/reports.tsv'
	run:
		with open(str(output),'wt') as fh:
			print("\t".join(['Sample','Asked Coverage','Given Genome Size','Total Length','Cigar base mapped', 'Sequenced Coverage','Mapped coverage' ]), file=fh)
			for file in input:
				match = re.search('subsampling\/(\w+)\/(\w+)-(\d+)x\.bam\.stats',file)
				sample = match.group(1)
				cov = match.group(3)
				size = get_Size(samples.loc[sample, "size"])
				cigar = get_Cigar(file)
				TotLen = get_TotLen(file)
				print("\t".join([sample,str(cov),str(size),str(TotLen),str(cigar), str(format(cigar/size,'0.4f')), str(format(TotLen/size,'0.4f')) ]), file=fh)
