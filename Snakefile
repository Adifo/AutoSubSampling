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
	pattern = "subsampling/{sample}/{sample}-covx.bam"
	samplesPart_list = [ s.replace('cov','{cov}') for s in expand(pattern, sample=samples.index)]
	pprint(samplesPart_list)
	
	full = list()
	for i in range(len(samples)):
		cov_list = samples.loc[samples.index[i], 'cov'].split(',')
		full = full + expand(samplesPart_list[i], cov=cov_list)
		
	pprint(full)
	return full

def get_TotLen(input):
	pattern1 = 'SN'
	pattern2 = 'total length.*?(\d+)'
	SN_list = list()
	SN_pass = False
	with open(input.stats,'rt') as fh:
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

def get_Cigar(input):
	pattern1 = 'SN'
	pattern2 = '\(cigar\).*?(\d+)'
	SN_list = list()
	SN_pass = False
	with open(input.stats,'rt') as fh:
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
	pprint(value)
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
		complexExpand

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

rule base_stat:
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
		get_threads('base_stat',4)
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
		TotLen = lambda wildcards, input: get_TotLen(input),
		Cigar = lambda wildcards, input: get_Cigar(input)
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
