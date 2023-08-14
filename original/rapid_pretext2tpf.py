#!/usr/bin/env python
#Alan Tracey 2021


'''
MIT License
 
Copyright (c) 2021 Genome Research Ltd.
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


'''

import sys
import argparse
import pyfastaq
import subprocess
from datetime import datetime

'''
#rapid_pretext2tpf_XL.py has now been written for large fragmented genomes!!  For small genomes with few gaps this program is fine and quicker to use.
#A known issue is that we are filtering away small agp fragments.  If we have a large highly fragmented genome
#requiring lots of breaking, we may legitimately create small fragments which could get discarded by this method.
#For full curation, we know the precise break coordinates which enables us to solve this problem with certainty.  
#In rapid we know where breaks might occur, then we have to figure out where the most likely coordinate is but
#we don't know with certainty.  There is always the 
#possibility that we legitimately create a small fragment in a scaffold which has approx the same size as an 
#artefact.  
#Attempted to address this by calculating the smallest genuine fragment (contig) lengths for
#each scaffold.  Take agp lines for anything that is
#1) not small (>x*texel)
#2) > 0.25 texel length smaller than the smallest contig size.
#Any discarded agp lines are reported to help curators address any errors like this.  It is recommended that
#rapid_pretext2tpf.py is not used on large, highly fragmented genomes that need to be broken into many very small
#fragments.  Many breaks should be ok, but many breaks to create many very tiny pieces may need some manual intervention.
#It is possible that in scenarios above, very small tpf chunks may remain unbroken (ie they may remain in situ in their original
#placement and not follow the intended rearrangement in pretext.  These would then need to be manually found (by generating a 
#new pretext and inspecting it.  It is hoped such cases would be rare.


'''

#Global parameters
#=================

netsize=3	#throw a wide net to find tpf coords matching agp breaks - but not too wide (in testing >4 misses breaks in highly fragmented genomes - this just means that the tpfchunks stay together rather than splitting fully to match the agp)
lowcutoff=1.5	#how many texels smaller than the smallest contig in this scaffold we are prepared to go down to look for small contigs (if too small we may be misled by agp fragment artifacts).  Had to raise to 1.2 for idMelMell2_1
#NB there is no wholesale texel length cutoff, lowcutoff performs this role more intelligently on a per scaffold basis
sex=["X","Y","Z","W"]
prefix="R"
borderlen=80
errors={}
	

def append_dict(k,v,d):

	if k not in d:
		d[k]=[v]
	else:
		d[k].append(v)
	return d


#Takes dictionary as input
def genome_size(tpfdict):

	total=0
	totals={}
	for k,v in tpfdict.items():
		for line in v:
			x=line.split()
			if not "GAP" in line:
				scaff=x[1].split(":")[0]
				hi=int(x[1].split("-")[1])
				lo=int(x[1].split(":")[1].split("-")[0])
				amount=(hi-lo)+1
				if scaff not in totals:
					totals[scaff]=amount
				else:
					totals[scaff]+=amount
				
	for k,v in totals.items():
		total+=v
		
	texel=small=int(round(total/32768,0))

	return total,texel


#Takes line list as input	
def genome_size2(lines):

	total=0
	totals={}
	for line in lines:
		x=line.split()
		if not "GAP" in line:
			scaff=x[1].split(":")[0]
			hi=int(x[1].split("-")[1])
			lo=int(x[1].split(":")[1].split("-")[0])
			amount=(hi-lo)+1
			if scaff not in totals:
				totals[scaff]=amount
			else:
				totals[scaff]+=amount
				
	for k,v in totals.items():
		total+=v

	return total


def make_date_string():
	v=time.ctime()	#Human month eg Oct
	v=str(datetime.now())	#date time string	
	r=v.split()	
	return r[0].replace("-","")


def complement_scaffold(scaff_line_list):
	x=scaff_line_list[::-1]
	newlines=[]
	for line in x:
		if "PLUS" in line:
			newlines.append(line.replace("PLUS","MINUS"))	
		elif "MINUS" in line:
			newlines.append(line.replace("MINUS","PLUS"))
		else:
			newlines.append(line)	

	return newlines
	

def scaffs_from_agp(agp,fragsize,scafflens,ctg_lengths,texel):

	sscaffdict={}
	discards={}
	agplines={}
	tagdict={}
	sexchrm={}
	allsex=[]
	nohap_sex=""
	scaff_with_haplo=[]
	with open(agp, 'r') as f:
		for line in f:
			#ßprint(line)
			if line[0] != "#"and line !="\n":
				x=line.strip().split()
				superscaff=x[0]
				cumulativelow=x[1]
				cumulativehigh=x[2]
				frag=int(cumulativehigh)-int(cumulativelow)
				scaff=x[5]
				orientation=x[8]
				entry=x[4]
				if entry != "U":	#is not a gap
					low=int(x[6])
					high=int(x[7])
					vals=[scaff,low,high,orientation,cumulativelow,cumulativehigh]
					append_dict(superscaff,line,agplines)
					##print(superscaff)
					if frag > 10*texel:	#always take agp frags above this size - never get an artefact bigger than 10 texels
						append_dict(superscaff,vals,sscaffdict)
					else:	#If fragment is small...
						if frag > min(ctg_lengths[scaff])-lowcutoff*texel:	#take if small but bigger than smallest contig by a margin
							#print(f'{frag} > {min(ctg_lengths[scaff])}-{lowcutoff}*{texel}')
							append_dict(superscaff,vals,sscaffdict)
						else:
							#print('Here!')
							append_dict(scaff,low,discards)		#discard everything else
					#print(sscaffdict)
					#Get the tags
					s=[]
					if len(x)>9:
						unqkey="/".join([scaff,str(low),str(high)])
						##print(unqkey, x[8:])
						for i in range(len(x[9:])):
							s.append(x[9:][i].upper())
						if "HAPLOTIG" in s and "UNLOC" in s:
							msg=scaff+" is both \'Haplotig\' and \'unloc\' - needs to be one or the other\n"
							append_dict("hap_unloc",msg,errors)
						if "HAPLOTIG" in s:
							scaff_with_haplo.append(superscaff)
						if "W" in s:
							nohap_sex="W"
						if "Y" in s:
							nohap_sex="Y"	
						for sx in sex:
							if sx in s:
								if sx not in allsex:
									allsex.append(sx)
								if sx not in sexchrm:
									sexchrm[sx]=[superscaff]	#Starting a check to see if eg Z chr is referenced in >1 chrm
								else:
									if superscaff not in sexchrm[sx]:
										sexchrm[sx].append(superscaff)	
						
						append_dict(unqkey,s,tagdict)

	if len(allsex)>2:
		msg = "Too many sex chromosomes - you have "+" ".join(allsex)
		append_dict("xssex",msg,errors)
	elif len(allsex)==2:
		if "Z" in allsex and "X" in allsex:
			msg = "Heterogametic sex but bad combination found "+" ".join(allsex)
			append_dict("badsex",msg,errors)
		elif "W" in allsex and "X" in allsex:
			msg = "Heterogametic sex but bad combination found "+" ".join(allsex)
			append_dict("badsex",msg,errors)
		elif "Z" in allsex and "Y" in allsex:
			msg = "Heterogametic sex but bad combination found "+" ".join(allsex)
			append_dict("badsex",msg,errors)
		else:
			for k,v in sexchrm.items():
				if v[0] in scaff_with_haplo:
					append_dict("hetero_sex_haplo",k,errors)

	##print(tagdict)						
	multisex2={}
	for k,v in sexchrm.items():
		if len(v)>1:
			length=str(len(v))
			msg=k+" chromosome is referenced in "+length+" AGP chromosomes:\n"
			append_dict("multiple_sex",msg,errors)
			for i in v:
				append_dict("multiple_sex",i,errors)
		for i in v:
			append_dict(i,k,multisex2)
	
	for k,v in multisex2.items():
		if len(v)>1:
			msg=k+" chromosome is referenced as both "+v[0]+" and "+v[1]+" in AGP:\n"
			append_dict("multiple_sex2",msg,errors)

	print(discards)
	return sscaffdict, discards, agplines, tagdict, sexchrm


def tpf_sanity(tpfout):

	cmd="perl /software/grit/projects/vgp_curation_scripts/test_tpf_sanity.pl -scafflevel "+tpfout
	pyfastaq.utils.syscall(cmd)


def location(f):

	cmd="readlink -f "+f
	b=subprocess.getoutput(cmd)
	
	return b

	
#Gets tpf lines and checks for typos in coordinates (some coordinates may be added manually so we check them)
def parse_tpf(tpf):

	prevscaff=""
	result={}
	coordtest={}	#Checking for obvious coord typos in the input tpf
	errors2=[]
	with open(tpf,'r') as f:

		for line in f:
			if not "gap" in line.lower():
				x=line.split()
				scaff=x[2]
				lo=x[1].split(":")[1].split("-")[0]
				hi=x[1].split(":")[1].split("-")[1]
				if scaff in coordtest:
					if int(lo)<=max(coordtest[scaff]):	#If coord is less than previous line
						if line not in errors2:
							errors2.append(line)
							
					if int(hi)<=max(coordtest[scaff]):
						if line not in errors2:
							errors2.append(line)
							
					if int(lo)>=int(hi):
						if line not in errors2:
							errors2.append(line)
							
				append_dict(scaff,int(hi),coordtest)
				if scaff not in result:
					result[scaff]=[line]
					prevscaff=scaff
				else:
					result[scaff].append(line)
					prevscaff=scaff
			else:
				result[prevscaff].append(line)

	return result, errors2


def report_errors(errors2):

	if len(errors2)>0:
		#print("coordinate errors detected amongst the following tpf lines - presumed typos:\n")
		for e in errors2:
			print(e)
		sys.exit()


def check_components(lines):

	check=[]
	for l in lines:
		if not "GAP" in l:
			comp=l.split()[1]
			check.append(comp)
			
	return check


def components_from_dict(tpf):

	check=[]
	for k,v in tpf.items():
		for l in v:
			if not "GAP" in l:
				comp=l.split()[1]
				check.append(comp)

	return check


#The coordinates where our agp fragments terminate
def agp_dividers(agpdict):

	tdivs={}
	result={}
	for k,v in agpdict.items():
		for i in v:
			tscaff=i[0]
			tmax=int(i[2])
			append_dict(tscaff,tmax,tdivs)

	for k, v in tdivs.items():
		result[k]=sorted(v)
	#print(result)
	return result


#dividers is agp dividing coordinates - here we add the agp coord into our results key with scaff name, then add the closest tpf coord that meets our parameterised requirements
def nearest(tpfdict,dividers,fragsize,discards,scafflen):
	#waypoint
	closest={}	#scaff:agpdiv key - val is closest tpf max
	results={}
	for k,v in dividers.items():
		for div in v:
			for line in tpfdict[k]:
				if not "GAP" in line:
					tmax=int(line.split()[1].split("-")[1]) #scaff_end
					##print(k,line.strip(),div,tmax)	
					pre=k+":"+str(div)

					if abs(tmax-div) < netsize*fragsize:			#throw a wide net but not too wide (in testing >4 misses breaks in highly fragmented genomes - this just means that the tpfchunks stay together rather than splitting fully to match the agp)
						if pre not in closest:
							#print(f'{tmax} - {div} | {netsize} * {fragsize}')
							#print(pre, tmax)
							closest[pre]=tmax
						else:
							if abs(tmax-div)<abs(tmax-closest[pre]):	#but get the closest
								#print(f'{abs(tmax-div)}<{abs(tmax-closest[pre])}')
								#print(k)
								closest[pre]=tmax
	#[print(i) for i in closest]
	#waypoint
	agptpfdiscrep=[]	#What remains in this list are elements in tpf that need breaking
	for k,v in dividers.items():
		for div in v:
			if k+":"+str(div) not in closest:
				##print(k, div, "not found!")
				agptpfdiscrep.append(k+"\t"+str(div))
				if k in discards:
					for d in discards[k]:
						if div+1 == int(d):
							##print(k,div,"It's ok, this was discarded as fragment artefact")
							agptpfdiscrep.remove(k+"\t"+str(div))
	#[#print(x, y) for x, y in closest.items() if x.startswith('scaffold_1:')]
	corrected = add_break(agptpfdiscrep, tpfdict) # Introduce breaks per item in agptpfdiscrep
	[print(i, v) for i, v in closest.items() if i == 'scaffold_1:72267081']
	merged = {**closest, **corrected}

	for k,v in merged.items():
		base=k.split(":")[0]
		if v!=scafflen[base]:	#Filtering out unbroken scaffolds
			results[k]=v

	#for k,v in results.items():
		##print(k,v)	
	return results				


def add_break(to_break, tpfdict):
	
	corrected = {}
	[print(f'WARNING - Break missing from input tpf? {i}') for i in to_break]
	for x in to_break:
		item = x.split()
		corrected[':'.join(item)] = int(x.split()[1]) + 1 # Add gap component
		#print(f'\tADDING BREAK @ {x} -- {int(x.split()[1]) + 1}')

	return corrected

def lens(tpfdict):

	lens={}
	for k,v in tpfdict.items():
		for i in v:
			##print(k,i)
			if not "GAP" in i:
				hi=int(i.split()[1].split("-")[1])
				if k not in lens:
					lens[k]=hi
				elif hi > lens[k]:
					lens[k]=hi

	return lens


def breaktpf(tpfdict,breakpoint):

	unused=[]
	breaknew={}
	results={}
	for k,v in breakpoint.items():
		base=k.split(":")[0]
		if base not in breaknew:
			breaknew[base]=[v]
		else:
			breaknew[base].append(v)

	#for k,v in breaknew.items():
		#for coord in v:
			##print(k,coord)

	for k,v in tpfdict.items():
		for line in v:
			if not "GAP" in line:
				if line not in unused:
					unused.append(line)


	#Here we set up a dictionary of tpfchunks (all the places we can find to break the tpf from comparing to the agp chunks)
	for tscaff,v in tpfdict.items():
		iteration=1
		for line in v:
			##print(k,line.strip())
			if not "GAP" in line:
				pre=line.split()[1].split(":")[0]+"%"+str(iteration)
				tmax=int(line.split()[1].split("-")[1])
				if tscaff in breaknew:	#if scaff is broken
					for bscaff, breaks in breaknew.items():
						if tscaff==bscaff:
							for coord in breaks:
								##print(tscaff,coord,tmax)
								if coord==tmax:
									##print(tscaff,coord,tmax)
									if pre not in results:
										if line in unused:
											results[pre]=[line]
											unused.remove(line)
										iteration+=1
										pre=line.split()[1].split(":")[0]+"%"+str(iteration)	
									else:
										if line in unused:
											results[pre].append(line)
											unused.remove(line)
										iteration+=1
										pre=line.split()[1].split(":")[0]+"%"+str(iteration)
							
								else:
									if pre not in results:
										if line in unused:
											results[pre]=[line]
											unused.remove(line)
									else:
										if line not in results[pre]:
											if line in unused:
												results[pre].append(line)
												unused.remove(line)
										elif "GAP" in line:
											results[pre].append(line)	
			
				else:	#if scaff not broken
					iteration=1	#Capturing unbroken scaffolds
					##print("a",line)
					pre=line.split()[1].split(":")[0]+"%"+str(iteration)
					if pre not in results:
						results[pre]=[line]
						unused.remove(line)
					else:
						results[pre].append(line)
						unused.remove(line)
			else:
				if pre not in results:
					if line in unused:
						results[pre]=[line]
						unused.remove(line)
				else:
					if line in unused:
						results[pre].append(line)	#gap line
						unused.remove(line)
					elif "GAP" in line:
						results[pre].append(line.strip())	

	for k,v in results.items():
		if "GAP" in v[0]:
			del v[0]
		if "GAP" in v[-1]:
			del v[-1]

	test=[]

	for k,v in results.items():
		##print(k,v)
		for line in v:
			if not "GAP" in line:
				if line not in test:
					test.append(line)
				else:
					print("duplicate",line)

	return results	#scaffold_21%1 ['?\tscaffold_21:1-320171\tscaffold_21\tPLUS\n', 'GAP\tTYPE-2\t23\n', '?\tscaffold_21:320195-1446380\tscaffold_21\tPLUS\n']
	
	
def outputlist(tpfchunks,agpdict,tagdict):

	#Setting up agp scaff to tpfchunk key (agp,orientation,tpfchunk)
	#Scaffold_30:+:1#scaffold_32_ctg1%1
	#Scaffold_30:+:2#scaffold_52_ctg1%1
	#Scaffold_30:+:3#scaffold_51_ctg1%1

	check=[]
	joins=0

	##for k,v in tpfchunks.items():
		#if k not in check:
			#vals=[k,v]
			#check.append(vals)
	##print(tagdict)


	tmp=[]
	results_new={}
	outlines=[]
	tagged={}
	
	iteration=1
	

	for agk,v in agpdict.items():
		for i in v:
			##print("A",i)
			ascaff=i[0]
			ornt=i[3]
			alo=int(i[1])
			ahi=int(i[2])
			tgdictk = i[0]+"/"+str(i[1])+"/"+str(i[2])
			tmp.append(ascaff)
			agprefix=agk+":"+ornt+":"+ascaff+"%"+str(tmp.count(ascaff))
			##print(scaff,ornt,alo,ahi,prefix)
			for tchunk,tlines in  tpfchunks.items():
				##print(tchunk,tlines)
				##print(agprefix)
				tbase=tchunk.split("%")[0]
				tlo=int(tlines[0].split()[1].split(":")[1].split("-")[0])	#first line in the chunk
				##print(tlines[-1])
				thi=int(tlines[-1].split()[1].split(":")[1].split("-")[1])	#last line in the chunk
				if tbase==ascaff:
					##print(ascaff,ornt,alo,ahi,agprefix,tchunk,tlo,thi)
					##print(tlines)
					factor=0.7
					tlen=round((thi-tlo)*factor,2)
					adj=tlen*factor
					adjtlo=round(tlo+adj,2)
					adjthi=round(thi-adj,2)
					##print(tchunk,adjtlo,adjthi,alo,ahi,tlo,thi)
					##print(tlines)
					if adjtlo > alo and adjthi < ahi:
						t2akey=agk+":"+ornt+"#"+tchunk
						if tgdictk in tagdict:
							append_dict(t2akey,tagdict[tgdictk][0],tagged)	#attaching tags to the tpfchunk/agp unq key
						##print(k,tchunk,adjtlo,adjthi,alo,ahi,tlo,thi,ornt)
						if agk not in results_new:	#The first time we see this chromosome
							##print(k,tchunk,adjtlo,adjthi,alo,ahi,tlo,thi,ornt)
							if ornt=="-":
								results_new[agk]=[]	#set an empty list to add to
								for aaa in complement_scaffold(tlines):
									results_new[agk].append(aaa)
									
								#results_new[agk].append(complement_scaffold(tlines))	#complement that line and add it as first line in that chr
								
							else:
								results_new[agk]=[]	#set an empty list to add to
								for aaa in tlines:
									results_new[agk].append(aaa)	#don't complement that line but add it as the first line for that chr
									
						else:	
							gap="GAP\tTYPE-2\t200"	#before adding the next line, always add a gap
							results_new[agk].append(gap)
							joins+=1		#count joins
							##print(tchunk,tlines)
							if ornt=="-":
								for line in complement_scaffold(tlines):
									results_new[agk].append(line)
									
									##print(k,tchunk,line.strip())

							else:
								for line in tlines:
									results_new[agk].append(line)

	#for k,v in tpfchunks.items():
		##print(k,v)
	#for k,v in results_new.items():
		##print(k,v)
	#for k,v in tagged.items():
		##print(k,v)					
								
	iteration=0				
	for k,v in results_new.items():
		iteration+=1
		##print(k,v)
		for line in v:
			x=line.split()
			if not "GAP" in line:
				x[2]=prefix+str(iteration)
				newline="\t".join(x)
					##print(newline)
				outlines.append(newline)
			else:
				outlines.append(line)	#Putting original gaps back in

	return outlines, joins, tagged			


#Adds any tpf lines (eg shrapnel) missing from agp back in
def reinstate_lines(tpfdict,outlines,checkin,checkout):

	tmp=[]
	dups=[]
	finaloutlines=[]
	for i in checkin:	#i is eg scaffold_22:1693718-3211153 - checkin is all input regions
		if i not in checkout:
			for k,v in tpfdict.items():
				for line in v:
					if not "GAP" in line:
						if i in line:
							tmp.append(line.strip().split()[2])	#tmp == shrapnel

						
	shrap = sorted(list(set(sorted(tmp))))
	shrapnums=[]	#Going to use this to sort the final output
	for scaff in shrap:
		num=int(scaff.split("_")[1])
		if num not in shrapnums:
			shrapnums.append(num)
	sorted_shrapnums=sorted(shrapnums)
	for num in sorted_shrapnums:
		for scaff in shrap:
			scfnum=int(scaff.split("_")[1])
			if num==scfnum:
				for i in tpfdict[scaff]:
					region=i.strip().split()[1]
					if region in checkout:
						dups.append(region)	#Find duplicates we're going to be rarely adding back in
					outlines.append(i.strip())	#Shrapnel is added back in after sorting it on scaffnum
			

	for scaff in shrap:
		for i in tpfdict[scaff]:
			for d in dups:
				if d in i:
					for o in outlines:
						if i.strip() in o:
							outlines.remove(i.strip())	#Ensure these rare duplicates are removed
	

	final={}	#Build a new dict with sole purpose of removing trailing gap lines left over from removal of rare shrap duplicates
	prev=""
	for line in outlines:
		x=line.strip().split()
		if not "GAP" in line:
			scaff=x[2]
			append_dict(scaff,line,final)
			prev=scaff
		else:
			append_dict(scaff,line,final)
	
	for k,v in final.items():	#After we've removed duplicate shraps, we have trailing gaps to remove
		while "GAP" in v[-1]:
			del v[-1]
		while "GAP" in v[0]:
			del v[0]	
		for line in v:
			x=line.split()
			scaff=x[2]
			for d in dups:
				if scaff in d:
					x[2]=scaff+"_1"	#Add a suffix to say the scaffold has changed from the original scaffold
			newline="\t".join(x)
			finaloutlines.append(newline)
			
	#So far this hasn't happened, but it's conceivable that shrapnel dups could be removed from start and end of a scaffold, leaving them reinstated with multiple gaps inbetween
	#The below 50 or so lines of code inserted to hopefully deal with this, although no test case has yet come up.  Program runs fine
	#with or without this code.  For all test-cases so far, the below code is not called.  However it should work if this rare
	#situation ever crops up

	gaps=[]
	prevscaff=""
	table={}
	keystoupdate=[]
	iteration=1
	inter={}
	linenumber=1
	for line in finaloutlines:
		x=line.strip().split()
		if not "GAP" in line:
			scaff=x[2]
			gaps=[]
			append_dict(scaff,linenumber,inter)
			linenumber+=1
			
		else:
			gaps.append(x)
			linenumber+=1
			if len(gaps)>1:
				if scaff not in keystoupdate:
					keystoupdate.append(scaff)
					
	#keystoupdate are all scaffs that contain consecutive gap lines - places where broken scaffs have been reinstated incorrectly

	new={}
	for k,v in inter.items():
		a = [v[i+1]-v[i] for i in range(len(v)-1)]
		
		for i in a:
			if i>2:
				if k not in new:
					new[k]=[0]+a	#Here we're collecting distances between non-gap line numbers
					
	breakers={}	#The scaffs we want to assign a new ID to as they're fragments of the original, with line number of the new fragment startpoint			
	for k, v in new.items():
		for e, i in enumerate(v):
			if i >2:		#These are consecutive gap lines we are trying to catch - ie distance between non-gap lines is greater than 2
				if k not in breakers:
					breakers[k]=[inter[k][e]]
				else:
					breakers[k].append(inter[k][e])
			
	iteration=1
	newkeys=[]
	linenumber=0	
	lastline=[]
	finalfinal=[]	
	for line in finaloutlines:
		linenumber+=1
		x=line.strip().split()
		if not "GAP" in line:
			scaff=x[2]
			if scaff in breakers:
				newscaff=scaff+"_"+str(iteration)
				if linenumber in breakers[scaff]:
					iteration+=1
					newscaff=scaff+"_"+str(iteration)	#Start of a new fragment
					x[2]=newscaff
					newline="\t".join(x)
					finalfinal.append(newline)
					lastline.append(line.strip())
				else:
					x[2]=newscaff
					newline="\t".join(x)
					finalfinal.append(newline)
					lastline.append(line.strip())
			else:
				finalfinal.append(line.strip())
				lastline.append(line.strip())
		elif not "GAP" in lastline[-1]:			#We don't want to output multiple gaps - just take one
			finalfinal.append(line.strip())
			lastline.append(line.strip())

	return finalfinal


def write_output_tpf(outlinesfull,outfile):
 
	with open(outfile,'w') as fout:
		for k,v in outlinesfull.items():
			for line in v:
				fout.write(line+"\n")


def compare_scaff(tpfdict,agp):

	scaffs=[]
	probs=[]
	with open(agp, 'r') as f:
		for line in f:
			if line[0] != "#" and line !="\n":
				x=line.strip().split()
				scaff=x[5]
				entry=x[4]
				if entry != "U":	#is not a gap
					scaffs.append(scaff)

	for s in scaffs:
		if s not in tpfdict:
			probs.append(s)
	if len(probs) > 0:
		#print("\nagp and tpf not in sync, ",probs[0],"not in tpf for example\n")
		sys.exit()
		

def commas(number):

	n=list(str(number))
	for i in range(len(n))[::-3]:
			n.insert(i+1,",")
	if n[-1]==",":
		del n[-1]
	if n[0]==",":
		del n[0]
	
	return "".join(n)
	
	
def contig_lens(tpfdict):

	ctg_lengths={}
	for k,v in tpfdict.items():
		for i in v:
			if not "GAP" in i:
				x=i.strip().split()
				hi=int(x[1].split(":")[1].split("-")[1])
				lo=int(x[1].split(":")[1].split("-")[0])
				length=hi-lo
				append_dict(k,length,ctg_lengths)
	return ctg_lengths
	
def report_agp_discards(discards,agplines):
	print(f'HERE {discards}')
	for k,v in agplines.items():
		for line in v:
			i=line.strip().split()
			scaff=i[5]
			coord=i[6]
			if scaff in discards:
				##print(scaff,coord)
				##print(type(discards[scaff][0]),type(coord))
				if discards[scaff][0]==int(coord):
					print(line.strip())


def report_discreps(discards,agpdict,tpfchunks,outlinesfull):

	tpfcounts={}				#How many tpfchunks we've mapped to the agp (eg,scaffold_15 2)
	agpcounts={}
	results={}
	for k,v in tpfchunks.items():
		tpfbase=k.split("%")[0]
		if tpfbase not in tpfcounts:
			tpfcounts[tpfbase]=1
		else:
			tpfcounts[tpfbase]+=1
	for k,v in agpdict.items():
		for line in v:
			agpbase=line[0]
			if agpbase not in agpcounts:
				agpcounts[agpbase]=1
			else:
				agpcounts[agpbase]+=1

	for k,v in agpcounts.items():
		if v!= tpfcounts[k]:
			vals=[v,tpfcounts[k]]
			append_dict(k,vals,results)
			##print(k,v,tpfcounts[k])

	if len(results.keys())==1:
		print("\nDiscrepancy between agp and tpf - likely a break has been missed for some reason:\n")
	if len(results.keys())>1:
		print("\nDiscrepancies between agp and tpf - likely breaks have been missed for some reason:\n")
		
	for k,v in results.items():
		print(k+":\tagp count:\t"+str(v[0][0])+"\ttpf count:\t"+str(v[0][1]))


def report_stats(agpdict,outlinesfull,breakpoint,gsize,gsize2,texel,tpfout,tpf,discards,agplines,joins,named_haps,sexchrm):

	breaks=0

	scaffs=[]
	for k,v in agpdict.items():
		for line in v:
			scaff=line[0]
			if scaff not in scaffs:
				scaffs.append(scaff)	#First time we find a scaff, put it in scaffs
			else:
				breaks+=1		#Subsequent occurences must be breaks

	if len(discards)>0:
		print("\n"+100*"="+"\n")
		print("\n\n\tUnused AGP lines - either snap-mode artefacts we don't care about, or something\n\telse like a missing break.  Check these and if necessary fix these manually\n\n")
		report_agp_discards(discards,agplines)
		print("\n"+100*"="+"\n")	
	if len(sexchrm.keys())==1:
		print("Sex chromosome:\t\t\t"+"".join(sexchrm.keys()))
	elif len(sexchrm.keys())==2 and "W" in sexchrm.keys():
		print("Sex chromosomes:\t\tZW")
	elif len(sexchrm.keys())==2 and "Y" in sexchrm.keys():
		print("Sex chromosomes:\t\tXY")	
	else:
		print("No sex chromosomes defined")	

	if len(tpf) > len(tpfout):
		a=len(tpf)-len(tpfout)
		b=tpfout+(a*" ")
		print("length "+tpf+"\t\t"+str(commas(gsize))+" bp\t(input)")
		print("length "+b+"\t\t"+str(commas(gsize2))+" bp\t(output)")
	elif len(tpfout) > len(tpf):
		a=len(tpfout)-len(tpf)
		b=tpf+(a*" ")
		print("length "+b+"\t\t"+str(commas(gsize))+" bp\t(input)")
		print("length "+tpfout+"\t\t"+str(commas(gsize2))+" bp\t(output)")
	else:
		print("length "+tpf+"\t\t"+str(commas(gsize))+" bp\t(input)")
		print("length "+tpfout+"\t\t"+str(commas(gsize2))+" bp\t(output)")

	print("\nTexel length:\t\t\t"+commas(texel)+"bp\t\t(texels in map: 32,767)")
	#not sure if any difference between 2 break counting methods - check they're the same:
	if breaks==len(breakpoint):
		print("Break count:\t\t\t"+str(commas(breaks))) 
		#print("Break count:t\t"+str(commas(len(breakpoint)))) 
	else:
		print("Break count is either\t\t"+str(commas(breaks))+" or "+str(commas(len(breakpoint))))
	print("Join count:\t\t\t"+str(commas(joins)))	
	print("Haps count:\t\t\t"+str(commas(len(named_haps))))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  			
	chng=breaks+joins 
	changes=round(chng/(gsize/1000000000),1)
	##print("Interventions per Gb:\t\t"+str(round(changes/(gsize/1000000000),1)))
	if changes==0:
		print("No manual interventions!\n")
	elif changes >=1 and changes <= 49:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t(very low)\n")
	elif changes >=30 and changes <= 49:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t(below average)\n")
	elif changes >=50 and changes <= 199:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t(around average)\n")
	elif changes >=200 and changes <= 349:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t(above average)\n")
	elif changes >=350 and changes <= 399:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t(moderately high)\n")
	elif changes >=400 and changes <= 549:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t(high)\n")
	elif changes >=550 and changes <= 949:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t(very high)\n")
	elif changes >=950:
		print("Interventions per Gb:\t\t"+str(changes)+"\t\t(unusually high)\n")	


#We can use this file to see where rapid_pretext wants to break the genome based on the AGP.  This is useful should we need to switch to the XL version of the script		
def write_dividers(dividers):

	outfile="dividers.tsv"
	with open(outfile,'w') as fout:
		fout.write("#scaffold\tAGP\tTPF\n")
		for k,v in dividers.items():
			scaff=k.split(":")[0]
			agp=k.split(":")[1]
			tpf=str(v)
			new=[scaff,agp,tpf]
			a="\t".join(new)
			fout.write(a+"\n")


def tag_tpfchunks(tagdict,outlinesfull,tpfchunks):	
		
	tpfchunktags={}
	for k,v in tpfchunks.items():
		##print(k,v)
		for a,tags in tagdict.items():
			achunk=a.split("#")[1]
			if achunk == k:
				append_dict(k,tags[0],tpfchunktags)

	##print(tpfchunktags)

	#{'scaffold_38_ctg1%1': [['W']], 'scaffold_35_ctg1%1': [['W']], 'scaffold_60_ctg1%1': [['W']], 'scaffold_33_ctg1%1': [['W']], 'scaffold_49_ctg1%1': [['UNLOC', 'Z']], 'scaffold_34_ctg1%1': [['HAPLOTIG']]}
				
	return tpfchunktags

	
#def get_chunk_len(tpfchunks):

#	results={}
#	for k,v in tpfchunks.items():
#		length=0
#		for line in v:
#			#print(line)
		
	
#get haplines and removes them from output	
def get_haps(tpfchunktags,tpfchunks,output):

	scaff_dict={}
	hapcomponents=[]
	haptpfchunklens={}
	final={}
	
	##print(tpfchunktags)
	for k,v in tpfchunktags.items():
		##print(k,v)
		for tag in v[0]:
			##print(k,tag)
			if "HAPLOTIG" in tag:
				##print(k,tag)
				total=0
				for line in tpfchunks[k]:
					if not "GAP" in line:
						##print(line)
						lo=int(line.split(":")[1].split("-")[0].split("\t")[0])
						hi=int(line.split(":")[1].split("-")[1].split("\t")[0])
						length=(hi-lo)+1
						total+=length

				haptpfchunklens[k]=total
				
	for k in haptpfchunklens.keys():
		for line in tpfchunks[k]:
			if not "GAP" in line:
				comp=line.split()[1]
				hapcomponents.append(comp)

	prev=""
	chrm=""
	for line in output:
		if line.split()[1] not in hapcomponents:
			if not "GAP" in line:
				chrm=line.split()[2]
				##print(chrm,line)
				append_dict(chrm,line,scaff_dict)
				prev=chrm
			else:
				##print(line)
				if chrm !="":
					append_dict(prev,line,scaff_dict)
						
	prev=""
	for k,v in scaff_dict.items():
		for line in v:
			if not "GAP" in line:
				chrm=line.strip().split()[2]
				append_dict(k,line,final)
				prev=chrm
			elif "GAP" not in final[k][-1]:
				append_dict(k,line,final) 
				
	#for k,v in final.items():
	#	for line in v:
	#		#print(k,line)

	#scaff_dict is output minus haps
	#haptpfchunklens contains the hap tpfchunk key and its length so we can name haps by size
	return final, haptpfchunklens


def get_unlocs(tpfchunktags,tpfchunks,haplessoutput,haptpfchunklens, sex_chrms,tagged):

	scaff_dict={}
	unloctpfchunklens={}
	chrms_containing_unlocs=[]
	check={}
	
	for k,v in tpfchunktags.items():
		##print(k,v)
		total=0
		for tag in v[0]:
			if "UNLOC" in tag:
				for line in tpfchunks[k]:
					if not "GAP" in line:
						lo=int(line.split(":")[1].split("-")[0].split("\t")[0])
						hi=int(line.split(":")[1].split("-")[1].split("\t")[0])
						length=(hi-lo)+1
						total+=length
				unloctpfchunklens[k]=total

		
	named_unlocs = name_unlocs(unloctpfchunklens,tpfchunktags,sex_chrms,tagged)
	##print(named_unlocs)

	return named_unlocs

	
def name_haps(haptpfchunklens):

	sizes=[]
	named_haps={}
	iteration=0
	for k,size in haptpfchunklens.items():
		if size not in sizes:
			sizes.append(size)
	sortedsizes=sorted(sizes, reverse=True)
	for size in sortedsizes:
		for k,v in haptpfchunklens.items():
			if size==v:
				iteration+=1
				name="H_"+str(iteration)
				named_haps[k]=name


	return named_haps


def name_unlocs(unloctpfchunklens,tpfchunktags,sex_chrms,tagged):

	sex2={}
	sizes_={}	#chrm key to size list in descending order of unlocs
	sizes={}
	named_unlocs={}
	iteration=0
	for k,v in sex_chrms.items():
		for i in v:
			if i not in sex2:
				sex2[i]=k

	for k, v in tagged.items():
		#Scaffold_2:3#scaffold_2_ctg1%3
		aroot=k.split(":")[0]
		troot=k.split("#")[1]
		if troot in unloctpfchunklens:
			size=unloctpfchunklens[troot]
		##print(size)
			if aroot not in sizes_:
				sizes_[aroot]=[size]
			else:
				sizes_[aroot].append(size)

	for k, v in sizes_.items():
		sv=sorted(v)[::-1]
		sizes[k]=sv

	for k,v in sizes.items():
		##print(k,v)
		for a, b in unloctpfchunklens.items():
			##print(k,a,v)
			for i in v:
				if i==b:	#if it's the right size
					##print(k,sex2)
					if k not in sex2:
						##print(k,a,i,str(v.index(i)+1), k.replace("Scaffold","PRTXT")+"_unloc_"+str(v.index(i)+1))
						named_unlocs[a]=k.replace("Scaffold_",prefix)+"_unloc_"+str(v.index(i)+1)
					else:
						##print(k,a,i,str(v.index(i)+1), sex2[k]+"_unloc_"+str(v.index(i)+1))
						named_unlocs[a]=sex2[k]+"_unloc_"+str(v.index(i)+1)
				#else:
					##print(k,a,i,b)
	
	#{'scaffold_49_ctg1%1': 'Z_unloc_1', 'scaffold_54_ctg1%1': 'Z_unloc_2', 'scaffold_57_ctg1%1': 'Z_unloc_3'}		
	return named_unlocs


def prepare_haps_tpf(named_haps, tpfchunks):

	outlines={}

	for k,v in named_haps.items():
		for line in tpfchunks[k]:
			if not "GAP" in line:
				x=line.strip().split()
				x[2]=v
				newline="\t".join(x)
				append_dict(v,newline,outlines)
			else:
				append_dict(v,line.strip(),outlines)
			
	return outlines

			
def write_hap_tpf(hapoutlines,tpfout):

	haptpf="haps_"+tpfout
	with open(haptpf,'w') as fout:
		for k,v in hapoutlines.items():
			for line in v:
				fout.write(line+"\n")
				

#Ensuring that all members of a chr are labelled with the sex even if all members not labelled in AGP				
def sex_components(tpfchunks,tpfchunktags,tagged,sex_chrms):

	comp_sex={}

	for k,v in tagged.items():
		chrm=k.split(":")[0]
		tpc=k.split("#")[1]
		for a,b in sex_chrms.items():
			if chrm==b[0]:
				for line in tpfchunks[tpc]:
					if not "GAP" in line:
						comp=line.split()[1]
						comp_sex[comp]=a
								
	return comp_sex
	
	
def apply_sex(hapoutlines,comp_sex,sex_chrms):
	
	sex_prtxt={}
	for k,v in sex_chrms.items():
		s=prefix+v[0].split("_")[1]
		sex_prtxt[s]=k
	##print(sex_prtxt)
	
	outlines={}
	for k,v in hapoutlines.items():
		for line in v:
			x=line.strip().split()
			scaff=x[2]
			if scaff in sex_prtxt:
				x[2]=sex_prtxt[scaff]
				newline="\t".join(x)
				append_dict(k,newline,outlines)
			else:
				append_dict(k,line,outlines)

	return outlines


#all the components that belong to unloc scaffolds 	
def get_unloc_comps(named_unlocs,tpfchunks):

	results={}
	for k,v in named_unlocs.items():
		##print(tpfchunks[k],v)
		for line in tpfchunks[k]:
			comp=line.split()[1]
			append_dict(comp,v,results)
			
			
	return results


def update_output_unlocs(unloc_comps,output,sex_chrms):

	results={}
	check={}
	prob_comps={}

	for k,v in output.items():
		for line in v:
			##print(k,line)
			if not "GAP" in line:
				x=line.strip().split()
				comp=x[1]
		
				if comp in unloc_comps:
					##print(line,unloc_comps[comp])
					newline=line.replace(x[2],unloc_comps[comp][0])
					append_dict(k,newline,results)
				else:
					append_dict(k,line,results)
			else:
				append_dict(k,line,results)

	#Making a check on results to ensure no internal unlocs
	for k,v in results.items():
		for line in v:
			if not "GAP" in line:
				append_dict(k,line,check)

	#Remove legitimate terminal unlocs...
	for k,v in check.items():
		while "unloc" in v[0]:
			del v[0]
		while "unloc" in v[-1]:
			del v[-1]

	#...then look to see if any remain
	for k,v in check.items():
		##print(k,v)	
		#if "unloc" in v:
		for line in v:
			if "unloc" in line:
				##print(k,line)
				comp=line.split()[1]
				append_dict(k,comp,prob_comps)

	#if len(prob_comps)>0:
	for k,v in prob_comps.items():
		while len(v)>2:
			del v[1]

	##print(prob_comps)
	for k,v in output.items():
		if k in prob_comps:
			if len(prob_comps[k])>0:
				msg=[k,prob_comps[k]]
				##print(msg)
				append_dict("internal_unloc",msg,errors)

	return results


#remove original painted scaff key and replace keys with new keys derived from the new scaffold name (ie unloc keys)	
def update_chr_keys(outlines):

	updated={}
	chrm=""
	for k,v in outlines.items():
		for line in v:
			if not "GAP" in line:
				newkey=line.split()[2]
				chrm=newkey
				append_dict(newkey,line.strip(),updated)
			else:
				append_dict(chrm,line.strip(),updated)
				
	remove_excess_gaps(updated)

	return updated


def remove_excess_gaps(outlines):
	#Fixing up excess GAP lines
	
	#Terminal gap lines
	for k,v in outlines.items():	
		##print(k,v)
		#removing gaps from start/end of chromosomes
		#if "GAP" in v:
		while "GAP" in v[0]:
			del v[0]
		while "GAP" in v[-1]:
			del v[-1]
			
	#internal gap duplicates
	lastline=""
	for k,v in outlines.items():
		for e,line in enumerate(v):
			if line==lastline and "GAP" in line:
				##print("here", line)
				lastline=line
				del v[e]
			else:
				lastline=line

	return outlines


def parse_errors():

	if len(errors.keys())>0:
		#print("\n"+borderlen*"=")
		#print("\n\t\t\t\tERRORS\n")
		for k,v in errors.items():
			if k=="2nd_match":
				if len(set(v[1:]))==1:
					print(v[0])
				else:
					print("These occur only once in TPF - Either break missing or pretext editing error\n")
				for i in set(v[1:]):
					print(i+"  occurs"+str(v.count(i)+1)+"  times in AGP")
			if k=="hap_unloc":
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()
			if k=="multiple_sex":
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()
			if k=="multiple_sex2":
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()
			if k=="internal_unloc":
				#print("\nUnloc scaffs are internal to painted chromosomes.  Please move and rerun\n")
				for i in v:
					if len(i[1])==1:
						print(i[0],i[1][0])
					else:
						print(i[0],i[1][0]+" to "+i[1][1])
				#print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()	
			if k=="hetero_sex_haplo":
				#print("Doesn't make sense - haplotig painted into heterogametic sex chromosome:\n")
				for i in v:
					print(i[0])
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()
			if k=="xssex":
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")
				sys.exit()	
			if k=="badsex":
				for i in v:
					print(i)
				print("\n\t\t>>> PLEASE FIX AGP TAGS AND RERUN <<<\n")	
				sys.exit()				
			else:
				if len(v)>1:
					for i in v:
						print(i)
					print("\n\t\t>>> PLEASE FIX BREAKS FILE AND RERUN <<<\n")
					sys.exit()	###IMPORTANT - turn this back on once program is written	
				else:
					print(v[0])
					print("\n\t\t>>> PLEASE FIX BREAKS FILE AND RERUN <<<\n")
					sys.exit()	###IMPORTANT - turn this back on once program is written


def check_componentsH(finalout2,hapoutlines,checkin):


	finalout=[]
	for k,v in finalout2.items():
		for line in v:
			if not "GAP" in line:
				x=line.split()
				comp=x[1]
				finalout.append(comp)
	for k,v in hapoutlines.items():
		for line in v:
			if not "GAP" in line:
				x=line.split()
				comp=x[1]
				finalout.append(comp)		

	for c in checkin:
		if c in finalout:
			finalout.remove(c)

	if len(finalout) !=0:
		#print("\n\\t>>>   PROBLEM - components missing from out files!!   <<<n")
		for f in finalout:
			print(f)


def main():

	parser = argparse.ArgumentParser(description='Designed to take pretext generated AGP and fit your assembly TPF to it.') 

	#positional args
	parser.add_argument('tpf', metavar='tpf', type=str, help='assembly TPF with gaps as needed to allow rearrangement to match the edited PretextView map.')
	parser.add_argument('agp', metavar='agp', type=str, help='Pretext agp')
	#parser.add_argument('breaks', metavar='breaks', type=str, help='breaks file')
	#parser.add_argument('fasta', metavar='fasta', type=str, help='original assembly fasta')

	#display help when misusage
	if len(sys.argv) <2: 
		parser.print_help()

	args = parser.parse_args()  #gets the arguments
	start_time = datetime.now()

	#print("\n")

	#print("\nChecking "+args.tpf+" sanity...\n\n")
	tpf_sanity(args.tpf)


	tpfdict, errors2=parse_tpf(args.tpf)
	
	compare_scaff(tpfdict,args.agp)

	ctg_lengths=contig_lens(tpfdict)
	gsize,texel=genome_size(tpfdict)
	scafflens=lens(tpfdict)
	#for k,v in scafflens.items():
		##print(k,v)
	fragcutoff=1*texel
	agpdict,discards,agplines, tagdict, sex_chrms = scaffs_from_agp(args.agp,fragcutoff,scafflens,ctg_lengths,texel)	#All the agp order and orientation information
	#print(agpdict)
	#for k,v in agplines.items():
	#	#print(k,v)
	
	##print(tagdict)

	report_errors(errors2)
	checkin=components_from_dict(tpfdict)

	dividers=agp_dividers(agpdict)
	
	#for k,v in dividers.items():
		##print(k,v)
	
	breakpoint=nearest(tpfdict,dividers,fragcutoff,discards,scafflens)
	#print(breakpoint)

	write_dividers(breakpoint)	#Produce output which we can parse to create input for the XL versin of the script (if a curator runs this version of the script instead of the XL version by mistake).

	tpfchunks=breaktpf(tpfdict,breakpoint)

	outlines,joins, tagged = outputlist(tpfchunks,agpdict,tagdict) #joins is join count

	#tagged eg:
	#{'Scaffold_1:-#scaffold_141%1': [['Z', 'UNLOC']], 'Scaffold_1:-#scaffold_32%1': [['Z']], 'Scaffold_9:-#scaffold_68%1': [['HAPLOTIG']], 'Scaffold_9:-#scaffold_81%1': [['UNLOC']], 'Scaffold_9:+#scaffold_8%1': [['UNLOC']], 'Scaffold_12:+#scaffold_15%1': [['W']]}
	##print(tagged)
	
	#for l in outlines:
	#	#print(l)

	checkout=check_components(outlines)
	outlinesfull=reinstate_lines(tpfdict,outlines,checkin,checkout)


	tpfchunktags = tag_tpfchunks(tagged,outlinesfull,tpfchunks)
	##print(tpfchunktags)

	haplessoutput, haptpfchunklens = get_haps(tpfchunktags,tpfchunks,outlines)

	named_unlocs = get_unlocs(tpfchunktags,tpfchunks,haplessoutput,haptpfchunklens, sex_chrms,tagged)
	##print(named_unlocs)
	unloc_comps = get_unloc_comps(named_unlocs,tpfchunks)
	unlochaplessoutput = update_output_unlocs(unloc_comps, haplessoutput, sex_chrms)

	named_haps = name_haps(haptpfchunklens)
	hapoutlines = prepare_haps_tpf(named_haps, tpfchunks)
	
	comp_sex = sex_components(tpfchunks,tpfchunktags,tagged, sex_chrms)
	
	sexedunlochaplessoutput = apply_sex(unlochaplessoutput,comp_sex, sex_chrms)
	
	finalout1 = update_chr_keys(sexedunlochaplessoutput)
	
	finalout2 = remove_excess_gaps(finalout1)
	
	##print(errors)
	parse_errors()
	
	gsize2=genome_size2(outlinesfull)
		
	tpfout="rapid_prtxt.tpf"

	write_output_tpf(finalout2,tpfout)
	#if len(hapoutlines)>0:
	write_hap_tpf(hapoutlines,tpfout)
	
	check_componentsH(finalout2,hapoutlines,checkin)	#Belt and braces checking again that final naming based on tags hasn't lost any components	
	
	report_stats(agpdict,outlinesfull,breakpoint,gsize,gsize2,texel,tpfout,args.tpf,discards,agplines,joins,named_haps,sex_chrms)
	
	report_discreps(discards,agpdict,tpfchunks,outlinesfull)

	#print("\nChecking "+tpfout+" sanity...\n")
	
	tpf_sanity(tpfout)
	
	#print("Written:\n")
	#print(location(tpfout))
	#print(location("haps_"+tpfout))
	end_time=datetime.now()
	#print('\n\nFINISHED:\t{}'.format(end_time - start_time)+"\n\n")


if __name__ == '__main__':
	main()
