# coding: utf-8

import sys
import os 
import re
import pandas as pd


###########
## Parameters:
##########
#HMM pwd
HMM = ''
#Pipeline peptidase pwd
PP = ''


# Proteome fasta 
proteomefile=sys.argv[1]


##########
## 1st step : process hmmsearch with all hmm peptidase profil
#########

try:
	os.stat(PP+'temp')
except:
	os.mkdir(PP+'temp')


os.chdir(HMM)

def RunHMMpep(PEP):
	for f in os.listdir(PP+'HMMMotif/'+PEP):
		os.system('hmmsearch '+PP+'HMMMotif/'+PEP+'/'+f+' '+proteomefile+' > '+PP+'temp/'+f.strip('hmmprofil.txt')+'hmmresult.txt')

RunHMMpep('PreP')
RunHMMpep('OOP')
RunHMMpep('MPP')
RunHMMpep('SPP')

##########
## 2st step : hmm result to 1 csv 
#########

def OutHMMstoOneCSV(hmmfilelist,outcsv,filtrCov='on',threshold='------ inclusion threshold ------'):

        DF=pd.DataFrame(columns=['Index','Evalue','score','bias','EvalueD','scoreD','biasD','exp','N','Find in'])
        p = re.compile('^\s+(?P<Evalue>[.\S]+)\s+(?P<score>[.\S]+)\s+(?P<bias>[.\S]+)'+
                                        '\s+(?P<EvalueD>[.\S]+)\s+(?P<scoreD>[.\S]+)\s+(?P<biasD>[.\S]+)'+
                                        '\s+(?P<exp>[.\S]+)\s+(?P<N>\d)\s+(?P<Sequence>[.\S]+)'+
                                        '\s+(?P<Description>.+)')
        SUPR=[]

        for hmmfile in hmmfilelist:
                f_CP=open(PP+'temp/'+hmmfile,'r')
                l=f_CP.readline()
                while threshold not in l and 'Internal pipeline statistics summary:' not in l:
                        m=re.search(p,l)
                        if m is not None:
                                i=len(DF)
                                DF.at[i,'Index']=m.group('Sequence')#.split('|')[1]
                                DF.at[i,'Name']=m.group('Sequence')#.split('|')[2]
                                #DF.at[i,'TaxID']=m.group('Description').split(' ')[0][3:]
                                DF.at[i,'Evalue']=m.group('Evalue')
                                DF.at[i,'score']=m.group('score')
                                DF.at[i,'bias']=m.group('bias')
                                DF.at[i,'EvalueD']=m.group('EvalueD')
                                DF.at[i,'scoreD']=m.group('scoreD')
                                DF.at[i,'biasD']=m.group('biasD')
                                DF.at[i,'exp']=m.group('exp')
                                DF.at[i,'N']=m.group('N')
                                DF.at[i,'Find in']=hmmfile.split('/')[-1]
                        l=f_CP.readline()
                f_CP.close()
                DF,SUPR=Add_Couv(DF,hmmfile,filtrCov,SUPR)
        DF=DF.drop(SUPR)
        DF.to_csv(outcsv,sep='\t',index=False)


def Add_Couv(DF,hmmfile,filtrCov,SUPR):
        p = re.compile('^\s+(?P<N>[.\S]+)\s+(?P<sign>[.\S]+)\s+(?P<score>[.\S]+)\s+(?P<bias>[.\S]+)\s+(?P<cEval>[.\S]+)\s+(?P<iEval>[.\S]+)\s+(?P<hmmfrom>[.\S]+)\s+(?P<hmmto>[.\S]+)\s+(?P<ali>[.\S]+)\s+(?P<alifrom>[.\S]+)\s+(?P<alito>[.\S]+)')
        f_CP=open(PP+'temp/'+hmmfile,'r')
        l=f_CP.readline()
        while l!='':
                if l[0:6]=='Query:':
                        LenHMM=l.split('[')[-1][2:-2]
                if l[0:2]=='>>':
                        for ind in DF[DF['Index']==l.split(' ')[1]].index :
                                if hmmfile.split('/')[-1] in DF.at[ind,'Find in']:
                                        l=f_CP.readline()
                                        l=f_CP.readline()
                                        l=f_CP.readline()
                                        Sign=[]#!:significatif ?:non significatif
                                        Cov=[]#Cal coverage on hmm
                                        iEval=[]#indep eval
                                        From=[]#pos deb hits sur seq
                                        To=[]#pos fin hits sur seq
                                        while l != '\n':
                                                m=re.search(p,l)
                                                if '#' not in l and '---' not in l and m is not None: 
                                                        Sign.append(m.group('sign'))
                                                        Cov.append( (float(m.group('hmmto')) - float(m.group('hmmfrom')) + 1)/float(LenHMM))
                                                        iEval.append(m.group('iEval'))
                                                        From.append(m.group('alifrom'))
                                                        To.append(m.group('alito'))

                                                l=f_CP.readline()
                                        DF.at[ind,'Dom_hits']=int(len(Sign))
                                        for n in range(len(Sign)):
                                                DF.at[ind,'Signif_hit'+str(n+1)]=Sign[n]
                                                DF.at[ind,'Coverage_hit'+str(n+1)]=Cov[n]
                                                DF.at[ind,'iEval_hit'+str(n+1)]=iEval[n]
                                                DF.at[ind,'aliFrom_hit'+str(n+1)]=From[n]
                                                DF.at[ind,'aliTo_hit'+str(n+1)]=To[n]

                                        # 70% covery
                                        if filtrCov=='on' and len([x for x in Cov if x >0.7])==0:
                                        	# DF.drop([ind] , inplace=True)
                                        	SUPR.append(ind)




                l=f_CP.readline()
        return(DF,SUPR)

OutHMMstoOneCSV(os.listdir(PP+'temp'),PP+'out.csv','on')

os.system('rm -f * '+PP+'temp/*')


##########
## 3rd step : Fasta candidate sequence
#########
def HMMcsv_to_fasta(Nom_fichier_fasta,Nom_fichier_HMM):
	DF=pd.read_csv(Nom_fichier_HMM,sep='\t')
	LIndex=set(DF[DF.columns[0]])
	NOK=0
	f=open(Nom_fichier_fasta,'w')
	fREF=open(proteomefile,'r')
	l=fREF.readline()
	while NOK<len(LIndex) and l!='':

		if l[0]=='>' and (l.split(' ')[0][1:] in LIndex or l[1:-1].strip(' ') in LIndex):
			NOK+=1
			f.write(l)
			l=fREF.readline()
			while(len(l)>0 and l[0]!='>') :
				f.write(l)
				l=fREF.readline()
		else:
			l=fREF.readline()
	fREF.close()
	f.close()

HMMcsv_to_fasta(PP+'out.fasta',PP+'out.csv')

##########
## 5rd step : motifs hmmsearch 
#########
os.chdir(HMM)

if  os.stat(PP+'out.fasta').st_size != 0:
	for f in os.listdir(PP+'HMMMotif/Pfamm'):
		os.system('hmmsearch '+PP+'HMMMotif/Pfamm/'+f+' '+PP+'out.fasta > '+PP+'temp/'+f.strip('hmmprofil.txt')+'_hmmresult.txt')


##########
## 6rd step : Add hmmresult to csv
#########
if  os.stat(PP+'out.fasta').st_size != 0:
	OutHMMstoOneCSV(os.listdir(PP+'temp'),PP+'outbis.csv','off')
os.system('rm -f * '+PP+'temp/*')

if  os.stat(PP+'out.fasta').st_size != 0:
	DF1=pd.read_csv(PP+'out.csv',sep='\t')
	DF2=pd.read_csv(PP+'outbis.csv',sep='\t')
	DF=pd.concat([DF1,DF2],ignore_index=True)
	DF.to_csv(PP+'Final_'+proteomefile.split('/')[-1]+'_out.csv',sep='\t',index=False)
else:
	DF=pd.read_csv(PP+'out.csv',sep='\t')
	DF.to_csv(PP+'Final_'+proteomefile.split('/')[-1]+'_out.csv',sep='\t',index=False)


##########
## 7rd step : Create final result file
#########

def ResumeResultMotifFiltre(DF,out_file):
	if out_file not in os.listdir(PP): 
		RESDATA=pd.DataFrame(columns=['ProteomeFile','Species','TaxID','#PreP','#SPP','#OOP','#MPP','listPreP','listSPP','listOOP','listMPP'])
	else:
		RESDATA=pd.read_csv(PP+out_file,sep='\t')
	IND=len(RESDATA)
	RESDATA.at[IND,'ProteomeFile']=proteomefile.split('/')[-1]
	PEP=['PreP','OOP','SPP','MPP']
	LPEP=[[],[],[],[]]
	for prot in set(DF['Index']):
		df=DF[DF['Index']==prot]
		motifs=[]
		pep=[]
		for i in df.index:
			FI=df.at[i,'Find in'].split('_')
			if FI[0] == 'Motif':
				for j in range(1,int(df.at[i,'Dom_hits'])+1):
					if df.at[i,'Signif_hit'+str(j)]=='!':
						motifs.append(FI[1])
			elif FI[1] not in pep:
				pep.append(FI[1])

		if 'PreP' in pep and 'M16' in motifs and 'M16C' in motifs and 'M16Cassoc' in motifs:
				LPEP[0].append(df.at[i,'Index'])
		elif 'OOP' in pep and 'M3' in motifs :
				LPEP[1].append(df.at[i,'Index'])
		elif 'SPP' in pep and 'M16' in motifs and len([x for x in motifs if x=='M16C']) > 1 :
				LPEP[2].append(df.at[i,'Index'])
		elif ('MPPA' in pep or 'MPPB' in pep )and 'M16' in motifs and 'M16C' in motifs :
				LPEP[3].append(df.at[i,'Index'])

	for pep in PEP:
		RESDATA.at[IND,'#'+pep]=len(LPEP[PEP.index(pep)])
		RESDATA.at[IND,'list'+pep]=set(LPEP[PEP.index(pep)])
	RESDATA.to_csv(PP+out_file,sep='\t',index=False)


ResumeResultMotifFiltre(DF,'ResumeFinding.csv')


os.system('rm -f * '+PP+'out*')

