import numpy as np
import pandas as pd
import regex
import re
from Bio import SeqIO
import Bio

###CHANGE PAM AND WINDOW INFO HERE###
PAM = 'NGG'
windowstart = 4
windowend = 8
###CHANGE PAM AND WINDOW INFO HERE###

def RC(seq):
	encoder = {'A':'T','T':'A','C':'G','G':'C','N':'N','R':'Y','Y':'R', 'M':'K', 'K':'M', 'S':'S', 'W':'W', 'H':'D', 'B':'V', 'V':'B', 'D':'H'}
	rc = ''
	for n in reversed(seq):
		rc += encoder[n]
	return rc

def create_PAM(pam):
	encoder = {'A':'A','T':'T','G':'G','C':'C','R':'[A|G]','Y':'[C|T]','N':'[A|T|C|G]','M':'[A|C]','K':'[G|T]','S':'[C|G]','W':'[A|T]','H':'[A|C|T]','B':'[C|G|T]','V':'[A|C|G]','D':'[A|G|T]'}
	enc_pam = {'f':'','r':''}
	rc_pam = RC(pam)
	for n,m in zip(pam, rc_pam):
		enc_pam['f'] += encoder[n]
		enc_pam['r'] += encoder[m]
	return enc_pam
enc_pam = create_PAM(PAM)

windowlen=windowend-windowstart+1
lenpam=len(PAM)


####################################

CFTR=pd.read_csv('CFTR2_com.csv', encoding = "ISO-8859-1")
CFTR_size = len(CFTR.index)

print CFTR_size

SNV = CFTR[CFTR['Variant cDNA name'].str.contains(">")]
SNV_RS = SNV[SNV['rsID'].str.contains('rs')]
SNV_RS.to_csv('RSlist.csv', sep='\t', encoding = 'utf-8')
SNV_RS['rsID'] = SNV_RS['rsID'].str.replace('rs',"")



SNV_RS.to_csv('SNV_strip.csv', sep='\t', encoding = 'utf-8')
BE_CtoT_editable = SNV_RS[SNV_RS['Variant cDNA name'].str.contains("T>C")]
BE_GtoA_editable = SNV_RS[SNV_RS['Variant cDNA name'].str.contains("A>G")]
ABE_TtoC_editable = SNV_RS[SNV_RS['Variant cDNA name'].str.contains("C>T")]
ABE_AtoG_editable = SNV_RS[SNV_RS['Variant cDNA name'].str.contains("G>A")]



handle = open('mart_export.txt' , "rU")
Flanks={}
for record in SeqIO.parse(handle, "fasta") :
	Flanks[int(record.id.split("|")[3].strip('rs'))]=regex.findall('.{25}[^A,T,C,G].{25}', record.seq.tostring())
handle.close()
F=pd.DataFrame({'rsID': list(Flanks.keys()), 'Flanks': [x for x in Flanks.values()]})

F['rsID'] =F['rsID'].astype(str)
SNV_RS['rsID'] =SNV_RS['rsID'].astype(str)
ABE_AtoG_editable = F.merge(ABE_AtoG_editable, left_on='rsID', right_on='rsID', how='right')
ABE_TtoC_editable = F.merge(ABE_TtoC_editable, left_on='rsID', right_on='rsID', how='right')
BE_CtoT_editable = F.merge(BE_CtoT_editable, left_on='rsID', right_on='rsID', how='right')
BE_GtoA_editable = F.merge(BE_GtoA_editable, left_on='rsID', right_on='rsID', how='right')


ABE_AtoG_editable['gRNAs']=None
ABE_AtoG_editable['gRNAall']=None
for i in range(len(ABE_AtoG_editable)):
    if type(ABE_AtoG_editable.iloc[i].Flanks)==list and ABE_AtoG_editable.iloc[i].Flanks!=[]:
        test=ABE_AtoG_editable.iloc[i].Flanks[0]
        # define a potential gRNA spacer for each window positioning
        gRNAoptions=[test[(26-windowstart-j):(26-windowstart-j+lenpam+20)] for j in range(windowlen)]
        #if there is an appropriate PAM placed for a given gRNA spacer
        #save tuple of gRNA spacer, and the position of off-target As in the window
        gRNA=[(gRNAoptions[k],[x.start()+1 for x in re.finditer('A',gRNAoptions[k]) if windowstart-1<x.start()+1<windowend+1]) for k in range(len(gRNAoptions)) if regex.match(enc_pam['f'], gRNAoptions[k][-lenpam:])]
        gRNAsingleA=[]
        for g,c in gRNA:
            #if the target A is the only A in the window save this as a single C site
            if g[windowstart-1:windowend].count('A')==0:
                gRNAsingleA.append(g)
            #OPTIONAL uncomment the ELIF statement if you are interest in filtered based upon position of off-target C
            #if the target C is expected to be editted more efficiently than the off-target Cs, also save as a single C Site
            #elif all([p<priority[x] for x in c]):
                #gRNAsingleC.append(g)
        
        ABE_AtoG_editable.gRNAs.iloc[i]=gRNAsingleA
        ABE_AtoG_editable.gRNAall.iloc[i]=[g for g,c in gRNA]


ABE_TtoC_editable['gRNAs']=None
ABE_TtoC_editable['gRNAall']=None
for i in range(len(ABE_TtoC_editable)):
    if type(ABE_TtoC_editable.iloc[i].Flanks)==list and ABE_TtoC_editable.iloc[i].Flanks!=[]:
        test=ABE_TtoC_editable.iloc[i].Flanks[0]
        gRNAoptions=[test[(25+windowstart+j-20-lenpam):(25+windowstart+j)] for j in range(windowlen)]
        gRNA=[(gRNAoptions[k],[20+lenpam-x.start() for x in re.finditer('T',gRNAoptions[k]) if windowstart-1<20+lenpam-x.start()<windowend+1]) for k in range(len(gRNAoptions)) if regex.match(enc_pam['r'], gRNAoptions[k][:lenpam])]
        gRNAsingleT=[]
        for g,c in gRNA:
            if g[20+lenpam-windowstart-windowlen+1:20+lenpam-windowstart+1].count('T')==0:
                gRNAsingleT.append(g)
            #elif all([p<priority[x] for x in c]):
                #gRNAsingleC.append(g)
        ABE_TtoC_editable.gRNAs.iloc[i]=gRNAsingleT
        ABE_TtoC_editable.gRNAall.iloc[i]=[g for g,c in gRNA]


BE_CtoT_editable['gRNAs']=None
BE_CtoT_editable['gRNAall']=None
for i in range(len(BE_CtoT_editable)):
    if type(BE_CtoT_editable.iloc[i].Flanks)==list and BE_CtoT_editable.iloc[i].Flanks!=[]:
        test=BE_CtoT_editable.iloc[i].Flanks[0]
        # define a potential gRNA spacer for each window positioning
        gRNAoptions=[test[(26-windowstart-j):(26-windowstart-j+lenpam+20)] for j in range(windowlen)]
        #if there is an appropriate PAM placed for a given gRNA spacer
        #save tuple of gRNA spacer, and the position of off-target Cs in the window
        gRNA=[(gRNAoptions[k],[x.start()+1 for x in re.finditer('C',gRNAoptions[k]) if windowstart-1<x.start()+1<windowend+1]) for k in range(len(gRNAoptions)) if regex.match(enc_pam['f'], gRNAoptions[k][-lenpam:])]
        gRNAsingleC=[]
        for g,c in gRNA:
            #if the target C is the only C in the window save this as a single C site
            if g[windowstart-1:windowend].count('C')==0:
                gRNAsingleC.append(g)
            #OPTIONAL uncomment the ELIF statement if you are interest in filtered based upon position of off-target C
            #if the target C is expected to be editted more efficiently than the off-target Cs, also save as a single C Site
            #elif all([p<priority[x] for x in c]):
                #gRNAsingleC.append(g)
        
        BE_CtoT_editable.gRNAs.iloc[i]=gRNAsingleC
        BE_CtoT_editable.gRNAall.iloc[i]=[g for g,c in gRNA]


BE_GtoA_editable['gRNAs']=None
BE_GtoA_editable['gRNAall']=None
for i in range(len(BE_GtoA_editable)):
    if type(BE_GtoA_editable.iloc[i].Flanks)==list and BE_GtoA_editable.iloc[i].Flanks!=[]:
        test=BE_GtoA_editable.iloc[i].Flanks[0]
        gRNAoptions=[test[(25+windowstart+j-20-lenpam):(25+windowstart+j)] for j in range(windowlen)]
        gRNA=[(gRNAoptions[k],[20+lenpam-x.start() for x in re.finditer('G',gRNAoptions[k]) if windowstart-1<20+lenpam-x.start()<windowend+1]) for k in range(len(gRNAoptions)) if regex.match(enc_pam['r'], gRNAoptions[k][:lenpam])]
        gRNAsingleG=[]
        for g,c in gRNA:
            if g[20+lenpam-windowstart-windowlen+1:20+lenpam-windowstart+1].count('G')==0:
                gRNAsingleG.append(g)
            #elif all([p<priority[x] for x in c]):
                #gRNAsingleC.append(g)
        BE_GtoA_editable.gRNAs.iloc[i]=gRNAsingleG
        BE_GtoA_editable.gRNAall.iloc[i]=[g for g,c in gRNA]



ABE_AtoG_hasPAM=ABE_AtoG_editable[[type(x)==list and x!=[] for x in ABE_AtoG_editable.gRNAall]]
ABE_AtoG_hasPAM.to_csv('ABE_AtoG_hasPAM.csv', sep='\t', encoding = 'utf-8')
ABE_AtoG_SingleA=ABE_AtoG_editable[[type(x)==list and x!=[] for x in ABE_AtoG_editable.gRNAs]]
ABE_AtoG_SingleA.to_csv('ABE_AtoG_SingleA.csv', sep='\t', encoding = 'utf-8')

ABE_TtoC_hasPAM=ABE_TtoC_editable[[type(x)==list and x!=[] for x in ABE_TtoC_editable.gRNAall]]
ABE_TtoC_hasPAM.to_csv('ABE_TtoC_hasPAM.csv', sep='\t', encoding = 'utf-8')
ABE_TtoC_SingleA=ABE_TtoC_editable[[type(x)==list and x!=[] for x in ABE_TtoC_editable.gRNAs]]
ABE_TtoC_SingleA.to_csv('ABE_TtoC_SingleA.csv', sep='\t', encoding = 'utf-8')

BE_CtoT_hasPAM=BE_CtoT_editable[[type(x)==list and x!=[] for x in BE_CtoT_editable.gRNAall]]
BE_CtoT_hasPAM.to_csv('BE_CtoT_hasPAM.csv')
BE_CtoT_SingleC=BE_CtoT_editable[[type(x)==list and x!=[] for x in BE_CtoT_editable.gRNAs]]
BE_CtoT_SingleC.to_csv('BE_CtoT_SingleC.csv')

BE_GtoA_hasPAM=BE_GtoA_editable[[type(x)==list and x!=[] for x in BE_GtoA_editable.gRNAall]]
BE_GtoA_hasPAM.to_csv('BE_GtoA_hasPAM.csv')
BE_GtoA_SingleG=BE_GtoA_editable[[type(x)==list and x!=[] for x in BE_GtoA_editable.gRNAs]]
BE_GtoA_SingleG.to_csv('BE_GtoA_SingleG.csv')

