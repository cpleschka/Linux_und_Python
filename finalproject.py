print('Experiment:\nRNA-seq of lung and ileum samples at 1 and 3 days post infection (dpi) from chickens infected with either low pathogenic (H5N2) or highly pathogenic (H5N1) avian influenza')

print('\nVariables\nday1 or day3\nlung or ileum\nH5N1 or H5N2 or none')

print('\nI checked for day 3 ileum infected with H5N1 against the control and day 3 ileum infected with the H5N2 strain against the control\n')

#import packages
import numpy as np 
import csv

#open file with experimental design and find Accession numbers that are relavant for my chosen groups
exp_file=open('experiment-design', 'r')
read_exp=csv.reader(exp_file, delimiter='\t')
next(read_exp)
counter=0
run_ids=dict()
infect=dict()
organ=dict()
time=dict()
for run in read_exp:
    counter+=1
    run_ids[counter]=run[0] #jedem run eine id zuteilen, die ein int ist
    infect[counter]=run[1]
    organ[counter]=run[5]
    time[counter]=run[7]
    

exp_file.close() 

#save all accession numbers for each group that I am looking at
none_d3_ileum_runs=[]
for id in run_ids:
    if infect[id]=='none' and organ[id]=='ileum' and time[id]=='3 day':
        none_d3_ileum_runs.append(run_ids[id])
#print(none_d3_ileum_runs)

H5N1_d3_ileum_runs=[]
for id in run_ids:
    if infect[id]=='H5N1' and organ[id]=='ileum' and time[id]=='3 day':
        H5N1_d3_ileum_runs.append(run_ids[id])
#print(H5N1_d3_ileum_runs)

H5N2_d3_ileum_runs=[]
for id in run_ids:
    if infect[id]=='H5N2' and organ[id]=='ileum' and time[id]=='3 day':
        H5N2_d3_ileum_runs.append(run_ids[id])
#print(H5N2_d3_ileum_runs)

#open file with raw counts and extract counts by group, calculate averages for multiples
counts_file=open('raw-counts', 'r')
read_counts=csv.reader(counts_file, delimiter='\t')
header=next(read_counts)
none_d3_ileum_cols=[]
H5N1_d3_ileum_cols=[]
H5N2_d3_ileum_cols=[]
for c in range(len(header)):
    if header[c] in none_d3_ileum_runs: none_d3_ileum_cols.append(c)
    if header[c] in H5N1_d3_ileum_runs: H5N1_d3_ileum_cols.append(c)
    if header[c] in H5N2_d3_ileum_runs: H5N2_d3_ileum_cols.append(c)


my_gene_ids=dict()
none_d3_ileum_avgs=dict()
H5N1_d3_ileum_avgs=dict()
H5N2_d3_ileum_avgs=dict()
counter=0
for gene in read_counts:
    counter+=1
    my_gene_ids[counter]=gene[0]
    gene_none_d3_ileum_avg=0
    for i in none_d3_ileum_cols:
        tmp=int(gene[i])
        gene_none_d3_ileum_avg+=tmp
    gene_none_d3_ileum_avg=gene_none_d3_ileum_avg/3
    none_d3_ileum_avgs[counter]=gene_none_d3_ileum_avg
    gene_H5N1_d3_ileum_avg=0
    for i in H5N1_d3_ileum_cols:
        tmp=int(gene[i])
        gene_H5N1_d3_ileum_avg+=tmp
    gene_H5N1_d3_ileum_avg=gene_H5N1_d3_ileum_avg/3
    H5N1_d3_ileum_avgs[counter]=gene_H5N1_d3_ileum_avg
    gene_H5N2_d3_ileum_avg=0
    for i in H5N2_d3_ileum_cols:
        tmp=int(gene[i])
        gene_H5N2_d3_ileum_avg+=tmp
    gene_H5N2_d3_ileum_avg=gene_H5N2_d3_ileum_avg/3
    H5N2_d3_ileum_avgs[counter]=gene_H5N2_d3_ileum_avg
    
#str type array containing the values and gene ids, by group
expr_arr=np.zeros((len(my_gene_ids)+1,4), dtype='U800')
for i in range(0,len(my_gene_ids)+1):
    if i==0:
        expr_arr[i,0]='gene id'
        expr_arr[i,1]='none d3 ileum avg'
        expr_arr[i,2]='H5N1_d3_ileum_avg'
        expr_arr[i,3]='H5N2_d3_ileum_avg'
    else:
        expr_arr[i,0]=my_gene_ids[i]
        expr_arr[i,1]=str(none_d3_ileum_avgs[i]) 
        expr_arr[i,2]=str(H5N1_d3_ileum_avgs[i])
        expr_arr[i,3]=str(H5N2_d3_ileum_avgs[i])

#print(expr_arr[0:11,:])

#calculate log2foldchanges and save up/ down regulated and unchanged genes into lists

up_genes_H5N1vsnone=[]
down_genes_H5N1vsnone=[]
no_genes_H5N1vsnone=[]
up_genes_H5N2vsnone=[]
down_genes_H5N2vsnone=[]
no_genes_H5N2vsnone=[]
#cal. diff. in expression for each gene
log2foldchange=np.zeros((len(my_gene_ids)+1,3), dtype='U800')
for i in range(1,len(my_gene_ids)+1):
    log2foldchange[i,0]=float(i)
    foldchange=float(expr_arr[i,2])-float(expr_arr[i,1]) #the gene ids are wrong here
    if foldchange>1.0: up_genes_H5N1vsnone.append(my_gene_ids[i])
    elif -1.0<=foldchange<=1.0: no_genes_H5N1vsnone.append(my_gene_ids[i])
    elif foldchange<-1.0: down_genes_H5N1vsnone.append(my_gene_ids[i])
    else: print('oh no')
    foldchange_abs=abs(foldchange)
    if foldchange!=0:log2foldchange[i,1]=np.log2(foldchange_abs) #is log2 correct?
    else: log2foldchange[i,1]=0
    
for i in range(1,len(my_gene_ids)+1):
    log2foldchange[i,0]=float(i)
    foldchange=float(expr_arr[i,3])-float(expr_arr[i,1]) #the gene ids are wrong here
    if foldchange>1.0: up_genes_H5N2vsnone.append(my_gene_ids[i])
    elif -1.0<=foldchange<=1.0: no_genes_H5N2vsnone.append(my_gene_ids[i])
    elif foldchange<-1.0: down_genes_H5N2vsnone.append(my_gene_ids[i])
    else: print('oh no')
    foldchange_abs=abs(foldchange)
    if foldchange!=0:log2foldchange[i,2]=np.log2(foldchange_abs) #is log2 correct?
    else: log2foldchange[i,2]=0
log2foldchange=np.delete(log2foldchange,0,0)
#print(log2foldchange[0:11,:])

#evaluate differential genes and print results
Patt1_count=0
for k in up_genes_H5N1vsnone:
    if k in up_genes_H5N2vsnone:
        Patt1_count+=1
print('Expression pattern 1: Up in H5N1 and H5N2, '+str(Patt1_count)+' Genes')

Patt2_count=0
for k in up_genes_H5N1vsnone:
    if k in down_genes_H5N2vsnone:
        Patt2_count+=1
print('Expression pattern 2: Up in H5N1 and down in H5N2, '+str(Patt2_count)+' Genes')

Patt3_count=0
for k in up_genes_H5N1vsnone:
    if k in no_genes_H5N2vsnone:
        Patt2_count+=1
print('Expression pattern 3: Up in H5N1 and not changed in H5N2, '+str(Patt3_count)+' Genes')
