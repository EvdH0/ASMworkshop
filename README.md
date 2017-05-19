## ASM Microbe '17 workshop
_Room 365_

Instructors:
- [Morten Sommer](https://twitter.com/moasommer)
- Lejla Imamovic
- [Eric van der Helm](https://twitter.com/EricvdHelm)
- Mostafa Ellabaan

# Timetable
- 12:30 Morten Sommer:  Introduction to antibiotic resistome
- 12:45 Lejla Imamovic:  Resistome profiling  (from the lab to databases, overview of the data and analysis to be used during the workshop)
- 13:00 Mostafa Ellabaan: Code examples part 1
- 13:30 break
- 13:45 Code examples part 2
- 14:45 break
- 15:30 Eric van der Helm:  Rapid resistome profiling using Nanopore sequencing
- 15:45 Morten Sommer:  Antibiotic resistance gene exchange networks

# Tasks
- Task 1.      Introduce basic command line and tools how to get the data on/off server
      Output 1:   Ability to use basic commands needed for the workshop (data, database locations etc)

- Task 2.        Assign ORFs in GeneMark for selected dataset ((MinION last run)
      Output 2:    Protein (pro) and nucleotide (nt) sequences in fasta format
      Question:    How many ORFs are predicted?
- Task 3.        Format databases (e.g. CARD nt and pro or Resfinder nucletide)
      Output 3:    Database for further use, knowledge on importance of databases
                           Also mention  that one can curate their own database  
                           Task 4A.     Annotate insert nt sequences with Resfinder online server (promote DTU)
- Task 4B.      Annotate ORFs nt sequences against Resfinder database (blastn)
      Output 4:    Manually curated horizontally transferred list of genes from inserts/contigs  and ORFs. 
      Questions:  How many inserts have antibiotic resistance gene?
                            Is there a difference between DTU web server (contig) (4A) and annotated ORFs (4B)?
                            - Task 5.       Annotate protein sequences against CARD protein database
      Output 5:   List of antibiotic resistance genes from ORFs. 
      Questions: How many inserts have antibiotic resistance gene? 
                           Input from participant of interesting hits?

- Task 6.       Annotate protein sequences using hmm model in Pfam 
      Output 6:   List of genes from ORFs. 
      Questions: What can we see from Pfam domain list?
If there is time left, annotate ORFs against optimal databases using plasmid database (to look for mobility). 
We can use megablast to map the nt ORFs to plasmid database.
# Content
- [1. Login on amazon cloud](#foo)
- [2. Genemark](#foo)
- [3. Resfinder](#foo)
- [4. CARD](#card)


# Workshop

## 1. Login in to the amazon cloud

## Prepare your directories
```shell
mkdir ASM cd ASM mkdir Databases Programs Data Scripts
ls
echo "export ASM=$(pwd)" >> ~/.bashrc source ~/.bashrc
tail -n 1 ~/.bashrc
ls $ASM
```


## 2. Finding Open Reading Frames (ORFs)
### Installation of Genemark
- Go to http://exon.gatech.edu/Genemark/ and search for programs and click on the link it will take you to http://exon.gatech.edu/Genemark/license_download.cgi 
- Select GeneMarkS v.4.30 and LINUX 64 and fill the form then click on the button "I agree to the terms of this license agreement" 
- Download both the program and the key

```shell
cd $ASM/Programs
mkdir GenMarkS
cd GenMarkS
## copy link address of the download program here
## and wget it 
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_ZH3mh/gm_key_64.gz
## copy link address of the key using 64_bit version
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_ZH3mh/genemark_suite_linux_64.tar.gz

## extract both the program and the key. 
tar -xvf genemark_suite_linux_64.tar.gz 
gzip -d gm_key_64.gz


## copy genmarks key to the home directory
cp gm_key_64 ~/.gm_key

## 

cd genemark_suite_linux_64/gmsuite
ls 
## you should observe the 

## In the scope of this work shop we will use this command gmhmmp mainly.

## Create an Alias to simplify the job

echo "export gmsuite=$(pwd)" >> ~/.bashrc
echo "alias getORFUsingGeneMarks=$(pwd)/gmhmmp" >> ~/.bashrc

## to see whether the alias is there or not. 
##tail -n 5 ~/.bashrc

## to make the effect run. 
source ~/.bashrc
```

### Running genemark
```shell
## to see options available just write this 
cd $ASM/Data

## put the orginal file on the web.

wget $address

getORFUsingGeneMarks -m $gmsuite/MetaGeneMark_v1.mod \
-A orfs.protein.fa -D orfs.nucleotide.fa -o all.orfs.result gud_np.fa

cat orfs.nucleotide.fa 
cat orfs.protein.fa 

## barcodeing the genes for both protein and nucleotide.

## coding genes


## to know how many orfs are there

grep ">"  orfs.protein.fa  | wc -l 

## or 

grep ">"  orfs.nucleotide.fa | wc -l


## OPTIONAL 

cat orfs.nucleotide.fa | sed 's/%//g' | 
awk -F"\n" '{if( index($1,">") > 0) {printf "\n"$1"\t"} else {printf $1}}' | 
awk -F" " '{print substr($1,2,length($1))"\t"$0}' | awk '{print $0}' | 
awk -F"\t" 'BEGIN {print "ORFID\tGENMarksORFID\tDescription\tsequence"}
{if(FNR>1) print "ORFID"FNR-1"\t"$0}' > orfs.nucleotide.tab

cat  orfs.protein.fa | sed 's/%//g' | 
awk -F"\n" '{if( index($1,">") > 0) {printf "\n"$1"\t"} else {printf $1}}' | 
awk -F" " '{print substr($1,2,length($1))"\t"$0}' | awk '{print $0}' | 
awk -F"\t" ' BEGIN {print "ORFID\tGENMarksORFID\tDescription\tsequence"}
{if(FNR>1) print "ORFID"FNR-1"\t"$0}' > orfs.protein.tab
```


## 3. The first database, Resfinder

### Installation of Resfinder
```shell
## to download the most updated version of resfinder
mkdir $ASM/Databases/Resfinder
cd $ASM/Databases/Resfinder
git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
cd resfinder_db


## make sure all files have a linux format

ls *fsa | while read file; do dos2unix $file; done 


## compile them all in one file.

cat *.fsa > ../All.resfinder.fsa



##ls *fsa | while read file; do cat $file | awk -F"\t" '{print } END {print "\n"}'; done  > ../All.resfinder.fsa

## OPTIONAL
 ls *fsa | while read file ; do f=$(echo "$file" | cut -f1 -d.); grep ">" $file | 
 sed 's/>//g' | awk -F"\t" -v class="$f" '{print $1"\t"class}'; done | 
 awk -F"\t" 'BEGIN {print "Gene\tClass"} {print $1"\t"$2 }'  > ../Resfinder.gene.class

cd ..


## OPTIONAL
## translate the NA resfinder to protein
$EmbossFolder/transeq -sequence All.resfinder.fsa -outseq All.resfinder.AA.fa -frame 1

head All.resfinder.prot.fa
sed 's/_[1-6]$//g' All.resfinder.AA.fa  | head 

sed 's/_[1-6]$//g' All.resfinder.AA.fa -i

mkdir blastNA  ## blastAA

## OPTIONAL
## compiling blast Protein database for Resfinder
cd blastAA
$blastExecFolder/makeblastdb -in ../All.resfinder.AA.fa -title RESFINDERProt \
-out RESFINDERProt -input_type fasta   -hash_index -dbtype prot

echo "export RESFINDERProt=$(pwd)/RESFINDERProt" >> ~/.bashrc
tail -n 1 ~/.bashrc
source ~/.bashrc

cd ..
## compiling blast Nucleotides database for Resfinder


## install binary blast 

cd blastNA

$blastExecFolder/makeblastdb -in ../All.resfinder.fsa -title RESFINDERNucl \
-out RESFINDERNucl -input_type fasta   -hash_index -dbtype nucl

## REMOVE IT
echo "export RESFINDERNucl=$(pwd)/RESFINDERNucl" >> ~/.bashrc
tail -n 1 ~/.bashrc
source ~/.bashrc
```

### Install the BLAST tool
To download blast visit page https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

```shell
cd $ASM/Programs
mkdir blast
cd blast


## download latest version of blast
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz
## download md5 signature to confirm that you download the blast completely and nothing is missing
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz.md5


## check the signature if it is work
md5sum --check ncbi-blast-2.6.0+-src.tar.gz.md5 
## you should see "ncbi-blast-2.6.0+-src.tar.gz: OK"
```

```shell
cd  ncbi-blast-2.6.0+-src/c++
./configure --prefix=$(pwd)
time make
time make install
```
### Run Resfinder
```shell
cd $ASM/Data
mkdir Resfinder
cd Resfinder
## nucleotide versus Resfinder nucleotide

## replace  RESFINDERNucl with right path.
time $blastExecFolder/blastn -query ../orfs.nucleotide.fa -db $RESFINDERNucl -outfmt 6 \
 -max_target_seqs 1 -evalue 1E-50 -word_size 6 -num_threads 1 -out orf.resfinder.NA.versus.NA.tab 
 
 
## protein versus Resfinder nucleotide
$blastExecFolder/tblastn -query ../orfs.protein.fa -db $RESFINDERNucl -outfmt 6 \
 -max_target_seqs 10 -evalue 1E-50 -word_size 3 -num_threads 28 -out orf.resfinder.AA.versus.NA.tab 
 
## protein versus Resfinder nucleotide
$blastExecFolder/tblastx -query ../orfs.nucleotide.fa -db $RESFINDERNucl -outfmt 6 \
 -max_target_seqs 10 -evalue 1E-50 -word_size 3 -num_threads 28 -out orf.resfinder.trans.NA.versus.NA.tab 
 
## protein versus Resfinder protein
$blastExecFolder/blastp -query  ../orfs.protein.fa -db $RESFINDERProt -outfmt 6 \
 -max_target_seqs 10 -evalue 1E-50 -word_size 6 -num_threads 28  -out orf.resfinder.AA.versus.AA.tab 
 
## nucleotide versus Resfinder protein
$blastExecFolder/blastx -query ../orfs.nucleotide.fa -db $RESFINDERProt -outfmt 6 \
 -max_target_seqs 10 -evalue 1E-50 -word_size 6 -num_threads 28  -out orf.resfinder.NA.versus.AA.tab
 ```
 
 ### Next?
```shell
 ## fruther annotation
## maping gene to class

## 
Resfinder.gene.class 
JoinTwoFilesBasedOnKeys.sh 2 1 $ASM/Databases/Resfinder/Resfinder.gene.class orf.resfinder.NA.versus.NA.tab 

## how many gene classes are there?
## how many betalactamases are there?
## how many inserts have resistance gene?
```


## 4. CARD

### Install the CARD database
```shell
# getting and installing
cd $ASM/Databases
mkdir CARD
cd CARD

mkdir RawData
cd RawData

wget https://card.mcmaster.ca/download/0/broadstreet-v1.1.8.tar.gz

tar -xvf broadstreet-v1.1.8.tar.gz
```

```shell
cat nucleotide_fasta_*pro*homo* > ../All.CARD.NA.fa
cat protein_fasta_* > ../All.CARD.AA.fa
mkdir blastNA blastAA


## compiling blast protein database for CARD
cd blastAA

$blastExecFolder/makeblastdb -in ../All.CARD.AA.fa -title CARDProt \
-out CARDProt -input_type fasta   -hash_index -dbtype prot

echo "export CARDProt=$(pwd)/CARDProt" >> ~/.bashrc
tail -n 1 ~/.bashrc
source ~/.bashrc

cd ..

## compiling blast Nucleotides database for CARD
cd blastNA

$blastExecFolder/makeblastdb -in ../All.CARD.NA.fa -title CARDNucl \
-out CARDNucl -input_type fasta   -hash_index -dbtype nucl

echo "export CARDNucl=$(pwd)/CARDNucl" >> ~/.bashrc
tail -n 1 ~/.bashrc
source ~/.bashrc
```

### Run the CARD database
```shell
## Run CARD databses

cd $ASM/Data
mkdir CARD
cd CARD
## nucleotide versus card nucleotide
$blastExecFolder/blastn -query ../orfs.nucleotide.fa -db $CARDNuc -outfmt 6 \
 -max_target_seqs 10 -evalue 1E-50 -word_size 6 -num_threads 28 -out orf.card.NA.versus.NA.tab 
 
 
## protein versus card nucleotide
$blastExecFolder/tblastn -query ../orfs.protein.fa -db $CARDNuc -outfmt 6 \
 -max_target_seqs 10 -evalue 1E-50 -word_size 6 -num_threads 28 -out orf.card.AA.versus.NA.tab 
 
 
## protein versus card protein
$blastExecFolder/blastp -query ../orfs.protein.fa -db $CARDProt -outfmt 6 \
 -max_target_seqs 10 -evalue 1E-50 -word_size 6 -num_threads 28  -out orf.card.AA.versus.AA.tab 
 
## nucleotide versus card protein
## translaste sequence
## $Emboss/transeq -in  ../orfs.nucleotide.fa -frame [1,6] ../orfs.protein.t1.fa

$blastExecFolder/blastp -query  ../orfs.protein.t1.fa -db $CARDProt -outfmt 6 \
 -max_target_seqs 10 -evalue 1E-50 -word_size 6 -num_threads 28  -out orf.card.AA.versus.AA.tab
 ```

## Combine data
TBD
