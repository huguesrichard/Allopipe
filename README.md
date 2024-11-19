# AlloPipe

The AlloPipe tool is a computational workflow which imputes<br/>
&nbsp;&nbsp;&nbsp;&nbsp;(i) **directional amino acid mismatches** and <br/>
&nbsp;&nbsp;&nbsp;&nbsp;(ii) the related **minor histocompatibility antigens** candidates <br/>
within a pair of annotated human genomic datasets.

<br/>
The product is provided free of charge, and, therefore, on an "as is" basis, without warranty of any kind.
<br/>
AlloPipe is also available as a web application: https://www.allogenomics.com

<br/>

## In a nutshell

**The AlloPipe tool is divided into two sequential modules: Allo-Count, then Allo-Affinity**
<br/>
![AP](https://github.com/user-attachments/assets/ce944565-18ff-4e98-98bb-5d522e46cb21)

<br/>

&nbsp;&nbsp;&nbsp;&nbsp; **(1) Allo-Count imputes the directional amino acid mismatches from two genomic datasets**

Allo-Count reformats the relevant data from the variant-annotated .VCF file(s), performs a stringent data cleaning and computes the **directional comparison** of the genomic sequences. Allo-Count returns: <br/>
- **a quantitative output** called **Allogenomic Mismatch Score (AMS)** which is a discrete quantitative variable numbering the directional amino acid mismatches
- **a qualitative output** stored in the mismatch table, providing information about the non-synomous SNP contributing to the AMS
  
<br/>  

> **Directional comparison**
> 
> The sample comparison is directional: and accounts for either polymorphisms that are present in the donor but absent in the recipient (*donor-to-recipient*) or that are present in the recipient but absent in the donor (*recipient-to-donor*).\
> **_Donor-to-recipient_** accounts for polymorphisms present by the donor but absent by the recipient, i.e. triggerring the recipient's immune system after **solid organ transplantation**.\
> **_Recipient-to-donor_** accounts for polymorphisms present by the recipient but absent by the donor, i.e. triggerring the donor's immune system after **allogeneic haematopoietic cell transplantation**.

<br/>


&nbsp;&nbsp;&nbsp;&nbsp; **(2) Allo-Affinity imputes the candidates minor histocompatibility antigens**<br/>

Allo-Affinity **reconstructs peptides** of requested length around the polymorphisms, then returns their affinity towards HLA molecules using [NetMHCpan softwares](https://pubmed.ncbi.nlm.nih.gov/32406916/). Allo-Affinity returns: <br/>
- **a quantitative output** called **affinity-AMS (af-AMS)** which is a discrete quantitative variable numbering the candidates minor histocompatibility antigens
- **a qualitative output** stored in the af-AMS table, providing information about the peptides contributing to the af-AMS

*4-digits HLA typing has to be provided by the user for the HLA molecules of interest.*\
*Be careful with the terms of use of [NetMHCpan softwares](https://pubmed.ncbi.nlm.nih.gov/32406916/)*


<br/>
<br/>

**There are two modes of operation for each module: as single pair or as multiple pairs**
	
- **Single pair**\
Run as 'single pair mode' if you aim to compute AMS and/or af-AMF for one pair at a time. \
You need to provide one variant-annotated .VCF file per individual.
 
- **Multiple pairs**  
Run as 'multiple pairs mode' if you aim to compute AMS and/or af-AMF for more than one pair at a time.\
You need to provide one unique variant-annotated .VCF file containing the genotype of all individuals you want to analyse - i.e. a merged .VCF file - and the [example.csv](./data/example.csv) of the pairs you want to process.

---

# Table of contents

1. [Before getting started](#before)
	1. [Requirements](#requirements)
	2. [AlloPipe installation](#install)
	3. [Preprocess the input data: variant annotation](#vep)
    
<br/>

2. [Run the AlloPipe workflow](#run)
	1. [Launch Allo-Count](#ams_run)
		1. [Single pair](#single_ams)
		2. [Multiple pairs](#multi_ams)
		3. [Exploring the AMS table](#ams_table)
	2. [Launch Allo-Affinity](#aams_run)
		1. [Single pair](#single_aams)
		2. [Multiple pairs](#multi_aams)
		2. [Exploring the af-AMS table](#aams_table)

<br/>
     
  3. [Tutorial](#tuto)
---

## Before getting started <a name="before"></a>

### Requirements<a name="requirements"></a>

AlloPipe installation specifically requires
1. [Python](https://www.python.org/downloads/) ≥ 3.6 (developed on 3.9)

2. [Conda](https://docs.anaconda.com/free/working-with-conda/) installed in the suitable version for your operating system and python version, as we recommend to install the [dependencies](https://github.com/huguesrichard/Allopipe/blob/main/requirements.txt) in an dedicated environment.
  
3. [NetMHCpan](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1) and [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) downloaded as command line tools.\
*Make sure you use NetMHCpan in accordance with their user licence.*
   
<br/>

### AlloPipe installation <a name="install"></a>

1) Clone the repository from git\
*You might be requested to create a token for you to log in. See the [GitHub tutorial](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-personal-access-token-classic)*

2) Create a conda environment dedicated to the AlloPipe workflow
3) Install the AlloPipe requirements within the conda environment
   
The following command lines will perform the steps 1-2-3:

		git clone https://github.com/huguesrichard/Allopipe.git
		cd Allopipe
 		conda create --name Allopipe python=3.9
  		conda activate Allopipe
		python -m pip install -r requirements.txt

<br/>

### Preprocess the input data: variant annotation <a name="vep"></a>

<br/>

**AlloPipe input file(s) must be variant-annotated .VCF file(s). We highly recommend to perform the variant annotation with the most recent version of VEP using the command line installation and all the arguments specified below.**

<br/>

>*Any variant annotator could be used at this step, but keep in mind that AlloPipe has been developed with .VCF v4.2 files annotated with VEP command line installation for versions 103, 104, 110 and 111.*

<br/>

To install the VEP command line tool, follow the installation tutorial available [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#download).\
		During the installation, you will be asked if you want to download **cache** files, **FASTA** files and **plugins**.
   - We **recommend to download the cache files** for the assembly of your VCF files to be able to run VEP offline.\
	Download the VEP cache files which correspond to your Ensembl VEP installation and genome reference!
   - We **recommend to download the FASTA files** for the assembly of your VCF files to be able to run VEP offline.\
	Download the FASTA files which correspond to your Ensembl VEP installation and genome reference!
   - We **do not recommend to download any plugin**
      
 We then recommend to **add VEP to your PATH** by adding the following line to your ```~/.profile``` or ```~/.bash_profile```:

		export PATH=%%PATH/TO/VEP%%:${PATH}

<br/>

Run the following command to annotate you VCF file(s) with VEP.\
**All specified options are mandatory, with the exception of the assembly if you only downloaded one cache file.**  
	
		vep --fork 4 --cache --assembly <GRChXX> --offline --af_gnomade -i <PATH-TO-FILE-TO-ANNOTATE/FILE-TO-ANNOTATE>.vcf -o <PATH-TO-ANNOTATED-FILE/ANNOTATED-FILE>.vcf --coding_only --pick_allele --use_given_ref --vcf 

Where:\
```<GRChXX>``` is the version of the genome used to align the sequences.\
```<PATH-TO-FILE-TO-ANNOTATE/FILE-TO-ANNOTATE>.vcf``` is the path to your file to annotate.\
```<PATH-TO-ANNOTATED-FILE/ANNOTATED-FILE>``` is the path to the directory and the name of the ouput annotated file.\

This command line works for individual .VCF files or multi-VCF files, whether compressed (.gvcf) or not (.vcf). 
Run this command for every file you want to input in AlloPipe.

**Once the variant-annotation of your file(s) is(are) complete, you are now ready to launch your first AlloPipe run!**

---

## Run the AlloPipe workflow <a name="run"></a>


### Launch Allo-Count  <a name="ams_run"></a>

**What does Allo-Count perform?**

From variant-annotated .VCF file(s), data are first reformated to obtain one dataframe per individual.\
Those dataframe are then filtered considering a set of quality metrics (defaults values):
- minimal depth per position (20x)
- maxiaml depth per position (400x)
- minimal allelic depth (5x)
- homozygosity threshold (0.8)
- GnomADe allele frequency threshold (0.01)
- genotype quality threshold (20)
- maximal indels length (3)

The curated dataframes are then queried to assess the **directional mismatches** between samples.

<br/>

> **Directional comparison**
> 
> The sample comparison is **directional** and accounts for either polymorphisms that are present in the donor but absent in the recipient (*donor-to-recipient*) or that are present in the recipient but absent in the donor (*recipient-to-donor*).\
> **_Donor-to-recipient_** accounts for polymorphisms present by the donor but absent by the recipient, i.e. triggerring the recipient's immune system after **solid organ transplantation**.\
> **_Recipient-to-donor_** accounts for polymorphisms present by the recipient but absent by the donor, i.e. triggerring the donor's immune system after **allogeneic haematopoietic cell transplantation**.
> 
<br/>


#### Launch Allo-Count for a single pair <a name="simple_ams"></a>

Once the variant-annotation is complete, go to the root of the AlloPipe directory to run the following commands in the terminal \
*Do not forget to activate your conda environment!* 

		cd src/
		python ams_pipeline.py -n <NAME-RUN> -p <NAME-PAIR> <PATH-TO-DONOR-ANNOTATED-FILE/ANNOTATED-FILE>.vcf <PATH-TO-RECIPIENT-ANNOTATED-FILE/ANNOTATED-FILE>.vcf <DIRECTION OF THE MISMATCH>


Where :\
```<NAME-RUN>``` is the name of the run\
```<NAME-PAIR>``` is the name of the pair\
```<PATH-TO-DONOR-ANNOTATED-FILE/ANNOTATED-FILE>.vcf``` is the path to the donor's annotated VCF \
```<PATH-TO-RECIPIENT-ANNOTATED-FILE/ANNOTATED-FILE>.vcf``` is the path to the recipient's annotated VCF \
```<DIRECTION OF THE MISMATCH>``` = dr, present in the donor but absent in the recipient ; rd, present in the recipient but absent in the donor.

<br/>
A complete helper function is provided

		python ams_pipeline.py --help
  

<br/>



#### Launch Allo-Count for multiple pairs <a name="multi_ams"></a>

It is possible to launch Allo-Count from an annotated .VCF merged file containing the genomic data of interest.\
In that case, you need to upload a [example.csv](./data/example.csv) specifying the donor/recipient pairs.

		cd src/
		python multiprocess_ams.py -n <NAME-RUN> <PATH-TO-THE-MERGED-ANNOTATED-FILE>.vcf <PATH-TO-THE-PAIR-LIST>.csv <DIRECTION OF THE MISMATCH>

Where:\
```<NAME-RUN>``` is the name of the run\
```<PATH-TO-THE-MERGED-ANNOTATED-FILE>.vcf``` is the path to the annotated merged VCF file\
```<PATH-TO-THE-PAIR-LIST>.csv``` is the path to the list pairing the sample (template provided in the tutorial)\
```<DIRECTION OF THE MISMATCH>``` is the direction of the mismatch as previously described

*Only one directional comparison is accepted within the same command line.*

We provide a complete helper function

		python multiprocess_ams.py --help


<br/>

> **Normalisation**
> 
> To avoid artefacts related to the quality of the sequencing that might lead to AMS lower or higher than expected, we provide to the user the ref/commun ratio.
> 

<br/>


### Exploring the AMS table <a name="ams_table"></a>

After the run is complete, have look at the **output/runs/<NAME-RUN>/** directory that was created.  
The directory is structured as followed :  
1. the **AMS/** directory contains a subdirectory created for these run parameters specifically, the AMS value contained in a .csv file.  
2. the **plots/** subdirectory
3. the **run_tables/** subdirectory contains all the tables created during the run. 

In the **run_tables/** directory, you can find the mismatches table that will give you direct information on the mismatched positions.  
In this table, you can find the following information :  

1. **VCF information**  
	1. **CHROM (str)**: Chromosome of the variant
	2. **POS (int)**: Position on the chromosome
	3. **ID_{x, y} (str)**: Reference SNP cluster ID for the donor (x) or recipient (y)
	4. **REF, ALT (str)**: REF and ALT alleles at the given position
	5. **QUAL_{x, y} (float)**: Phred-scaled quality score for the assertion made in ALT
	6. **FILTER_{x, y} (str)**: PASS if this position has passed all filters
	7. **FORMAT_{x, y} (list)**: Format of the sample column post AlloPipe processing
	8. **Sample_{x, y} (str)**: Sample information regarding the position. Note that the column name is the one provided in the original VCF
        - In the case of transplantation, Sample_x is the donor and Sample_y is the recipient
 

2. **Sample information**
	1. **GT_{x, y} (str)**: Predicted genotype of the sample
	2. **GQ_{x, y} (float)**: Score of quality of the predicted genotype
	3. **AD_{x, y} (str)**: Allelic depth 
	4. **FT_{x, y} (str)**: Sample genotype filter indicating if this genotype was “called”
	5. **phased_{x, y} (str)**: Predicted genotype containing phased information (if provided in the sample column)
	6. **DP_{x, y} (int)**: Sequencing Depth at position
	7. **TYPE_{x, y} (str)**: type of genotype (homozygous, heterozygous)
3. **VEP information**
	1. **consequences_{x, y} (int)**: All the columns with a consequence with the number of times it is recorded in transcripts for the variant
	2. **transcripts_{x, y} (str)**: Transcripts recorded for the variant
	3. **genes_{x, y} (str)**: Genes recorded for the variant
	4. **aa_REF, aa_ALT (str)**: Amino-acid for REF and ALT alleles for the variant
	5. **gnomADe_AF_{x, y} (float)**: Frequency of existing variant in gnomAD exomes combined population
	6. **aa_ref_indiv_{x, y}, aa_alt_indiv_{x, y} (str)**: REF and ALT amino-acids recorded for the sample (x and y)
	7. **aa_indiv_{x, y} (str)**: REF and ALT amino-acids combined in one column
4. **AlloPipe information**
	1. **diff (str)**: difference between the amino-acids of both samples
	2. **mismatch (int)**: number of mismatches in the diff field
	3. **mismatch_type (str)**: type of mismatch (homozygous, heterozygous)

<br/>

### Launch Allo-Affinity <a name="aams_run"></a>

<br/>

**What does Allo-Affinity perform?**

<br/>

From previously generated files that are the TABLE-MISMATCH and the TRANSCRIPT-TABLE, **Allo-Affinity reconstructs the set of peptides that are different between the donor and the recipient** for a given length defined by the user using the principle of the sliding window.

**The directionality of the mismatch is kept**, meaning that if Allo-Count has been run within the *donor-to-recipient* direction, only peptides exhibiting a polymorphism present by the donor but absent from the recipient will be reconstructed.\
In the same way, if Allo-Count has been run within the *recipient-to-donor direction*, only peptides exhibiting a polymorphism present by the recipient but absent from the donor will be reconstructed.

In the mean of assessing the peptides potentially trigerring the alloreactivity, the latter need to be embedded in HLA class I and/or class II molecules, therefore being of 9-AA long (9-mers) or 15-AA long (15-mers), even though any length is tolerated by Allo-Affinity.
Allo-Affinity eventually prepares the .fasta files containing the peptides and the appended files that are required to run by
 - [NetMHCpan4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
 - [NETMHCIIpan4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) 
 - tbu


<br/>

> To run NetMHCpan, you will need the following files:\
> **BEWARE! The number of the release has to match your VEP version!**\
> Do not forget to select the reference genome used to perfom the alignment.
> - <REFERENCE-GENOME>.cdna.all.fa.gz [https://ftp.ensembl.org/pub/release-XXX/fasta/homo_sapiens/](https://ftp.ensembl.org/pub/release-XXX/fasta/homo_sapiens/)
> - <REFERENCE-GENOME>.pep.fa.gz [https://ftp.ensembl.org/pub/release-XXX/fasta/homo_sapiens/](https://ftp.ensembl.org/pub/release-XXX/fasta/homo_sapiens/)
> - <REFERENCE-GENOME>.refseq.tsv.gz [https://ftp.ensembl.org/pub/release-XXX/tsv/homo_sapiens/](https://ftp.ensembl.org/pub/release-XXX/tsv/homo_sapiens/) \
> 
>*Replace XXX by your VEP version in the previous links*
> We provide the latest versions of those files [here](./data/Ensembl)

<br/>

Each of these tool imputes the affinity of the reconstructed peptides towards the HLA peptide grooves, therefore output **potential minor histocompatibility antigens (mHAgs)**
**Please note that the HLA typing has to be known before running the command line**, as the AlloPipe tool does not impute the HLA typing from genomic data.

<br/>

#### Simple pair <a name="simple_aams"></a>

Once the AMS run is complete, go back to the AlloPipe root directory and run this second set of commands:  


		cd src
		gzip -d <PATH-TO-REFERENCE-GENOME>.cdna.all.fa.gz
		gzip -d <PATH-TO-REFERENCE-GENOME>.pep.fa.gz
		gzip -d <PATH-TO-REFERENCE-GENOME>.refseq.tsv.gz
		python aams_pipeline.py 
		-M <PATH-TO-MISMATCH-TABLE>.tsv 
		-T <PATH-TO-TRANSCRIPT-TABLE>.tsv
		-E <PATH-TO-REFERENCE-GENOME>.cdna.all.fa
		-P <PATH-TO-REFERENCE-GENOME>.pep.fa
		-R <PATH-TO-REFERENCE-GENOME>.refseq.tsv
		-n <TEST-RUN> -p <TEST-PAIR> -l <PEP-LENGTH> --el_rank <THRESHOLD-FOR-EL> 
		-a <HLA-TYPING> 

Where:\
	```<PATH-TO-REFERENCE-GENOME>.cdna.all.fa.gz``` is the path to the gzip .cdna.fasta file\
	```<PATH-TO-REFERENCE-GENOME>.pep.fa.gz``` is the path to the gzip .pep.fa file\
	```<PATH-TO-REFERENCE-GENOME>.refseq.tsv.gz``` is the path to the gzip .refseq.tsv file\
 These three first lines should be run only the first time you launch Allo-Affinity
	```<PATH-TO-MISMATCH-TABLE>.tsv``` is the path to the mismatch table generated by Allo-Count\
	```<PATH-TO-TRANSCRIPT-TABLE>.tsv``` is the path to the transcript table generated by Allo-Count\
	```<NAME-RUN>``` is the name of the run. **It has to be consistent with the ```<NAME-RUN>``` used in the Allo-Count part** \
	```<NAME-PAIR>``` is the name of the pair \
 	```<PEP-LENGTH>``` is the length of peptided to be imputed \
	```<HLA-TYPING>``` is the HLA typing e.g. ```HLA-A*01:01,HLA-A*02:01,HLA-B*08:01,HLA-B*27:05,HLA-C*01:02,HLA-C*07:01```


#### Multiple pairs <a name="multi_aams"></a>

To be implemeted


### Getting your affinity-AMS (af-AMS) <a name="aams_results"></a>

This second step of AlloPipe uses the AMS information of the first step.  
You will find 3 new subdirectories in the **test_run/** directory :  
1.	the **AAMS/** directory contains a subdirectory created for these run parameters specifically, the AAMS value contained in a csv file.
2.	the **netMHCpan_out/** subdirectory contains all tables generated during the netMHCpan step.
3.	the **aams_run_tables/** subdirectory contains all the other tables created during the run


The AAMS value obtained with VEP v107 and netMHCpan4.1 is 34.

### Exploring the af-AMS table <a name="aams_mismatches"></a>

If you want more in-depth information on the mismatches contributing to the AAMS, you will find a mismatches table in the **aams_run_tables/** directory.  
It contains the mismatches information from the AMS run along with information provided by netMHCpan :
1. **NetMHCpan information**
	1. **hla_peptides (str)**: Potential ligand peptide built from VEP information and Ensembl information
	2. **Gene_id (str)**: Ensembl Gene ID
	3. **NB (int)**: Number of Weak Binding/Strong Binding peptides accross given HLA
	4. **EL-score (float)**: Raw prediction score
	5. **EL_Rank (float)**: Rank of the predicted EL-score compared to a set of random natural peptides
	6. **BA-score (float)**: Binding-Affinity score
	7. **BA_Rank (float)**: Rank of the predicted BA-score
	8. **HLA (str)**: Specified MHC molecule / Allele name
	9. **Transcript_id (str)**: Ensembl Transcript ID
	10. **Peptide_id (str)**: Ensembl Peptide ID


You can now get started with your files, check the [documentation](#docs/documentation.pdf) if you want more control over the filters that we implemented.  

<br/>

## Tutorial <a name="tuto"></a>

We provide a couple of example data in /tutorial, i.e. tutorial/donor_to_annotate.vcf and tutorial/recipient_to_annotate.vcf *(those files correspond to human chr6)*.\
To test your VEP installation, run the following command at the root of the AlloPipe directory :  
	
		vep --fork 4 --cache --assembly GRCh38 --offline --af_gnomade -i tutorial/donor_to_annotate.vcf -o tutorial/donor_annotated_VEP.vcf --vcf
		vep --fork 4 --cache --assembly GRCh38 --offline --af_gnomade -i tutorial/recipient_to_annotate.vcf -o tutorial/recipient_annotated_VEP.vcf --vcf 


Once the VEP annotation is complete, go to the root of the AlloPipe directory to run the following commands in the terminal :  

	cd src/
	python ams_pipeline.py -n test_run -p test_pair ../tutorial/donor_annotated_VEP.vcf ../tutorial/recipient_annotated_VEP.vcf rd

 If your AMS returns 44, congrats ! You successfully generated your first Allogenomic Mismatch Score (AMS) and related tables !

 Finally, to get your af-AMS and related table, run:
 
  	gzip -d ../data/Ensembl/Homo_sapiens.GRCh38.cdna.all.103.fa.gz
 	gzip -d ../data/Ensembl/Homo_sapiens.GRCh38.pep.all.103.fa.gz
	gzip -d ../data/Ensembl/Homo_sapiens.GRCh38.103.refseq.tsv.gz
	
	python aams_pipeline.py \
	-M ../output/runs/test_run/run_tables/test_pair_test_run_mismatches_20_400_5_gq_20_0.8_bl_3.tsv \
	-T ../output/runs/test_run/run_tables/test_pair_test_run_transcripts_pair_codons_20_400_5_gq_20_0.8_bl_3.tsv \
	-E ../data/Ensembl/Homo_sapiens.GRCh38.cdna.all.103.fa \
	-P ../data/Ensembl/Homo_sapiens.GRCh38.pep.all.103.fa \
	-R ../data/Ensembl/Homo_sapiens.GRCh38.103.refseq.tsv \
	-n test_run -p test_pair \
	-l 9 --el_rank 2 \
	-a HLA-A*01:01,HLA-A*02:01,HLA-B*08:01,HLA-B*27:05,HLA-C*01:02,HLA-C*07:01

Ir your af-AMS returns 33, you are all set !

You can now enjoy AlloPipe. We will be happy of any feedback !
