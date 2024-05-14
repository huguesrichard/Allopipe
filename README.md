# AlloPipe

The AlloPipe tool is a computational workflow imputing directional amino acid mismatches then related minor histocompatibility antigens candidate from human genomic datasets.

--- 
## In a nutshell
<br/>

**THE ALLOPIPE TOOL IS DIVIDED INTO TWO DIFFERENT MODULES:**

- **ALLO-COUNT:** performs a stringent data cleaning and a directional comparison of the samples' amino acid sequences.

From two variant-annotated VCF files, variants are first filtered considering a set of quality metrics then constrained to high-confidence calling regions provided in a BED file (GIAB by default). The curated VCF file is then queried for the amino acid information to assess the amino acid mismatches.

**Samples' comparison is directional** and counts either the amino acids that are present in the donor but absent in the recipient (*donor-to-recipient*) or the other way around (*recipient-to-donor*: present by the recipient but absent in the donor.

This step returns the **Allogenomic Mismatch Score (AMS)** which is a discrete quantitative variable measuring the amino acid mismatches in the requested direction, and related information stored in the **AMS-table**.

<br/>

- **ALLO-AFFINITY:** reconstructs peptides around the amino acid changes then return the affinity mHAgs candidates towards HLA molecules thanks to [NetMHCpan softwares](https://pubmed.ncbi.nlm.nih.gov/32406916/)

Allo-Affinity generates a set of candidate minor histocompatibility antigens around each previously assessed directional amino acid mismatches using sliding window. The user defines the length of the potentially HLA-embedded peptides, usually 9-mers for HLA class I and 15-mers for HLA class II molecules. The affinity values are computed using [NetMHCpan4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) and [NETMHCIIpan4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/), respectively.

This step returns the **affinity-AMS (af-AMS)** which is a discrete quantitative variable measuring the amino acid mismatches in the requested mismatch direction, and related information stored in the **AMS-table**. Please note that the 2-field HLA typing has to be provided by the user.

<br/>
<br/>

**THE ALLOPIPE TOOL CAN BE RUN AS SIMPLE- AND MULTI-PROCESSING**
	
- **SIMPLE-PROCESSING**\
To compute AMS et af-AMF just for one pair, mostly if you do at at small scale
 
- **MULTI-PROCESSING**  
It is possible to run both modules from one unique file containing information of more than one pair, i.e. from merged vcf containing information of a batch of samples that represent a cohort.

Command lines are given for each processing mode.

---

## Table of contents

1. [Before getting started](#before)
	1. [Requirements checking](#requirements)
	2. [AlloPipe installation](#install)
	3. [VEP annotation](#vep)
    
<br/>

2. [Run the AlloPipe workflow](#run)
	1. [Launch Allo-Count](#ams_run)
		1. [Getting your Allogenomic Mismatch Score (AMS)](#ams_results)
		2. [Exploring the AMS table](#ams_mismatches)
	2. [Launch Allo-Affinity](#aams_run)
		1. [Getting your affinity-AMS (af-AMS)](#aams_results)
		2. [Exploring the af-AMS table](#aams_mismatches)

<br/>
     
  3. [Tutorial](#tuto)
---

## Before getting started <a name="before"></a>

### Requirements checking <a name="requirements"></a>

AlloPipe specifically requires
1. [Python](https://www.python.org/downloads/) >=3.6 (developed on 3.9)
   
2. [VEP annotation tool](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#download) >=v103
> **_VEP annotation: On-line or command line installation_**\
> VEP annotation can be done using the online tool (if the VCF are smaller than 50 MB) or by downloading the command line tool.
> 
>  - To use the web interface, follow this [link](https://www.ensembl.org/Tools/VEP).
>
> 
>  - To install the command line tool, follow the installation tutorial available [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#download).\
>		During the installation, you will be asked if you want to download **cache** files, **FASTA** files and **plugins**.
>    - We **recommend to download the cache file** for the assembly of your VCF files to be able to run VEP offline.\
	Use the VEP cache version which corresponds to your Ensembl VEP installation and genome reference !
>    - We **recommend to download the FASTA file** for the assembly of your VCF files to be able to run VEP offline.\
	Use the FASTA version which corresponds to your Ensembl VEP installation and genome reference !
>    - We **don't recommend to download any plugin**
>      
> We then recommend to **add VEP to your PATH** by adding the following line to your ```~/.profile``` or ```~/.bash_profile```:\
> 	```export PATH=%%path/to/vep%%:${PATH}```\
>*If you are on Windows, you can follow this [tutorial](https://medium.com/@kevinmarkvi/how-to-add-executables-to-your-path-in-windows-5ffa4ce61a53) to add VEP to your PATH.*
>
>For complete insights on VEP, see [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)

3. [NetMHCpan](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1) and [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) downloaded as a command line tool.\

*Note: you need an academic affiliation (ensured by an email address) to download NetMHCpan softwares from the DTU Health Tech website.*

As we recommend to create a conda environment to ensure a robust installation of AlloPipe, [conda](https://docs.anaconda.com/free/working-with-conda/) should be installed in the suitable version for your operating system and python version.

### AlloPipe installation <a name="install"></a>

To download and install the AlloPipe workflow, first clone the repository from git.\
*(You might be requested to create a token for you to log in. See the [GitHub tutorial](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-personal-access-token-classic))*

We then recommend to create a conda environment dedicated to the AlloPipe workflow. The dependencies specified in the requirements.txt are needed for AlloPipe to run and should be installed in the AlloPipe environment.

The following command lines will perform the above-mentioned steps:

		git clone https://github.com/huguesrichard/Allopipe.git
		cd Allopipe
 		conda create --name Allopipe python=3.9
  		conda activate Allopipe
		python -m pip install -r requirements.txt


### VEP annotation <a name="vep"></a>

AlloPipe input files must be variant-annotated files.\
*We tailored the AlloPipe code based on the VEP annotation architecture, but any other annotation tool could be used after code adjustments.*

Run the following command to annotate you VCF file(s) with VEP.\
**All specified options are mandatory, with the exception of the assembly if you only downloaded one cache file.**  \
	
		vep --cache --assembly GRCh38 --offline --af_gnomade -i %%PATH-TO-FILE-TO-ANNOTATE/FILE-TO-ANNOTATE%%.vcf -o %%PATH-TO-ANNOTATED-FILE/ANNOTATED-FILE%%.vcf --vcf 

This command line works for individual VCF as well as multi-VCF. 
Run this command for every file you want to input in AlloPipe: for individual VCF you will need to run the command twice (once for the donor's VCF and once for the recipient's VCF)


**Once the VEP annotation of your file(s) is complete, you are now ready to launch your first AlloPipe run !**

---

## Run the AlloPipe workflow <a name="run"></a>


### Launch Allo-Count  <a name="ams_run"></a>

Once the VEP annotation is complete, go to the root of the AlloPipe directory to run the following commands in the terminal *(don't forget to activate your conda environment)* :  

**SIMPLE-PROCESSING**

		cd src/
		python ams_pipeline.py -f -n %%NAME-TEST%% -p %%NAME-OF-THE-PAIR%% %%PATH-TO-DONOR-ANNOTATED-FILE/ANNOTATED-FILE%%.vcf %%PATH-TO-RECIPIENT-ANNOTATED-FILE/ANNOTATED-FILE%%.vcf %%DIRECTION OF THE MISMATCH%%

<br/>

Where :\
%%NAME-RUN%% = name of the run\
%%NAME-OF-THE-PAIR% = name of the pair\
%%PATH-TO-DONOR-ANNOTATED-FILE/ANNOTATED-FILE%%.vcf = path to the donor's annotated VCF \
%%PATH-TO-RECIPIENT-ANNOTATED-FILE/ANNOTATED-FILE%%.vcf = path to the recipient's annotated VCF \
%%DIRECTION OF THE MISMATCH%% = 'rd' or 'dr', depending on the direction of the mismach \

>Direction of the mismatch
(see Appendix1: calculation method of AlloPipe)

**MULTI-PROCESSING**

		cd src/
		python multiprocess_ams.py %%PATH-TO-THE-MERGED-ANNOTATED-FILE%%.vcf %%PATH-TO-THE-PAIR-LIST%%.csv %%DIRECTION OF THE MISMATCH%%

Where :
%%PATH-TO-THE-MERGED-ANNOTATED-FILE%%.vcf 
%%PATH-TO-THE-PAIR-LIST%%.csv 
%%DIRECTION OF THE MISMATCH%%


### Getting your ALlogenomic Mismatch Score (AMS) <a name="ams_results"></a>

After the run is complete, have look at the **output/runs/test_run/** directory that was created.  
The directory is structured as followed :  
1. the **AMS/** directory contains a subdirectory created for these run parameters specifically, the AMS value contained in a csv file.  
2. the **plots/** subdirectory
3. the **run_tables** subdirectory contains all the tables created during the run. 


  
We also annotated the tutorial file with VEP for the donor and recipient, `donor.vcf.gz` and `recipient.vcf.gz` respectively.

Run the following command to estimate the AMS for this pair :  

	python ams_pipeline.py -f -n allopipe_run -p test_pair ../tutorial/donor.vcf.gz ../tutorial/recipient.vcf.gz rd

If you have the same AMS for both runs (it should be 49), it means your VEP annotation worked as expected !  

### Exploring the AMS table <a name="ams_mismatches"></a>

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
        - In the case of oncogenetics, Sample_x is the tumor and Sample_y is the constitutional DNA

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


### Launch Allo-Affinity <a name="aams_run"></a>

Once the AMS run is complete, provided you have the class I HLA typing of your samples, you can run a second set of commands to get a filtration of the peptides contributing to the score, using NetMHCpan.  
To run the AAMS pipeline of AlloPipe on the previous example, go to its root directory and run the following commands in the terminal :  

```
	cd src/
	gzip -d ../tutorial/Ensembl/Homo_sapiens.GRCh38.cdna.all.103.fa.gz
 	gzip -d ../tutorial/Ensembl/Homo_sapiens.GRCh38.pep.all.103.fa.gz
  	gzip -d ../tutorial/Ensembl/Homo_sapiens.GRCh38.103.refseq.tsv.gz
	python aams_pipeline.py -M ../output/runs/test_run/run_tables/test_pair_test_run_mismatches_20_400_5_gq_20_0.8_bl_3.tsv \
 	-T ../output/runs/test_run/run_tables/test_pair_test_run_transcripts_pair_codons_20_400_5_gq_20_0.8_bl_3.tsv\
  	-E ../tutorial/Ensembl/Homo_sapiens.GRCh38.cdna.all.103.fa \
   	-P ../tutorial/Ensembl/Homo_sapiens.GRCh38.pep.all.103.fa \
	-R ../tutorial/Ensembl/Homo_sapiens.GRCh38.103.refseq.tsv \
	-n test_run -p test_pair -l 9 --el_rank 2\
	 -a HLA-A*01:01,HLA-A*02:01,HLA-B*08:01,HLA-B*27:05,HLA-C*01:02,HLA-C*07:01
```


### Getting your affinity-AMS (af-AMS) <a name="aams_results"></a>

This second step of AlloPipe uses the AMS information of the first step.  
You will find 3 new subdirectories in the **test_run/** directory :  
1.	the **AAMS/** directory contains a subdirectory created for these run parameters specifically, the AAMS value contained in a csv file.
2.	the **netMHCpan_out/** subdirectory contains all tables generated during the netMHCpan step.
3.	the **aams_run_tables/** subdirectory contains all the other tables created during the run

>![tree_aams_run](tutorial/tree_aams_run.png)

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

## Tutorial <a name="tuto"></a>

We provide a couple of example data in /tutorial, i.e. tutorial/donor_to_annotate and tutorial/donor_to_annotate *(those files correspond to human chr6)*.\
To test your VEP installation, run the following command at the root of the AlloPipe directory :  
	
		vep --cache --assembly GRCh38 --offline --af_gnomade -i tutorial/donor_to_annotate.vcf -o tutorial/donor_annotated_VEP.vcf --vcf --force_overwrite
		vep --cache --assembly GRCh38 --offline --af_gnomade -i tutorial/recipient_to_annotate.vcf -o tutorial/recipient_annotated_VEP.vcf --vcf --force_overwrite 


Once the VEP annotation is complete, go to the root of the AlloPipe directory to run the following commands in the terminal :  

	cd src/
	python ams_pipeline.py -f -n test_run -p test_pair ../tutorial/donor_annotated_VEP.vcf ../tutorial/recipient_annotated_VEP.vcf rd

 If your AMS returns 49, congrats ! You successfully generated your first Allogenomic Mismatch Score (AMS) and related tables !
