# Parallelised Fine-mapping Pipeline

Scripts for running fine-mapping QC steps and fine-mapping via SuSiE or CARMA.

### Pipeline steps (ran in all scripts):
1. Get intersecting SNPs between summary statistics and LD matrix
2. Check allele flipping
3. DENTIST outlier detection
4. Fine-mapping method (SuSiE or CARMA)

Scripts:
- `susie_single_locus.R` - Runs pipeline steps for an individual locus.
- `susie_loci_for.R` - Runs pipeline steps iteratively for a list of loci one by one.
- `susie_loci_paralle.R` and `CARMA_loci_paralle.R` - Runs pipeline steps iteratively for batches of loci.

### Requirements:
- `loci_list.txt` - a file that contains all the lead SNPs for the loci. This is used as a file name identifier.
    - The scripts expect the sumstats and LD matrix file names to have the lead SNP ID in them (e.g., loci_list has an ID column for lead SNPs like "1_1234_C_G" so we'd have files like "1_1234_C_G_LD_matrix.txt.gz" and "1_1234_C_G_sumstats.txt.gz") then it will output results like "1_1234_C_G_susie_results.txt.gz"

- PLINK installed via an anaconda environment (install via bioconda or use provided: https://github.com/hlnicholls/slalom-susie/blob/main/finemapping_env.yml).

- Working directories need to be changed in the scripts before running.

- Check column namnes referred to in the script match the column names of the input data (e.g. 'chromosome' versus 'chrom' etc.)

- Additional considerations: sample size is set for DENTIST and susie. DENTIST thresholds are currently those used in SLALOM.