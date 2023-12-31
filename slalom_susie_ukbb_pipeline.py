import pandas as pd
import numpy as np
import subprocess
import re
import scipy as sp
import warnings
warnings.filterwarnings("ignore")


class FinemappingPipeline:
    """Fine-mapping pipeline for an input locus"""
    def __init__(self, gwas_file_path, target, target_chrom, target_pos, lead_snp_ID,
                  n_sample, outlier_method, r2_threshold, nlog10p_dentist_s_threshold, window_size=500000):
        self.gwas_file_path = gwas_file_path
        self.target = target
        self.target_chrom = target_chrom
        self.target_pos = target_pos
        self.start_pos = target_pos - window_size
        self.end_pos = target_pos + window_size
        self.lead_snp_ID = lead_snp_ID
        self.n_sample = n_sample
        self.outlier_method = outlier_method
        self.r2_threshold = r2_threshold
        self.nlog10p_dentist_s_threshold = nlog10p_dentist_s_threshold
    """
    Args:
        gwas_file_path (str): Path to the GWAS summary statistics file.
        target (str): Locus name (expected in name of all files for a locus)
        target_chrom (int): Chromosome number where the target is located.
        target_pos (int): Position of the target on the chromosome.
        lead_snp_ID (str): ID of the lead SNP in the locus - use ":" separator.
        n_sample (int): Sample size used in the GWAS.
        outlier_method (str): The method to use for outlier detection ('DENTIST' or 'CARMA').
        r2_threshold (float): R-squared threshold for DENTIST outlier detection.
        nlog10p_dentist_s_threshold (float): Threshold for DENTIST outlier detection.
        window_size (int, optional): Size of the window around the target position. Default is 500kb.

    Methods:
        get_ld(): Get raw LD data for a locus (currently from UKBB using PLINK).
        get_ld_matrix_cor(): Calculate LD matrix correlation (pearsons R2).
        get_sumstats(): Get locus sumstats and rename columns (need specific names for SuSiE-inf package).
        match_snps(): Get only SNPs in sumstats that are also present in the LD matrix and vice versa.
        allele_flip_check(): Check if allele1 and allele2 order matches between sumstats and LD matrix. If not then flip sumstat Z-score sign.
        get_lead_ld(): Get LD for all variants with lead SNP (needed for DENTIST outlier dection).
        outlier_detection(): Run either DENTIST or CARMA outlier detection (CARMA work in progress).
        run_finemapping(): Run fine-mapping method (SuSiE-inf).
    """

    def get_ld(self):
        """Runs a shell command to generate Linkage Disequilibrium (LD) data using PLINK"""
        command = f"plink --bfile /Users/hn9/Documents/Analysis/FM-comparison/ukb_v3_downsampled10k/ukb_v3_chr{self.target_chrom}.downsampled10k --allow-extra-chr --recode A --chr {self.target_chrom} --from-bp {self.start_pos} --to-bp {self.end_pos} --maf 0.001 --out {self.target}_locus_UKBB.txt"
        subprocess.run(command, shell=True)

    def get_ld_matrix_cor(self):
        """Calculates the LD matrix correlation based on the LD data from PLINK"""
        # Calculate LD correlation
        ld_data = pd.read_csv(f"{self.target}_locus_UKBB.txt.raw", delim_whitespace=True)
        ld_data = ld_data.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
        print('Calculating LD correlation matrix...')
        ld_matrix = ld_data.corr(method='pearson')
        print('Finished calculating LD correlation matrix')
        self.ld_data = ld_data
        self.ld_matrix = ld_matrix
        return

    def get_sumstats(self):
        """Filters and renames the GWAS summary statistics columns (needed for SuSiE)."""
        gwas_file_path = self.gwas_file_path
        target_chrom = self.target_chrom
        start_pos = self.start_pos
        end_pos = self.end_pos
        # Preparing sumstats for filtering and input into SuSiE
        sumstat = pd.read_csv(gwas_file_path, sep="\t")
        sumstat = sumstat[(sumstat['hm_chrom'] == target_chrom) & (sumstat['hm_pos'] >= start_pos) & (sumstat['hm_pos'] <= end_pos)]
        sumstat = sumstat.dropna(subset=['hm_chrom'])
        sumstat['z'] = sumstat['beta'] / sumstat['standard_error']
        cols_to_rename = {
            'hm_variant_id': 'ID',
            'hm_rsid': 'rsid',
            'hm_chrom': 'chromosome',
            'hm_pos': 'position',
            'hm_other_allele': 'allele1',
            'hm_effect_allele': 'allele2',
            'hm_effect_allele_frequency': 'maf',
            'standard_error': 'se',
            'p_value': 'p'
        }
        sumstat.rename(columns=cols_to_rename, inplace=True)
        sumstat = sumstat[['ID', 'rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'p', 'beta', 'se', 'z']]
        self.sumstat = sumstat
        return

    def match_snps(self):
        """Matches SNPs between summary statistics and LD matrix and filters them accordingly."""
        sumstat = self.sumstat
        ld_matrix = self.ld_matrix
        ld_data = self.ld_data
        # Getting only SNPs in sumstats that are also in the LD matrix
        # Get SNP IDs from ld matrix to compare with sumstat SNP IDs
        pattern = re.compile(r"(^\d+)|(?<=:)\d+(?=:|$)")
        df1_transpose = ld_data.T.reset_index()
        df1_transpose.columns = ['SNP'] + list(df1_transpose.columns[1:])
        df1_transpose['position'] = df1_transpose['SNP'].apply(lambda x: re.search(pattern, x).group())
        df1_transpose['ID'] = df1_transpose['SNP'].str.replace(r'[:,_]', '_').str.replace(r'_[^_]+$', '')
        concordance_test = pd.merge(sumstat, df1_transpose, on='ID')

        # Filter sumstat and LD matrix for ID matches only
        sumstat_filtered = sumstat[sumstat['ID'].isin(concordance_test['ID'])]
        sumstat_filtered.reset_index(drop=True, inplace=True)
        ld_matrix_filtered = ld_matrix.loc[concordance_test['SNP'], concordance_test['SNP']]
        self.sumstat_filtered = sumstat_filtered
        self.ld_matrix_filtered = ld_matrix_filtered
        self.concordance_test = concordance_test
        return

    def allele_flip_check(self):
        """Checks for allele flips between summary statistics and LD matrix and adjusts z-scores accordingly."""
        sumstat_filtered = self.sumstat_filtered
        ld_matrix_filtered = self.ld_matrix_filtered
        target = self.target
        # Create DataFrame for allele flip check
        concordance_test2 = self.concordance_test.copy()

        # Extract alleles
        allele_df = concordance_test2['SNP'].str.extract(r'[:,_]([ACGT]+)[:,_]([ACGT]+)')
        concordance_test2['allele1_LD'] = allele_df[0]
        concordance_test2['allele2_LD'] = allele_df[1]

        # Indices need to be aligned
        sumstat_filtered.index = concordance_test2.index

        # Flip z-scores if alleles are discordant between sumstats and LD matrix
        sumstat_filtered['z'] = np.where(
            (sumstat_filtered['allele1'] != concordance_test2['allele1_LD']) | (sumstat_filtered['allele2'] != concordance_test2['allele2_LD']),
            -sumstat_filtered['z'], sumstat_filtered['z']
        )

        # Check for NaNs or Inf values (not accepted in fine-mapping)
        numeric_cols = sumstat_filtered.select_dtypes(include=[np.number])
        any_na_matrix = ld_matrix_filtered.isna().any().any()
        any_inf_matrix = np.isinf(ld_matrix_filtered.to_numpy()).any()
        any_na_sumstat = numeric_cols.isna().any().any()
        any_inf_sumstat = np.isinf(numeric_cols.to_numpy()).any()

        print(f"Correlation matrix has NaNs: {any_na_matrix}, Infs: {any_inf_matrix}")
        print(f"Sumstat has NaNs: {any_na_sumstat}, Infs: {any_inf_sumstat}")
        ld_matrix_filtered.to_csv(f"{target}_locus_ukbb_ld.txt.gz", sep='\t', index=False, header=False)
        sumstat_filtered.to_csv(f"{target}_locus_sumstat_flip_check.txt.gz", sep='\t', index=False)

        return

    def get_lead_ld(self):
        """Runs a shell command to get the LD with just the lead SNP using PLINK."""
        target = self.target
        target_chrom = self.target_chrom
        lead_snp_ID = self.lead_snp_ID
        lead_ld_command = f"""plink --bfile /Users/hn9/Documents/Analysis/FM-comparison/ukb_v3_downsampled10k/ukb_v3_chr{target_chrom}.downsampled10k \
            --allow-extra-chr \
            --r2 \
            --ld-snp {lead_snp_ID} \
            --ld-window-kb 1000 \
            --ld-window 99999 \
            --ld-window-r2 0 \
            --out /Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/ld/{target}_subset_for_ld_calculation
        """
        subprocess.run(lead_ld_command, shell=True)
        return

    def outlier_detection(self):
        """Performs outlier detection using the specified method."""
        sumstat = self.sumstat
        outlier_method = self.outlier_method
        target = self.target
        n_sample = self.n_sample
        lead_snp_ID = self.lead_snp_ID
        r2_threshold = self.r2_threshold
        nlog10p_dentist_s_threshold = self.nlog10p_dentist_s_threshold
        if outlier_method == 'DENTIST':
            # 1. Getting R2 column for sumstats
            # LD for all variants with lead SNP needs to be previously calculated by plink
            ld = pd.read_csv(f'/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/ld/{target}_subset_for_ld_calculation.ld', sep='\s+')
            lead_ld = ld[(ld['SNP_A'] == f'{lead_snp_ID}') | (ld['SNP_B'] == f'{lead_snp_ID}')]
            sumstat = pd.read_csv(f'/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/{target}_locus_sumstat_flip_check.txt.gz', sep='\t')
            sumstat['ID'] = sumstat['ID'].str.replace("_", ":")
            merged = pd.merge(lead_ld, sumstat[['ID']], left_on='SNP_B', right_on='ID')
            df = merged[['ID', 'R2']]
            df = pd.merge(sumstat, merged, on='ID', how='left')
            df['r'] = (n_sample * df['R2'].sum()) / (n_sample * df['R2'].count())
            lead_idx_snp = df.index[df['ID'] == lead_snp_ID].tolist()[0]

            # 2. Calculate 't_dentist_s' and 'dentist_outlier'
            # Copied and unaltered from https://github.com/mkanai/slalom/blob/master/slalom.py
            lead_z = (df.beta / df.se).iloc[lead_idx_snp]
            df["t_dentist_s"] = ((df.beta / df.se) - df.r * lead_z) ** 2 / (1 - df.r ** 2)
            df["t_dentist_s"] = np.where(df["t_dentist_s"] < 0, np.inf, df["t_dentist_s"])
            df["t_dentist_s"].iloc[lead_idx_snp] = np.nan
            df["nlog10p_dentist_s"] = sp.stats.chi2.logsf(df["t_dentist_s"], df=1) / -np.log(10)
            n_dentist_s_outlier = np.sum(
                        (df.R2 > r2_threshold) & (df.nlog10p_dentist_s > nlog10p_dentist_s_threshold)
                    )
            print('Number of DENTIST outliers detected:', n_dentist_s_outlier)
            # Identifying outliers
            df['dentist_outlier'] = np.where((df.R2 > r2_threshold) & (df.nlog10p_dentist_s > nlog10p_dentist_s_threshold), 1, 0)
            df = df.drop(['CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B'], axis=1)
            df.to_csv(f'{target}_locus_sumstat_with_dentist.txt.gz', sep='\t', index=False)

        elif outlier_method == 'CARMA':
            # unfinished work in progress
            from scipy import optimize
            import carmapy.carmapy_c
            sumstats = pd.read_csv(f"{target}_locus_sumstat_flip_check.txt.gz", sep='\t')
            ld = pd.read_csv(f"{target}_locus_ukbb_ld.txt.gz", sep='\t', header=None)
            outlier_tau = 0.04
            index_list = sumstats.index.tolist()
            z = sumstats['z'].tolist()
            ld_matrix = np.asmatrix(ld)
            modi_ld_S = ld_matrix

            def ridge_fun(x, modi_ld_S, index_list, temp_Sigma, z, outlier_tau, outlier_likelihood):
                temp_ld_S = x * modi_ld_S + (1 - x) * np.eye(modi_ld_S.shape[0])
                ld_matrix[index_list[:, None], index_list] = temp_ld_S
                return outlier_likelihood(index_list, ld_matrix, z, outlier_tau, len(index_list), 1)

            opizer = optimize(ridge_fun, interval=[0, 1], maximum=True)
            modi_ld = opizer['maximum'] * modi_ld_S + (1 - opizer['maximum']) * np.diag(np.diag(modi_ld_S))
            outlier_likelihood = carmapy.carmapy_c.outlier_Normal_fixed_sigma_marginal
            test_log_BF = outlier_likelihood(index_list, ld_matrix, z, outlier_tau, len(index_list), 1) - outlier_likelihood(index_list, modi_ld, z, outlier_tau, len(index_list), 1)
            test_log_BF = -abs(test_log_BF)
            print('Outlier BF:', test_log_BF)
            print('This is xi hat:', opizer)
        else:
            pass
        return

    def run_finemapping(self):
        """Runs the SuSiE fine-mapping algorithm from fine-mapping-inf package."""
        target = self.target
        n_sample = self.n_sample
        susieinf_command = f"""python /Users/hn9/Documents/GitHub/fine-mapping-inf/run_fine_mapping.py \
            --sumstats /Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/{target}_locus_sumstat_with_dentist.txt.gz \
            --beta-col-name beta \
            --se-col-name se \
            --ld-file /Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/{target}_locus_ukbb_ld.txt.gz \
            --n {n_sample} \
            --method susieinf \
            --save-tsv \
            --eigen-decomp-prefix {target}_locus \
            --output-prefix  {target}_locus """

        subprocess.run(susieinf_command, shell=True)
        #zip_comamnd = f"""gunzip -c {target}_locus.susieinf.bgz > {target}_locus.txt"""
        #subprocess.run(zip_comamnd, shell=True)
        return

    def run_pipeline(self):
        self.get_ld()
        self.get_ld_matrix_cor()
        self.get_sumstats()
        self.match_snps()
        self.allele_flip_check()
        self.get_lead_ld()
        self.outlier_detection()
        self.run_finemapping()

if __name__ == "__main__":
    pipeline = FinemappingPipeline(
        gwas_file_path="/Users/hn9/Documents/Analysis/FM-comparison/gwas-examples/PTK2B-Alzh/33589840-GCST90012877-EFO_0000249.h.tsv.gz",
        target="PTK2B",
        target_chrom=8,
        target_pos=27610986,
        lead_snp_ID="8:27610986:C:A",
        n_sample=472868,
        outlier_method='DENTIST',
        r2_threshold = 0.6,
        nlog10p_dentist_s_threshold = 1e-4
    )
    pipeline.run_pipeline()
