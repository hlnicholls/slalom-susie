import pandas as pd
import numpy as np
import subprocess
import re
import scipy as sp
import warnings
warnings.filterwarnings("ignore")


# Set input variables
gwas_file_path = "/Users/hn9/Documents/Analysis/FM-comparison/gwas-examples/APOE-LDL/24097068-GCST002222-EFO_0004611.h.tsv.gz"
target = "APOE_LDL"
target_chrom = 19
target_pos = 44908822
start_pos = target_pos - 500000
end_pos = target_pos + 500000
lead_snp_ID = f"{target_chrom}:{target_pos}:C:T"
n_sample = 94595

def get_ld():
    command = f"plink --bfile /Users/hn9/Documents/Analysis/FM-comparison/ukb_v3_downsampled10k/ukb_v3_chr{target_chrom}.downsampled10k --allow-extra-chr --recode A --chr {target_chrom} --from-bp {start_pos} --to-bp {end_pos} --maf 0.001 --out {target}_locus_UKBB.txt"
    subprocess.run(command, shell=True)
    return

def get_ld_matrix():
    # Calculate LD correlation
    ld_data = pd.read_csv(f"{target}_locus_UKBB.txt.raw", delim_whitespace=True)
    ld_data = ld_data.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    print('Calculating LD correlation matrix...')
    ld_matrix = ld_data.corr(method='pearson')
    print('Finished')
    return ld_data, ld_matrix

def get_sumstats(gwas_file_path, target_chrom, start_pos, end_pos):
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
    return sumstat


def match_snps(sumstat, ld_matrix, ld_data):
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
    return sumstat_filtered, ld_matrix_filtered, concordance_test


def allele_flip_check(concordance_test, sumstat_filtered, ld_matrix_filtered):
    # Create DataFrame for allele flip check
    concordance_test2 = concordance_test.copy()

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

    return sumstat_filtered

def get_lead_ld():
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


def outlier_detection(sumstat, method, r2_threshold = 0.6, nlog10p_dentist_s_threshold = 1e-4):
    if method == 'DENTIST':
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
        lead = df[df['ID'] == lead_snp_ID].iloc[0]
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
    elif method == 'CARMA':
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

def run_finemapping():
    # SuSiE fine-mapping with fine-mapping-inf package
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
    return

# Run functions
get_ld()
ld_data, ld_matrix = get_ld_matrix()
sumstat = get_sumstats(gwas_file_path, target_chrom, start_pos, end_pos)
sumstat_filtered, ld_matrix_filtered, concordance_test = match_snps(sumstat, ld_matrix, ld_data)
sumstat_filtered = allele_flip_check(concordance_test, sumstat_filtered, ld_matrix_filtered)
get_lead_ld()
outlier_detection(sumstat, method='DENTIST', r2_threshold = 0.6, nlog10p_dentist_s_threshold = 1e-4)
run_finemapping()
