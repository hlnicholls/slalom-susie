"""Step to run study locus fine-mapping."""

from __future__ import annotations

from dataclasses import dataclass

import re
import pyspark.sql.functions as F
from scipy import stats

from otg.common.session import Session
from otg.config import FinemappingStepConfig
from otg.dataset.study_locus import StudyLocus

@dataclass
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
        target (str): Target trait or gene for fine-mapping.
        target_chrom (int): Chromosome number where the target is located.
        target_pos (int): Position of the target on the chromosome.
        lead_snp_ID (str): ID of the lead SNP in the locus.
        n_sample (int): Sample size used in the GWAS.
        outlier_method (str): The method to use for outlier detection ('DENTIST' or 'CARMA').
        r2_threshold (float): R-squared threshold for outlier detection.
        nlog10p_dentist_s_threshold (float): Negative log10 p-value threshold for DENTIST method.
        window_size (int, optional): Size of the window around the target position. Default is 500000.
    """
 
    @classmethod
    def get_ld(self):
        """Runs a shell command to generate Linkage Disequilibrium (LD) data using PLINK"""
        command = f"plink --bfile /Users/hn9/Documents/Analysis/FM-comparison/ukb_v3_downsampled10k/ukb_v3_chr{self.target_chrom}.downsampled10k --allow-extra-chr --recode A --chr {self.target_chrom} --from-bp {self.start_pos} --to-bp {self.end_pos} --maf 0.001 --out {self.target}_locus_UKBB.txt"
        subprocess.run(command, shell=True)

    def get_ld_matrix(self):
        """Calculates the LD matrix based on the LD data from PLINK"""
        # Calculate LD correlation
        # TODO debug matrix correlation in PySpark
        # Doesn't finish running
        ld_data = pd.read_csv(f"{self.target}_locus_UKBB.txt.raw", delim_whitespace=True)
        ld_data = ld_data.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
        ld_matrix = ld_data.corr(method='pearson')
        ld_data_spark = spark.read.csv(f"{target}_locus_UKBB.txt.raw", sep=" ", header=True, inferSchema=True)
        drop_list = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
        ld_data_spark = ld_data_spark.drop(*drop_list)
        for col in ld_data_spark.columns:
            ld_data_spark = ld_data_spark.withColumn(col, ld_data_spark[col].cast("float"))
        column_names = ld_matrix.columns.tolist()
        schema = StructType([StructField(name, FloatType(), True) for name in column_names])
        self.ld_matrix_dict = ld_matrix.reset_index(drop=True).to_dict('records')
        self.ld_matrix_spark = spark.createDataFrame(ld_matrix_dict, schema=schema)
        return

    def get_sumstats(self):
        """Filters and renames the GWAS summary statistics columns."""
        gwas_file_path = self.gwas_file_path
        target_chrom = self.target_chrom
        start_pos = self.start_pos
        end_pos = self.end_pos
        # Preparing sumstats for filtering and input into SuSiE
        sumstat = spark.read.csv(gwas_file_path, header=True, sep="\t")
        sumstat = sumstat.filter(
            (col('chromosome') == target_chrom) &
            (col('position') >= start_pos) &
            (col('position') <= end_pos)
        )
        sumstat = sumstat.dropna(subset=['chromosome'])

        sumstat = sumstat.withColumn('z', col('beta') / col('standardError'))
        # making and renaming columns to match col names needed to run susieinf package
        # Use variantId to create allele columns
        cols_to_rename = {
            'hm_variant_id': 'ID',
            'hm_rsid': 'rsid',
            'hm_pos': 'position',
            'hm_other_allele': 'allele1',
            'hm_effect_allele': 'allele2',
            'hm_effect_allele_frequency': 'maf',
            'standardError': 'se',
            'p_value': 'p'
        }
        for old_name, new_name in cols_to_rename.items():
            sumstat = sumstat.withColumnRenamed(old_name, new_name)
        selected_cols = ['ID', 'rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'p', 'beta', 'se', 'z']
        sumstat = sumstat.select(selected_cols)
        self.sumstat = sumstat
        return

    def match_snps(self):
        """Matches SNPs between summary statistics and LD matrix and filters them accordingly."""
        sumstat = self.sumstat
        ld_matrix = self.ld_matrix
        ld_data = self.ld_data
        # Getting only SNPs in sumstats that are also in the LD matrix
        # Get SNP IDs from ld matrix to compare with sumstat SNP IDs
        # Regex pattern for extracting position from SNP
        pattern = re.compile(r"(^\d+)|(?<=:)\d+(?=:|$)")
        ld_data_spark = ld_data_spark.withColumn("position", F.regexp_extract(F.col("SNP"), pattern, 0))
        ld_data_spark = ld_data_spark.withColumn("ID", F.regexp_replace(F.col("SNP"), r'[:,_]', '_'))
        
        # Join sumstat with ld_data on 'ID'
        concordance_test = sumstat.join(ld_data_spark, 'ID', 'inner')
        
        # Filter sumstat and ld_matrix for ID matches only
        sumstat_filtered = sumstat.filter(F.col("ID").isin(concordance_test.select("ID").distinct().rdd.flatMap(lambda x: x).collect()))
        ld_matrix_filtered = ld_matrix_spark.filter(F.col("ID").isin(concordance_test.select("ID").distinct().rdd.flatMap(lambda x: x).collect()))

        return

    def allele_flip_check_pyspark(self):
        """Checks for allele flips between summary statistics and LD matrix and adjusts z-scores accordingly."""
        sumstat_filtered = self.sumstat_filtered  # Assuming this is a PySpark DataFrame
        ld_matrix_filtered = self.ld_matrix_filtered  # Assuming this is a PySpark DataFrame
        target = self.target  # Assuming this is a string

        # Create DataFrame for allele flip check
        concordance_test2 = self.concordance_test

        # Extract alleles
        allele_df = concordance_test2.withColumn("allele_split", F.split("SNP", r"[:,_]"))
        concordance_test2 = allele_df.withColumn("allele1_LD", allele_df["allele_split"].getItem(1)) \
                                    .withColumn("allele2_LD", allele_df["allele_split"].getItem(2))

        # Assumed that indices are the same in sumstats and concordance_test
        # Flip z-scores if alleles are discordant between sumstats and LD matrix
        condition = (sumstat_filtered["allele1"] != concordance_test2["allele1_LD"]) | \
                    (sumstat_filtered["allele2"] != concordance_test2["allele2_LD"])
        sumstat_filtered = sumstat_filtered.withColumn("z", F.when(condition, -sumstat_filtered["z"]).otherwise(sumstat_filtered["z"]))
        ld_matrix_filtered.write.csv(f"{target}_locus_ukbb_ld.txt.gz", sep='\t', header=False)
        sumstat_filtered.write.csv(f"{target}_locus_sumstat_flip_check.txt.gz", sep='\t', header=True)

        return

    def get_lead_ld(self):
        """Runs a shell command to get the LD of the lead SNP using PLINK."""
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


    def outlier_detection_pyspark(self):
        """Performs outlier detection using the specified method."""
        # Assuming these are already defined and initialized elsewhere in your code
        sumstat = self.sumstat
        outlier_method = self.outlier_method
        target = self.target
        n_sample = self.n_sample
        lead_snp_ID = self.lead_snp_ID
        r2_threshold = self.r2_threshold
        nlog10p_dentist_s_threshold = self.nlog10p_dentist_s_threshold

        if outlier_method == 'DENTIST':
            ld = spark.read.csv(f'/path/to/{target}_subset_for_ld_calculation.ld', sep='\s+', header=True)
            sumstat = spark.read.csv(f'/path/to/{target}_locus_sumstat_flip_check.txt.gz', sep='\t', header=True)

            # 1. Getting R2 column for sumstats
            lead_ld = ld.filter((ld['SNP_A'] == lead_snp_ID) | (ld['SNP_B'] == lead_snp_ID))
            sumstat = sumstat.withColumn('ID', F.regexp_replace('ID', "_", ":"))
            merged = lead_ld.join(sumstat.select('ID'), lead_ld.SNP_B == sumstat.ID)
            df = merged.select('ID', 'R2')
            df = sumstat.join(df, on='ID', how='left')

            # Calculate 'r'
            df = df.withColumn('r', (F.sum('R2') * n_sample) / (F.count('R2') * n_sample))

            lead_idx_snp_row = df.filter(df.ID == lead_snp_ID).collect()[0]
            lead_z = lead_idx_snp_row.beta / lead_idx_snp_row.se

            # 2. Calculate 't_dentist_s' and 'dentist_outlier'
            df = df.withColumn("t_dentist_s", ((df.beta / df.se - df.r * lead_z)**2) / (1 - df.r**2))
            df = df.withColumn("t_dentist_s", F.when(df["t_dentist_s"] < 0, float('inf')).otherwise(df["t_dentist_s"]))
            def calc_nlog10p_dentist_s(t_dentist_s):
                return stats.chi2.logsf(t_dentist_s, df=1) / -np.log(10)
            calc_nlog10p_dentist_s_udf = F.udf(calc_nlog10p_dentist_s, DoubleType())
            df = df.withColumn("nlog10p_dentist_s", calc_nlog10p_dentist_s_udf("t_dentist_s"))

            # Count the number of DENTIST outliers
            n_dentist_s_outlier = df.filter((df.R2 > r2_threshold) & (df.nlog10p_dentist_s > nlog10p_dentist_s_threshold)).count()
            print(f'Number of DENTIST outliers detected: {n_dentist_s_outlier}')
            df = df.withColumn('dentist_outlier', F.when((df.R2 > r2_threshold) & (df.nlog10p_dentist_s > nlog10p_dentist_s_threshold), 1).otherwise(0))
            df = df.drop('CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B')
            df.write.csv(f'{target}_locus_sumstat_with_dentist.txt.gz', sep='\t', header=True)

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
        return

    def run_pipeline(self):
        self.get_ld()
        self.get_ld_matrix()
        self.get_sumstats()
        self.match_snps()
        self.allele_flip_check()
        self.get_lead_ld()
        self.outlier_detection()
        self.run_finemapping()

if __name__ == "__main__":
    pipeline = FinemappingPipeline(
        gwas_file_path="/Users/hn9/Documents/Analysis/FM-comparison/gwas-examples/APOE-LDL/24097068-GCST002222-EFO_0004611.h.tsv.gz",
        target="APOE_LDL",
        target_chrom=19,
        target_pos=44908822,
        lead_snp_ID="19:44908822:C:T",
        n_sample=94595,
        outlier_method='DENTIST',
        r2_threshold = 0.6,
        nlog10p_dentist_s_threshold = 1e-4
    )
    pipeline.run_pipeline()
