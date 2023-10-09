"""Step to run study locus fine-mapping."""

from __future__ import annotations

from dataclasses import dataclass

import re
import pyspark.sql.functions as F
from scipy import stats
import hail as hl
from hail.linalg import BlockMatrix

from otg.common.session import Session
from otg.config import FinemappingStepConfig
from otg.dataset.study_locus import StudyLocus

from otg.common.utils import _liftover_loci, convert_gnomad_position_to_ensembl
from otg.datasource.gnomad import GnomADLDMatrix

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

    @staticmethod
    def lift_and_calculate_ld_gnomad(variant_id, chain_file_38_to_37, chain_file_37_to_38):
        bm = BlockMatrix.read("gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm")
        variant_table = hl.read_table("gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.nfe.common.adj.ld.variant_indices.ht")
        chain_file_38_to_37 = 'gs://hail-common/references/grch38_to_grch37.over.chain.gz'
        chain_file_37_to_38 = 'gs://hail-common/references/grch37_to_grch38.over.chain.gz'
        # Lift over the lead variant from 38 to 37
        locus = hl.parse_variant(variant_id)
        locus_37 = hl.get_reference('GRCh37').parse_locus(str(locus))
        lifted_locus = hl.liftover(locus_37, chain_file_38_to_37, include_strand=True)

        # Get locus from gnomad LD matrix
        window_size = 500000
        relevant_variants = variant_table.filter(
            (hl.abs(variant_table.locus.position - lifted_locus.result.position) <= window_size) &
            (variant_table.locus.contig == lifted_locus.result.contig)
        )
        indices = relevant_variants['index'].collect()
        sub_bm = bm.filter(indices, indices)

        # Lift over variants from 37 to 38
        lifted_variants = [
            hl.liftover(hl.locus(relevant_variants.locus.contig, x), chain_file_37_to_38).result
            for x in relevant_variants.locus.position.collect()
        ]

        # Get all SNP IDs in locus
        snp_ids = [
            f"{locus.contig}_{locus.position}_{alleles[0]}_{alleles[1]}"
            for locus, alleles in lifted_variants
        ]

        ld_df = hl.block_matrix_to_pandas(sub_bm)
        ld_pyspark_df = hl.pandas_to_spark(ld_df)

        return ld_pyspark_df, snp_ids


    def get_ld_matrix(self):
        """Get LD correlation matrix from hail"""
        bm = BlockMatrix.read("gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm")
        variant_table = hl.read_table("gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.nfe.common.adj.ld.variant_indices.ht")
        lead_snp = hl.parse_variant("1:90277796:G:A")
        window_size = 500000
        filtered_matrix = variant_table.filter(
            (variant_table.locus.contig == lead_snp.locus.contig) &
            (hl.abs(variant_table.locus.position - lead_snp.locus.position) <= window_size)
        )
        self.variant_ids = filtered_matrix.select(
            variant_id=hl.str(filtered_matrix.locus) + ":" + hl.delimit(filtered_matrix.alleles, ":")
        ).variant_id.collect()
        indices = filtered_matrix.idx.collect()
        subset_bm = bm.filter(indices, indices)
        mt_from_bm = subset_bm.to_matrix_table_row_major()
        mt_from_bm = mt_from_bm.entries()
        self.ld_matrix = mt_from_bm.to_spark(flatten=True)
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
        variant_ids = self.variant_ids
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
            # get study locus with columsn: r2 with lead snp, beta, se, z
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
