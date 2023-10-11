import pandas as pd
import numpy as np

import hail as hl
from hail.linalg import BlockMatrix
hl.init()

from pyspark.sql import SparkSession
spark = SparkSession.builder \
    .appName("Finemapping Pipeline") \
    .getOrCreate()

def get_gnomad_ld_matrix(lead_snp_ID):
    rg38 = hl.get_reference('GRCh38')
    rg37 = hl.get_reference('GRCh37')
    rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)
    rg38.add_liftover('gs://hail-common/references/grch38_to_grch37.over.chain.gz', rg37)
    bm = BlockMatrix.read("gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm")
    variant_table = hl.read_table("gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.nfe.common.adj.ld.variant_indices.ht")
    # Liftover lead SNP ID to GRCh37
    locus = hl.parse_variant(lead_snp_ID,reference_genome='GRCh38').locus
    temp_table = hl.Table.parallelize([hl.struct(locus=locus)])
    locus_values = temp_table.select('locus').collect()
    locus_value = locus_values[0].locus
    contig = locus_value.contig
    position = locus_value.position
    locus_37 = hl.liftover(hl.locus(contig, position, 'GRCh38'), 'GRCh37')

    # Get LD matrix
    window_size = 500000
    locus_variants = variant_table.filter(
        (hl.abs(variant_table.locus.position - locus_37.position) <= window_size) &
        (variant_table.locus.contig == locus_37.contig)
    )
    indices = locus_variants['idx'].collect()
    sub_bm = bm.filter(indices, indices)
    numpy_array = sub_bm.to_numpy()
    ld_df = pd.DataFrame(numpy_array)
    # need to change to iteritems due to old pandas version error
    ld_df.iteritems = ld_df.items

    # Get SNP IDs in GRCh38
    locus_variants = locus_variants.annotate(
        locus_38 = hl.liftover(locus_variants.locus, 'GRCh38')
    )
    locus_variants = locus_variants.annotate(
        snp_id_38 = hl.str(locus_variants.locus_38.contig) + "_" +
                    hl.str(locus_variants.locus_38.position) + "_" +
                    locus_variants.alleles[0] + "_" +
                    locus_variants.alleles[1]
    )
    snp_ids_38 = locus_variants['snp_id_38'].collect()
    ld_df.columns = snp_ids_38
    ld_pyspark_df = spark.createDataFrame(ld_df)

    return ld_pyspark_df, snp_ids_38


def get_matching_snps(ld_pyspark_df, snp_ids_38):
    study_snps = [row.snp_id_38 for row in StudyLocus.select("variantID").collect()]
    study_snps = set(study_snps)
    # Find the intersection of the SNPs
    common_snps = study_snps.intersection(snp_ids_38)
    
    # Filter StudyLocus to only include common SNPs
    filtered_StudyLocus = StudyLocus.filter(F.col("snp_id_38").isin(common_snps))
    
    # Filter LD matrix to only include the common SNPs
    selected_columns = [col for col in ld_pyspark_df.columns if col in common_snps]
    filtered_ld_pyspark_df = ld_pyspark_df.select(selected_columns)
    
    return filtered_ld_pyspark_df, filtered_StudyLocus

ld_pyspark_df, snp_ids_38 = get_gnomad_ld_matrix("chr8:27610986:C:A")