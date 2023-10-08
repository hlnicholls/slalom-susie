"""Config management for OTG."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from hydra.core.config_store import ConfigStore
from omegaconf import MISSING


@dataclass
class Config:
    """Configuration for otg."""

    defaults: List[Dict[str, str]] = field(
        default_factory=lambda: [{"step": "???"}, {"session": "session_config"}]
    )

    step: Any = MISSING
    session: Any = MISSING


@dataclass
class SessionConfig:
    """ETL config."""

    _target_: str = "otg.common.session.Session"
    app_name: str = "otgenetics"
    spark_uri: str = "local[*]"
    write_mode: str = "overwrite"
    hail_home: Optional[str] = None


@dataclass
class LDIndexStepConfig:
    """LD index step requirements.

    Attributes:
        pop_ldindex_path (str): Input population LD index file from gnomAD.
        ld_radius (int): Window radius around locus.
        grch37_to_grch38_chain_path (str): Path to GRCh37 to GRCh38 chain file.
        ld_index_path (str): Output LD index path.
    """

    _target_: str = "otg.ld_index.LDIndexStep"
    ld_index_raw_template: str = "gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.{POP}.common.ld.variant_indices.ht"
    ld_radius: int = 500_000
    grch37_to_grch38_chain_path: str = MISSING
    ld_index_template: str = MISSING
    ld_populations: List[str] = field(
        default_factory=lambda: [
            "afr",  # African-American
            "amr",  # American Admixed/Latino
            "asj",  # Ashkenazi Jewish
            "eas",  # East Asian
            "fin",  # Finnish
            "nfe",  # Non-Finnish European
            "nwe",  # Northwestern European
            "seu",  # Southeastern European
        ]
    )


@dataclass
class VariantIndexStepConfig:
    """Variant index step requirements.

    Attributes:
        variant_annotation_path (str): Input variant annotation path.
        study_locus_path (str): Input study-locus path.
        variant_index_path (str): Output variant index path.
    """

    _target_: str = "otg.variant_index.VariantIndexStep"
    variant_annotation_path: str = MISSING
    study_locus_path: str = MISSING
    variant_index_path: str = MISSING


@dataclass
class ColocalisationStepConfig:
    """Colocalisation step requirements.

    Attributes:
        study_locus_path (DictConfig): Input Study-locus path.
        coloc_path (DictConfig): Output Colocalisation path.
        priorc1 (float): Prior on variant being causal for trait 1.
        priorc2 (float): Prior on variant being causal for trait 2.
        priorc12 (float): Prior on variant being causal for traits 1 and 2.
    """

    _target_: str = "otg.colocalisation.ColocalisationStep"
    study_locus_path: str = MISSING
    study_index_path: str = MISSING
    coloc_path: str = MISSING
    priorc1: float = 1e-4
    priorc2: float = 1e-4
    priorc12: float = 1e-5


@dataclass
class VariantAnnotationStepConfig:
    """Variant annotation step requirements.

    Attributes:
        gnomad_genomes (str): Path to gnomAD genomes hail table.
        chain_38_to_37 (str): Path to GRCh38 to GRCh37 chain file.
        variant_annotation_path (str): Output variant annotation path.
        populations (List[str]): List of populations to include.
    """

    _target_: str = "otg.variant_annotation.VariantAnnotationStep"
    gnomad_genomes: str = MISSING
    chain_38_to_37: str = MISSING
    variant_annotation_path: str = MISSING
    populations: List[str] = field(
        default_factory=lambda: [
            "afr",  # African-American
            "amr",  # American Admixed/Latino
            "ami",  # Amish ancestry
            "asj",  # Ashkenazi Jewish
            "eas",  # East Asian
            "fin",  # Finnish
            "nfe",  # Non-Finnish European
            "mid",  # Middle Eastern
            "sas",  # South Asian
            "oth",  # Other
        ]
    )


@dataclass
class V2GStepConfig:
    """Variant to gene (V2G) step requirements.

    Attributes:
        variant_index_path (str): Input variant index path.
        variant_annotation_path (str): Input variant annotation path.
        gene_index_path (str): Input gene index path.
        vep_consequences_path (str): Input VEP consequences path.
        lift_over_chain_file_path (str): Path to GRCh37 to GRCh38 chain file.
        approved_biotypes (list[str]): List of approved biotypes.
        anderson_path (str): Anderson intervals path.
        javierre_path (str): Javierre intervals path.
        jung_path (str): Jung intervals path.
        thurnman_path (str): Thurnman intervals path.
        liftover_max_length_difference (int): Maximum length difference for liftover.
        max_distance (int): Maximum distance to consider.
        output_path (str): Output V2G path.
    """

    _target_: str = "otg.v2g.V2GStep"
    variant_index_path: str = MISSING
    variant_annotation_path: str = MISSING
    gene_index_path: str = MISSING
    vep_consequences_path: str = MISSING
    liftover_chain_file_path: str = MISSING
    anderson_path: str = MISSING
    javierre_path: str = MISSING
    jung_path: str = MISSING
    thurnman_path: str = MISSING
    liftover_max_length_difference: int = 100
    max_distance: int = 500_000
    v2g_path: str = MISSING
    approved_biotypes: List[str] = field(
        default_factory=lambda: [
            "protein_coding",
            "3prime_overlapping_ncRNA",
            "antisense",
            "bidirectional_promoter_lncRNA",
            "IG_C_gene",
            "IG_D_gene",
            "IG_J_gene",
            "IG_V_gene",
            "lincRNA",
            "macro_lncRNA",
            "non_coding",
            "sense_intronic",
            "sense_overlapping",
        ]
    )


@dataclass
class GWASCatalogStepConfig:
    """GWAS Catalog step requirements.

    Attributes:
        catalog_studies_file (str): Raw GWAS catalog studies file.
        catalog_ancestry_file (str): Ancestry annotations file from GWAS Catalog.
        catalog_sumstats_lut (str): GWAS Catalog summary statistics lookup table.
        catalog_associations_file (str): Raw GWAS catalog associations file.
        variant_annotation_path (str): Input variant annotation path.
        ld_populations (list): List of populations to include.
        min_r2 (float): Minimum r2 to consider when considering variants within a window.
        catalog_studies_out (str): Output GWAS catalog studies path.
        catalog_associations_out (str): Output GWAS catalog associations path.
    """

    _target_: str = "otg.gwas_catalog.GWASCatalogStep"
    catalog_studies_file: str = MISSING
    catalog_ancestry_file: str = MISSING
    catalog_sumstats_lut: str = MISSING
    catalog_associations_file: str = MISSING
    variant_annotation_path: str = MISSING
    min_r2: float = 0.5
    ld_matrix_template: str = MISSING
    ld_index_template: str = MISSING
    ld_populations: List[str] = field(
        default_factory=lambda: [
            "afr",  # African-American
            "amr",  # American Admixed/Latino
            "asj",  # Ashkenazi Jewish
            "eas",  # East Asian
            "fin",  # Finnish
            "nfe",  # Non-Finnish European
            "nwe",  # Northwestern European
            "seu",  # Southeastern European
        ]
    )
    catalog_studies_out: str = MISSING
    catalog_associations_out: str = MISSING


@dataclass
class GeneIndexStepConfig:
    """Gene index step requirements.

    Attributes:
        target_path (str): Open targets Platform target dataset path.
        gene_index_path (str): Output gene index path.
    """

    _target_: str = "otg.gene_index.GeneIndexStep"
    target_path: str = MISSING
    gene_index_path: str = MISSING


@dataclass
class GWASCatalogSumstatsPreprocessConfig:
    """GWAS Catalog Sumstats Preprocessing step requirements.

    Attributes:
        raw_sumstats_path (str): Input raw GWAS Catalog summary statistics path.
        out_sumstats_path (str): Output GWAS Catalog summary statistics path.
        study_id (str): GWAS Catalog study identifier.
    """

    _target_: str = (
        "otg.gwas_catalog_sumstat_preprocess.GWASCatalogSumstatsPreprocessStep"
    )
    raw_sumstats_path: str = MISSING
    out_sumstats_path: str = MISSING
    study_id: str = MISSING


@dataclass
class FinnGenStepConfig:
    """FinnGen study table ingestion step requirements.

    Attributes:
        finngen_phenotype_table_url (str): FinnGen API for fetching the list of studies.
        finngen_release_prefix (str): Release prefix pattern.
        finngen_sumstat_url_prefix (str): URL prefix for summary statistics location.
        finngen_sumstat_url_suffix (str): URL prefix suffix for summary statistics location.
        finngen_study_index_out (str): Output path for the FinnGen study index dataset.
    """

    _target_: str = "otg.finngen.FinnGenStep"
    finngen_phenotype_table_url: str = MISSING
    finngen_release_prefix: str = MISSING
    finngen_sumstat_url_prefix: str = MISSING
    finngen_sumstat_url_suffix: str = MISSING
    finngen_study_index_out: str = MISSING


@dataclass
class UKBiobankStepConfig:
    """UKBiobank study table ingestion step requirements.

    Attributes:
        ukbiobank_manifest (str): UKBiobank manifest of studies.
        ukbiobank_study_index_out (str): Output path for the UKBiobank study index dataset.
    """

    _target_: str = "otg.ukbiobank.UKBiobankStep"
    ukbiobank_manifest: str = MISSING
    ukbiobank_study_index_out: str = MISSING


@dataclass
class UKBiobankSumstatsPreprocessConfig:
    """UKBiobank Sumstats Preprocessing step requirements.

    Attributes:
        ukbiobank_sumstat_path (str): UKBiobank input file path for sumstats location.
        ukbiobank_study_id (str): UKBiobank study identifier.
        ukbiobank_sumstats_out (str): Output UKBiobank summary statistics path.
        ukbiobank_sumstats_dirs (list): List of directories containing the summary statistics files.
        ukbiobank_file_batch_size (int): Number of files to process in each iteration.
        merged_ukbiobank_sumstats_out (str): Output file path for the merged summary statistics.
    """

    _target_: str = "otg.ukbiobank_sumstat_preprocess.UKBiobankSumstatsPreprocessStep"
    ukbiobank_sumstat_path: str = MISSING
    #ukbiobank_study_id: str = MISSING
    ukbiobank_sumstats_out: str = MISSING
    ukbiobank_sumstats_dirs: List[str] = MISSING
    ukbiobank_file_batch_size: int = MISSING
    ukbiobank_sumstats_merged_out: str = MISSING

@dataclass
class FinemappingConfig:
    """Finemapping step requirements.

    Attributes:
        target (str): Target locus name.
        target_chrom (int): Target chromosome number.
        target_pos (int): Lead SNP position.
        lead_snp_ID (str): Lead SNP variant ID (must match ID format provided in LD matrix).
        n_sample (int): Study sample size.
        outlier_method (str): Outerlier detection method.
        r2_threshold (int): LD threshold for DENTIST outlier detection method.
        nlog10p_dentist_s_threshold (int): = DENTIST measure threshold for outlier detection.
        sumstat_path (str): Input file path for sumstats location.
        study_id (str): UKBiobank study identifier.
        ukbiobank_sumstats_out (str): Output UKBiobank summary statistics path.
        ukbiobank_sumstats_dirs (list): List of directories containing the summary statistics files.
    """

    _target_: str = "otg.pyspark_finmapping_pipeline.FinemappingStep"
    target: str = MISSING
    target_chrom: int = MISSING
    target_pos: int = MISSING
    lead_snp_ID: str = MISSING
    n_sample: int = MISSING
    outlier_method: str = MISSING
    r2_threshold: int = MISSING
    nlog10p_dentist_s_threshold: int = MISSING

# Register all configs
def register_configs() -> None:
    """Register configs."""
    cs = ConfigStore.instance()
    cs.store(name="config", node=Config)
    cs.store(name="session_config", group="session", node=SessionConfig)
    cs.store(name="gene_index", group="step", node=GeneIndexStepConfig)
    cs.store(name="ld_index", group="step", node=LDIndexStepConfig)
    cs.store(name="variant_index", group="step", node=VariantIndexStepConfig)
    cs.store(name="variant_annotation", group="step", node=VariantAnnotationStepConfig)
    cs.store(name="v2g", group="step", node=V2GStepConfig)
    cs.store(name="colocalisation", group="step", node=ColocalisationStepConfig)
    cs.store(name="gwas_catalog", group="step", node=GWASCatalogStepConfig)
    cs.store(name="finngen", group="step", node=FinnGenStepConfig)
    cs.store(name="ukbiobank", group="step", node=UKBiobankStepConfig)
    cs.store(
        name="gwas_catalog_sumstats_preprocess",
        group="step",
        node=GWASCatalogSumstatsPreprocessConfig,
    )
    cs.store(
        name="ukbiobank_sumstat_preprocess",
        group="step",
        node=UKBiobankSumstatsPreprocessConfig,
    )
    cs.store(
        name="finemapping",
        group="step",
        node=FinemappingConfig,
    )

