"""Test finemapping precossing dataset."""
from __future__ import annotations

from otg.dataset.study_locus import StudyLocus
from otg.dataset.summary_statistics import SummaryStatistics
from otg.dataset.fimeapping import LDMatrix
from otg.dataset.fimeapping import Outliers


def test_summary_statistics__creation(
    mock_summary_statistics: SummaryStatistics,
) -> None:
    """Test gene index creation with mock gene index."""
    assert isinstance(mock_summary_statistics, SummaryStatistics)

def test_ld_matrix_creation(
    mock_ld_matrix: LDMatrix,
) -> None:
    """Test gene index creation with mock gene index."""
    assert isinstance(mock_ld_matrix, SummaryStatistics)

def test_snp_matching(
    mock_summary_statistics_filtered: SummaryStatisticsFiltered,
    mock_ld_matrix_filtered: LDMatrixFiltered,
) -> None:
    """Test gene index creation with mock gene index."""
    assert isinstance()

def test_allele_flip(
    mock_summary_statistics_filtered: SummaryStatisticsFiltered,
    mock_ld_matrix_filtered: LDMatrixFiltered,
) -> None:
    """Test gene index creation with mock gene index."""
    assert isinstance()

def test_outliers_detection(
    mock_outliers: Outliers,
) -> None:
    """Test gene index creation with mock gene index."""
    assert isinstance(mock_outliers, SummaryStatistics)

def test_finemapping_results(
    mock_finemapping_results: FMResults,
) -> None:
    """Test gene index creation with mock gene index."""
    assert isinstance(mock_finemapping_results, FMResults)
