import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from aa_explorer.aa_data import AMINO_ACIDS, get_property
from aa_explorer.profiler import compare_point_mutation, compute_profile
from aa_explorer.sequence_io import clean_sequence, validate_sequence


def test_all_20_amino_acids_present():
    assert len(AMINO_ACIDS) == 20


def test_get_property():
    assert get_property("W", "hydropathy") == -0.9
    assert get_property("i", "hydropathy") == 4.5  # lowercase handled


def test_get_property_unknown_raises():
    try:
        get_property("X", "hydropathy")
        raise AssertionError("expected KeyError")
    except KeyError:
        pass


def test_clean_sequence_strips_whitespace_and_uppercases():
    assert clean_sequence("  mkt\nayi ak\n") == "MKTAYIAK"


def test_validate_sequence_flags_invalid_chars():
    assert validate_sequence("MKTX*AYI") == ["*", "X"]
    assert validate_sequence("MKTAYI") == []


def test_compute_profile_basic():
    profile = compute_profile("MKTAYIAK", label="test", window_size=3)
    assert len(profile.per_residue["hydropathy"]) == 8
    assert len(profile.smoothed["hydropathy"]) == 8
    assert set(profile.composition_mean.keys()) == {
        "hydropathy", "volume", "charge", "aromaticity", "mol_weight"
    }


def test_compute_profile_rejects_invalid_sequence():
    try:
        compute_profile("MKTXAYIAK")
        raise AssertionError("expected ValueError")
    except ValueError:
        pass


def test_compute_profile_rejects_even_window():
    try:
        compute_profile("MKTAYIAK", window_size=4)
        raise AssertionError("expected ValueError")
    except ValueError:
        pass


def test_compare_point_mutation():
    result = compare_point_mutation("MKTAYIAK", position=1, mutant_residue="A")
    assert result["wt_residue"] == "M"
    assert result["mutant_residue"] == "A"
    assert result["label"] == "M1A"
    assert result["deltas"]["hydropathy"]["wt"] == get_property("M", "hydropathy")
    assert result["deltas"]["hydropathy"]["mutant"] == get_property("A", "hydropathy")


def test_compare_point_mutation_out_of_range():
    try:
        compare_point_mutation("MKT", position=10, mutant_residue="A")
        raise AssertionError("expected ValueError")
    except ValueError:
        pass


if __name__ == "__main__":
    test_all_20_amino_acids_present()
    test_get_property()
    test_get_property_unknown_raises()
    test_clean_sequence_strips_whitespace_and_uppercases()
    test_validate_sequence_flags_invalid_chars()
    test_compute_profile_basic()
    test_compute_profile_rejects_invalid_sequence()
    test_compute_profile_rejects_even_window()
    test_compare_point_mutation()
    test_compare_point_mutation_out_of_range()
    print("All tests passed.")
