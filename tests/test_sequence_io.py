import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from aa_explorer.sequence_io import read_fasta

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"


def test_read_example_fasta():
    records = read_fasta(EXAMPLES_DIR / "example.fasta")
    assert "lysozyme_c_hen_egg_white" in records
    assert "myoglobin_sperm_whale" in records
    assert records["lysozyme_c_hen_egg_white"].startswith("KVFGRCELAAAMKRHG")


def test_read_fasta_multiline_sequence():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "test.fasta"
        p.write_text(">seq1\nMKTA\nYIAK\n>seq2\nGGGG\n")
        records = read_fasta(p)
        assert records["seq1"] == "MKTAYIAK"
        assert records["seq2"] == "GGGG"


def test_read_fasta_empty_raises():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "empty.fasta"
        p.write_text("")
        try:
            read_fasta(p)
            raise AssertionError("expected ValueError")
        except ValueError:
            pass


if __name__ == "__main__":
    test_read_example_fasta()
    test_read_fasta_multiline_sequence()
    test_read_fasta_empty_raises()
    print("All FASTA tests passed.")
