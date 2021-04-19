import pytest

from subprocess import run


@pytest.mark.parametrize(
    "cmd_line",
    [
        (["bin/tadrep"]),  # not provided
        (["bin/tadrep", '--genome']),  # missing path
        (["bin/tadrep", '--genome', '']),  # empty
        (["bin/tadrep", '--genome', 'foo.fasta']),  # not existing
        (["bin/tadrep", '--genome', 'foo.fasta', '']),  # not existing, empty
        (["bin/tadrep", '--genome', '', 'foo.fasta']),  # empty, not existing
    ]
)
def test_genome_failing(cmd_line):
    # test genome arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    "cmd_line",
    [
        (["bin/tadrep"]),  # not provided
        (["bin/tadrep", '--plasmids']),  # missing path
        (["bin/tadrep", '--plasmids', '']),  # empty
        (["bin/tadrep", '--plasmids', 'foo.fasta']),  # not existing
    ]
)
def test_plasmids_failing(cmd_line):
    # test plasmid arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.slow
def test_full(tmpdir):
    # all parameter OK
    proc = run(["bin/tadrep", '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir])
    assert proc.returncode == 0