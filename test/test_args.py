import os
from pathlib import Path
from subprocess import run


def test_genome():
    # test genome arguments

    # not provided
    proc = run(["bin/tadrep"])
    assert proc.returncode != 0

    # parameter missing path
    proc = run(["bin/tadrep", '--genome'])
    assert proc.returncode != 0

    # non-existing
    proc = run(["bin/tadrep", '--genome', 'foo.fasta'])
    assert proc.returncode != 0


def test_plasmids():
    # test database arguments

    # not provided
    proc = run(["bin/tadrep"])
    assert proc.returncode != 0

    # parameter missing path
    proc = run(["bin/tadrep", '--plasmids'])
    assert proc.returncode != 0

    # non-existing
    proc = run(["bin/tadrep", '--plasmids', 'foo.fasta'])
    assert proc.returncode != 0


def test_full(tmpdir):
    # all parameter OK
    proc = run(["bin/tadrep", '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir])
    assert proc.returncode == 0