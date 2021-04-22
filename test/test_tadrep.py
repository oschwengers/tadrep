import pytest

from subprocess import run


@pytest.mark.slow
def test_full(tmpdir):
    # all parameter OK
    proc = run(['bin/tadrep', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir])
    assert proc.returncode == 0