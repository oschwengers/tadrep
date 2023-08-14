from pathlib import Path
from subprocess import run

import pytest

from .conftest import FILES

def test_full(tmpdir):
    # all parameter OK
    proc = run(['bin/tadrep', 'detect', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir])
    
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
        assert output_path.stat().st_size > 0