import shutil as sh

from pathlib import Path
from subprocess import run

from .conftest import FILES

def test_full(tmpdir):
    tmpdir_path = Path(tmpdir).resolve()
    sh.copyfile(Path('test/data/db.json'), tmpdir_path.joinpath('db.json'))
    
    # all parameter OK
    proc = run(['bin/tadrep', '-v', '--output', tmpdir_path, 'detect', '--genome', 'test/data/draft.fna'])
    assert proc.returncode == 0
    
    for file in FILES:  # all files in place and not empty
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
        assert output_path.stat().st_size > 0