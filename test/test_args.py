import pytest

from subprocess import run


@pytest.mark.parametrize(
    "parameters",
    [
        ([]),  # not provided
        (['--genome']),  # missing path
        (['--genome', '']),  # empty
        (['--genome', 'foo.fasta']),  # not existing
        (['--genome', 'foo.fasta', '']),  # not existing, empty
        (['--genome', '', 'foo.fasta']),  # empty, not existing
        (['--genome', 'foo.fasta', '', 'test/data/draft.fna']),  # not existing, empty, OK
        (['--genome', '', 'foo.fasta', 'test/data/draft.fna'])  # empty, not existing, OK
    ]
)
def test_genome_failing(parameters, tmpdir):
    # test genome arguments
    cmd_line = ['bin/tadrep', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir] + parameters
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    "parameters",
    [
        ([]),  # not provided
        (['--plasmids']),  # missing path
        (['--plasmids', '']),  # empty
        (['--plasmids', 'foo.fasta']),  # not existing
    ]
)
def test_plasmids_failing(parameters, tmpdir):
    # test genome arguments
    cmd_line = ['bin/tadrep', '--genome', 'test/data/draft.fna', '--output', tmpdir] + parameters
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'arguments',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo']),  # string
        (['-1']),  # smaller than zero
        (['1.1']),  # float
        (['101'])  # larger than 100 %
    ]
)
def test_min_contig_coverage_failing(arguments, tmpdir):
    # test min-contig-coverage arguments
    cmd_line = ['bin/tadrep', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir, '--min-contig-coverage'] + arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'arguments',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo']),  # string
        (['-1']),  # smaller than zero
        (['1.1']),  # float
        (['101'])  # larger than 100 %
    ]
)
def test_min_contig_identity_failing(arguments, tmpdir):
    # test min-contig-identity arguments
    cmd_line = ['bin/tadrep', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir, '--min-contig-identity'] + arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'arguments',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo']),  # string
        (['-1']),  # smaller than zero
        (['1.1']),  # float
        (['101'])  # larger than 100 %
    ]
)
def test_min_plasmid_coverage_failing(arguments, tmpdir):
    # test min-plasmid-coverage arguments
    cmd_line = ['bin/tadrep', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir, '--min-plasmid-coverage'] + arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'arguments',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo']),  # string
        (['-1']),  # smaller than zero
        (['1.1']),  # float
        (['101'])  # larger than 100 %
    ]
)
def test_min_plasmid_identity_failing(arguments, tmpdir):
    # test min-plasmid-identity arguments
    cmd_line = ['bin/tadrep', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir, '--min-plasmid-identity'] + arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'arguments',
    [
        ([]),  # not provided
        (['']),  # empty
        (['foo']),  # string
        (['-1']),  # smaller than zero
        (['1.1']),  # float
    ]
)
def test_gap_sequence_length_failing(arguments, tmpdir):
    # test gap-sequence-length arguments
    cmd_line = ['bin/tadrep', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir, '--gap-sequence-length'] + arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    "arguments",
    [
        ([]),  # missing path
        (['']),  # empty
    ]
)
def test_prefix_failing(arguments, tmpdir):
    # test writable output dir
    cmd_line = ['bin/tadrep', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir, '--prefix'] + arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


def test_output_failing():
    # test writable output dir
    cmd_line = ['bin/tadrep', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', '/']
    proc = run(cmd_line)
    assert proc.returncode != 0


def test_tmpdir_failing(tmpdir):
    # test writable tmp dir
    cmd_line = ['bin/tadrep', '--genome', 'test/data/draft.fna', '--plasmids', 'test/data/plasmids.fna', '--output', tmpdir, '--tmpdir', '/']
    proc = run(cmd_line)
    assert proc.returncode != 0
