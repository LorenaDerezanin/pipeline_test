import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_trim_PE_reads():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/trim_PE_reads/data")
        expected_path = PurePosixPath(".tests/unit/trim_PE_reads/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("results/trimmed_reads/sample.R1.paired_val_1.fq.gz results/trimmed_reads/sample.R1.paired.fq.gz_trimming_report.txt results/trimmed_reads/sample.R2.paired_val_2.fq.gz results/trimmed_reads/sample.R2.paired.fq.gz_trimming_report.txt", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/trimmed_reads/sample.R1.paired_val_1.fq.gz results/trimmed_reads/sample.R1.paired.fq.gz_trimming_report.txt results/trimmed_reads/sample.R2.paired_val_2.fq.gz results/trimmed_reads/sample.R2.paired.fq.gz_trimming_report.txt",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
