import bz2
import gzip
import json
from pathlib import Path
from typing import Union
from multiprocessing import Pool, cpu_count
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqIO import parse
from itertools import repeat


def get_nt_with_begin_end(sequence, forward_seq, forward_end):
    """Get the nt library region that is stripped at primer

    Now it may not be the codin region because it may have begin and end barcodes that are not in coding region
    """
    return sequence[str(sequence).find(forward_seq) + len(forward_seq) : str(sequence).find(forward_end)]


def get_coding_nt(sequence, first_codon, end_codon):
    """Strip the begin and end down just the coding region to be translated"""
    return sequence[first_codon:end_codon]


class ReferenceLibrary:
    def __init__(self, forward_primer: str, reverse_primer: str, reference_nt: str, first_codon=None, end_codon=None):
        self.forward_primer = forward_primer
        self.end_primer = reverse_primer
        self.first_codon = first_codon
        self.end_codon = end_codon
        self.reference_nt = reference_nt
        self.stripped_ref = get_nt_with_begin_end(reference_nt, forward_primer, reverse_primer)
        self.reading_frame_ref = get_coding_nt(self.stripped_ref, self.first_codon, self.end_codon)
        self.reading_frame_ref_aa = str(Seq(self.reading_frame_ref).translate())
        self.skipped = 0
        self.ref_library = []

    def apply(self, fasta_file: Union[str, Path]) -> pd.DataFrame:
        if fasta_file.endswith("bz2"):
            file_handle = bz2.open(fasta_file, "rt")
        else:
            file_handle = open(fasta_file, "r")

        for sequence in parse(file_handle, "fastq"):
            sequence_with_begin_and_end = get_nt_with_begin_end(sequence.seq, self.forward_primer, self.end_primer)
            sequence_coding = get_coding_nt(sequence_with_begin_and_end, self.first_codon, self.end_codon)
            if len(sequence_coding) != len(self.reading_frame_ref):
                self.skipped += 1
                continue
            sequence_coding_aa = sequence_coding.translate()
            if "*" in sequence_coding_aa:
                self.skipped += 1
                continue
            for aa_index, (wt_seq, seq) in enumerate(zip(self.reading_frame_ref_aa, sequence_coding_aa)):
                # only want to handle mutations for our reference library
                if str(wt_seq) != str(seq):

                    # 3 x aa index = codon in
                    codon_index = aa_index * 3

                    # first handle upstream barcode
                    if aa_index == 0:
                        # if we are at codon 0, we can only go 1 way
                        upstream_barcode = None
                        us_reference_barcdoe = None
                    elif aa_index == 1:
                        # we are at codon position 1 meaing we only have 3 nt to work with
                        us_reference_barcdoe = self.reading_frame_ref[0:codon_index]
                    else:
                        # want to reach 6 back to the codon of interest
                        upstream_barcode = sequence_coding[codon_index - 6 : codon_index]
                        us_reference_barcdoe = self.reading_frame_ref[codon_index - 6 : codon_index]

                    # downstream is easier since it will automatically put in "" when we stretch past the end
                    downstream_barcode = sequence_coding[codon_index + 3 : codon_index + 9]
                    ds_reference_barcdoe = self.reading_frame_ref[codon_index + 3 : codon_index + 9]

                    # 3 letter codon
                    codon = sequence_coding[codon_index : codon_index + 3]

                    # we want to make sure that whatever barcode we are seeing is not the same as the reference barcode stretch
                    if (downstream_barcode == str(ds_reference_barcdoe)) and (
                        upstream_barcode == str(us_reference_barcdoe)
                    ):
                        continue

                    self.ref_library.append(
                        {
                            "codon_index": codon_index,
                            "aa_index": aa_index,
                            "wt_aa": wt_seq,
                            "mut_aa": seq,
                            "codon": str(codon),
                            "upstream_barcode": str(upstream_barcode),
                            "downstream_barcode": str(downstream_barcode),
                        }
                    )

        # this is our pre-normalized library, i.e before we have found the most common barcode per postion
        self.ref_library = pd.DataFrame(self.ref_library)

        # now we can groupby and sort by counts of paired barcodes to find the most frequent barcode and call that the correct one
        self.ref_library_normal = (
            pd.DataFrame(self.ref_library)
            .groupby(["codon_index", "aa_index", "wt_aa", "upstream_barcode", "downstream_barcode"])
            .count()
            .rename({"mut_aa": "count"}, axis=1)
            .drop("codon", axis=1)
            .reset_index()
            .sort_values("count")[::-1]
            .groupby(["codon_index", "aa_index", "wt_aa"])
            .head(1)
            .sort_values("aa_index")
            .replace("", "None")
            .reset_index(drop=True)
        )
        return self.ref_library_normal

    def to_json(self, filename: str) -> None:
        """Save the reference library to a json file if file name ends with .gz will compress to gzip"""
        _json = {
            "normal": self.ref_library_normal.to_json(orient="records", indent=4),
            "raw": self.ref_library.to_json(orient="records", indent=4),
            "forward_primer": self.forward_primer,
            "reverse_primer": self.end_primer,
            "first_codon": self.first_codon,
            "last_codon": self.end_codon,
            "reference_nt": self.reference_nt,
        }
        if filename.endswith(".gz"):
            with gzip.open(filename, "wt") as f:
                json.dump(_json, f)
        else:
            json.dump(_json, gzip.open(filename, "w"), indent=4)

    @staticmethod
    def from_json(filename: str) -> "ReferenceLibrary":
        """Load a reference library from a json file"""
        if filename.endswith(".gz"):
            with gzip.open(filename, "rt") as f:
                _json = json.load(f)
        _ref = ReferenceLibrary(
            _json["forward_primer"],
            _json["reverse_primer"],
            _json["reference_nt"],
            _json["first_codon"],
            _json["last_codon"],
        )
        _ref.ref_library_normal = pd.read_json(_json["normal"], orient="records", dtype={"codon_index": int})
        _ref.ref_library = pd.read_json(_json["raw"], orient="records", dtype={"codon_index": int})
        return _ref


def _multi(
    sequence: str,
    reference_seq: str,
    reference_seq_aa: str,
    forward_primer: str,
    end_primer: str,
    first_codon: int,
    end_codon: int,
    lookup_hash: dict,
):
    single_seq_library = []
    sequence_with_begin_and_end = get_nt_with_begin_end(sequence.seq, forward_primer, end_primer)
    sequence_coding = get_coding_nt(sequence_with_begin_and_end, first_codon, end_codon)

    # if its not the same lengh as ref sequence, continue
    if len(sequence_coding) != len(reference_seq):
        return

    # if it has stop codon continue
    sequence_coding_aa = str(Seq(sequence_coding).translate())
    if "*" in sequence_coding_aa:
        return

    for aa_index, (wt_seq, seq) in enumerate(zip(reference_seq_aa, sequence_coding_aa)):
        codon_index = aa_index * 3
        if aa_index == 0:
            upstream_barcode = None
        elif aa_index == 1:
            upstream_barcode = sequence_coding[0:codon_index]
        else:
            upstream_barcode = sequence_coding[codon_index - 6 : codon_index]
        downstream_barcode = sequence_coding[codon_index + 3 : codon_index + 9]
        codon = sequence_coding[codon_index : codon_index + 3]

        lookup = lookup_hash[codon_index]
        lookup_us_barcode = lookup["upstream_barcode"]
        lookup_ds_barcode = lookup["downstream_barcode"]
        break_us = False
        break_ds = False
        if lookup_us_barcode == "None" or (lookup_us_barcode == upstream_barcode):
            break_us = True

        if lookup_ds_barcode == "None" or (lookup_ds_barcode == downstream_barcode):
            break_ds = True

        if break_us and break_ds:
            single_seq_library.append(
                {
                    "codon_index": codon_index,
                    "aa_index": aa_index,
                    "wt_aa": wt_seq,
                    "mut_aa": seq,
                    "codon": str(codon),
                    "upstream_barcode": str(upstream_barcode),
                    "downstream_barcode": str(downstream_barcode),
                    "sequence_id": sequence.id,
                }
            )
    return pd.DataFrame(single_seq_library)


def analyze_library(ref: ReferenceLibrary, fastq_file: Union[Path, str], multi: bool = True) -> pd.DataFrame:
    """Given a fastq file, count and normalize nnk occurences relative to reference library"""
    if not isinstance(ref, ReferenceLibrary):
        raise ValueError("ref must be a ReferenceLibrary")
    if isinstance(fastq_file, str):
        fastq_file = Path(fastq_file)

    lookup_hash = ref.ref_library_normal.set_index("codon_index").to_dict(orient="index")

    if fastq_file.suffix == ".bz2":
        fastq_file = bz2.open(fastq_file, "rt")
    else:
        fastq_file = open(fastq_file, "r")
    sequences = list(parse(fastq_file, "fastq"))
    if multi:
        with Pool(cpu_count()) as pool:
            results = pool.starmap(
                _multi,
                zip(
                    sequences,
                    repeat(ref.reading_frame_ref),
                    repeat(ref.reading_frame_ref_aa),
                    repeat(ref.forward_primer),
                    repeat(ref.end_primer),
                    repeat(ref.first_codon),
                    repeat(ref.end_codon),
                    repeat(lookup_hash),
                ),
            )

        first_df = pd.concat(results)
    else:
        first_df = pd.concat(map(_multi, sequences))
    library_df = (
        pd.DataFrame(first_df)
        .groupby(["aa_index", "mut_aa", "wt_aa"])
        .count()["codon_index"]
        .reset_index()
        .rename({"codon_index": "count"}, axis=1)
    )
    naive_lib = (
        (library_df.groupby(["aa_index", "mut_aa", "wt_aa"]).apply(lambda x: x["count"] / library_df["count"].sum()))
        .to_frame()
        .rename({"count": "freq"}, axis=1)
    )
    naive_lib.index = naive_lib.index.droplevel(3)
    counts_and_freq_df = library_df.groupby(["aa_index", "mut_aa", "wt_aa"]).sum().join(naive_lib).reset_index()
    return counts_and_freq_df
