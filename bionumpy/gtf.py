from .datatypes import GFFEntry
from .bnpdataclass import bnpdataclass
from .strops import str_equal, split
from . import streamable


@bnpdataclass
class GFFGeneEntry(GFFEntry):
    gene_id: str


@bnpdataclass
class GFFTranscriptEntry(GFFGeneEntry):
    transcript_id: str


@bnpdataclass
class GFFExonEntry(GFFGeneEntry):
    transcript_id: str


def get_attributes(gtf_entries, attribute_names):
    gtf_entries.atributes[:, -1] = " "
    all_features = split(gtf_entries.atributes.ravel(), " ")
    keys = all_features[:-1:2]
    print(keys)
    values = all_features[1::2, 1:-2]
    return {name: values[str_equal(keys, name)] for name in attribute_names}


@streamable()
def get_genes(gtf_entries):
    genes = gtf_entries[str_equal(gtf_entries.feature_type, "gene")]
    attributes = get_attributes(genes, ["gene_id"])
    for key, val in attributes.items():
        assert len(val) == len(genes), (len(val), len(genes))
    return GFFGeneEntry(*genes.shallow_tuple(), **attributes)


def get_transcripts(gtf_entries):
    genes = gtf_entries[str_equal(gtf_entries.feature_type, "transcript")]
    attributes = get_attributes(genes, ["transcript_id", "gene_id"])
    return GFFTranscriptEntry(*genes.shallow_tuple(), **attributes)


def get_exons(gtf_entries):
    genes = gtf_entries[str_equal(gtf_entries.feature_type, "exon")]
    attributes = get_attributes(genes, ["transcript_id", "gene_id", "exon_id"])
    return GFFExonEntry(*genes.shallow_tuple(), **attributes)
