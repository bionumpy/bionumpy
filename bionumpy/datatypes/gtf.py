from ..bnpdataclass import bnpdataclass
from ..encodings import StrandEncoding
from ..io.strops import str_equal, split


@bnpdataclass
class GFFEntry:
    chromosome: str
    source: str
    feature_type: str
    start: int
    stop: int
    score: str
    strand: StrandEncoding
    phase: str
    atributes: str

    def _get_attributes(gtf_entries, attribute_names):
        gtf_entries.atributes[:, -1] = " "
        all_features = split(gtf_entries.atributes.ravel(), " ")
        keys = all_features[:-1:2]
        values = all_features[1::2, 1:-2]
        return {name: values[str_equal(keys, name)] for name in attribute_names}

    def get_genes(gtf_entries):
        genes = gtf_entries[str_equal(gtf_entries.feature_type, "gene")]
        attributes = genes._get_attributes(["gene_id"])
        for key, val in attributes.items():
            assert len(val) == len(genes), (len(val), len(genes))
        return GFFGeneEntry(*genes.shallow_tuple(), **attributes)

    def get_transcripts(gtf_entries):
        genes = gtf_entries[str_equal(gtf_entries.feature_type, "transcript")]
        attributes = genes._get_attributes(["transcript_id", "gene_id"])
        return GFFTranscriptEntry(*genes.shallow_tuple(), **attributes)
    
    def get_exons(gtf_entries):
        genes = gtf_entries[str_equal(gtf_entries.feature_type, "exon")]
        attributes = genes._get_attributes(["transcript_id", "gene_id", "exon_id"])
        return GFFExonEntry(*genes.shallow_tuple(), **attributes)


@bnpdataclass
class GFFGeneEntry(GFFEntry):
    gene_id: str


@bnpdataclass
class GFFTranscriptEntry(GFFGeneEntry):
    transcript_id: str


@bnpdataclass
class GFFExonEntry(GFFTranscriptEntry):
    exon_id: str
