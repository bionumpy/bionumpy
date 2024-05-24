from ..bnpdataclass import bnpdataclass, BNPDataClass
from ..encodings import StrandEncoding
from ..io.strops import str_equal, split, join
from ..io.regexp import match_regexp, match_regexp_string_array
from ..string_array import as_string_array
from ..typing import SequenceID


@bnpdataclass
class GTFEntry:
    chromosome: SequenceID
    source: str
    feature_type: SequenceID
    start: int
    stop: int
    score: str
    strand: StrandEncoding
    phase: str
    atributes: str

    def _get_attributes(self, attribute_names):
        group_thing = r''' \"(.*?)\"'''
        return {name: match_regexp_string_array(self.atributes.ravel(), name + group_thing)
                for name in attribute_names}
        self.atributes[:, -1] = " "
        all_features = split(self.atributes.ravel(), " ")
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
        transcripts = gtf_entries[str_equal(gtf_entries.feature_type, "transcript")]
        attributes = transcripts._get_attributes(["transcript_id", "gene_id"])
        return GFFTranscriptEntry(*transcripts.shallow_tuple(), **attributes)
    
    def get_exons(gtf_entries):
        genes = gtf_entries[str_equal(gtf_entries.feature_type, "exon")]
        attributes = genes._get_attributes(["transcript_id", "gene_id", "exon_id"])
        return GFFExonEntry(*genes.shallow_tuple(), **attributes)


class GFFEntry(GTFEntry):
    def _get_attributes(self, attribute_names):
        joined = join(self.atributes, ';')
        all_features = split(join(self.atributes, ';'), [";", '='])
        keys = all_features[:-1:2]
        values = all_features[1::2]
        return {name: as_string_array(values[str_equal(keys, name)]) for name in attribute_names}


@bnpdataclass
class GFFGeneEntry(GFFEntry):
    gene_id: SequenceID


@bnpdataclass
class GFFTranscriptEntry(GFFGeneEntry):
    transcript_id: SequenceID


@bnpdataclass
class GFFExonEntry(GFFTranscriptEntry):
    exon_id: SequenceID
