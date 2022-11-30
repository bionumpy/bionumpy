from ..streams import streamable
from ..datatypes import BamEntry, Bed6
from ..encoded_array import EncodedArray


@streamable()
def alignment_to_interval(alignment: BamEntry) -> Bed6:
    """Get the stranded interval that an alignment covers on the reference contig

    Parameters
    ----------
    alignment : BamEntry
        Aligments

    Returns
    -------
    Bed6
        Corresponding intervals

    """
    strand = alignment.flag & np.uint16(16)
    strand = EncodedArray(np.where(strand, ord("-"), ord("+"))[:, None])
    length = count_reference_length(alignment.cigar_op, alignment.cigar_length)
    return Bed6(alignment.chromosome,
                alignment.position,
                alignment.position+length,
                alignment.name,
                alignment.mapq,
                strand)
