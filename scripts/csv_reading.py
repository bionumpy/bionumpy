from importlib import reload
from npstructures import npdataclass
from bionumpy.delimited_buffers import get_bufferclass_for_datatype
import bionumpy as bnp
reload(bnp)


@npdataclass
class VCFAsTSV:
    chromosome: str
    position: int
    id: None
    ref: str
    alt: str


print(bnp.open("example_data/variants.vcf",
               mode="full",
               buffer_type=get_bufferclass_for_datatype(VCFAsTSV)))
