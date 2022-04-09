from .file_buffers import *
from .parser import *
from .delimited_buffers import *
from .chromosome_provider import *

def get_vcf_file(file_obj, mode=None):
    if mode is None:
        mode = "chromosome_stream"
    assert mode in ("chromosome_stream", "dict", "stream"), mode
    if mode=="chromosome_stream":
        stream = (buf.get_variants() for buf in BufferedNumpyParser(file_obj, VCFBuffer).get_chunks())
        return ChromosomeStreamProvider(stream)

def bnp_open(filename, mode=None):
    if filename.endswith("vcf"):
        return get_vcf_file(open(filename, "rb"), mode)

    
