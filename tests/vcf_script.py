import bionumpy as bnp

def method_name(outfile, infile):
    with bnp.open(outfile, 'w') as f:
        for chunk in bnp.open(infile).read_chunks():
            f.write(chunk[(chunk.info.LoF.lengths > 0) & chunk.info.ONCOGENE])


if __name__ == "__main__":
    import sys
    method_name(sys.argv[1], sys.argv[2])

