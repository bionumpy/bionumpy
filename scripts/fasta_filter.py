from bionumpy import bnp_open
import logging

logging.basicConfig(level="DEBUG")

def main(in_filename, out_filename):
    of = bnp_open(out_filename, "w")
    for sequence_entries in bnp_open(in_filename):
        mask = (sequence_entries.sequence == ord("C")).sum(axis=-1) > 50
        of.write(sequence_entries[mask])

if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])
