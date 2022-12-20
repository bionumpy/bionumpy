import time
t00 = time.perf_counter()
import pyranges
import sys
#from bionumpy.util.cli import run_as_commandline


#@run_as_commandline
def run(a: str, b: str, outfile: str):

    t0 = time.perf_counter()
    a = pyranges.read_bed(a)
    b = pyranges.read_bed(b)
    print("Reading took", (time.perf_counter()-t0))


    t0 = time.perf_counter()
    intersection = a.intersect(b, how="last")
    print("Intersection took", (time.perf_counter()-t0))


    t0 = time.perf_counter()
    intersection.to_bed(outfile)
    print("Writing to bed took", (time.perf_counter()-t0))
    print("Total time", (time.perf_counter()-t00))


if __name__ == "__main__":
    run(*sys.argv[1:])
