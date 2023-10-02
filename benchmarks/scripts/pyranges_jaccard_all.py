import pyranges
files = {i: pyranges.read_bed(i) for i in snakemake.input}
results = {}
n = 0
for a in files:
    for b in files:
        run_id = frozenset([a, b])
        if a == b or run_id in results:
            continue
        j = files[a].stats.jaccard(files[b])
        results[run_id] = j
        n += 1

print("%d pairs compared" % n)

with open(snakemake.output[0], "w") as f:
    f.write(str(list(results.values())) + "\n")
