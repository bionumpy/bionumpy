
chrom_sizes = bnp.open(input.chrom_sizes).read()
global_offset = GlobalOffset(chrom_sizes)
a = bnp.open(input[0], buffer_type=bnp.Bed6Buffer).read()
a = global_offset.from_local_interval(a, do_clip=True)
# a = bnp.open(input[0], buffer_type=bnp.Bed6Buffer).read_chunks()
b = bnp.open(input[0], buffer_type=bnp.Bed6Buffer).read()
b = global_offset.from_local_interval(b, do_clip=True)
# b = bnp.open(input[1], buffer_type=bnp.Bed6Buffer).read_chunks()
# ms = bnp.streams.MultiStream(chrom_sizes, a=a, b=b)
result = bnp.arithmetics.intersect(a, b)
result = global_offset.to_local_interval(result)
with bnp.open(output[0], "w", buffer_type=bnp.Bed6Buffer) as f:
    f.write(result)
