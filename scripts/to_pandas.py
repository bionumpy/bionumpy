import bionumpy as bnp
from tests.util import get_file_name


def main():
	f = bnp.open(get_file_name("example_data/big.fq.gz"))
	for chunk in f.read_chunks():
		print(chunk)
		df = chunk.topandas()
		print(df)


if __name__ == "__main__":
    main()
