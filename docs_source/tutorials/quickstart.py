import bionumpy as bnp
from bionumpy.sequence.string_matcher import StringMatcher, RegexMatcher

#For small exploration, we can define our list of sequences directly:
sequences = bnp.as_sequence_array(["ACGT", "AATGAT"], encoding=bnp.encodings.ACTGEncoding)

#Plain string matching is done with StringMatcher, using rolling_window to get match status for each possible start position:
matcher = StringMatcher("AT", encoding=bnp.encodings.ACTGEncoding)
matches = matcher.rolling_window(sequences)
print("Matches for AT: ", matches)

#One can also efficiently match regular expressions with fixed length (alternative letters and wildcards):
matcher = RegexMatcher("[AG].[AT]", encoding=bnp.encodings.ACTGEncoding)
matches = matcher.rolling_window(sequences)
print("Matches for [AG].[AT]: ", matches)

#For real usages, one would typically read sequences from a file, while matching works the same way:
#...