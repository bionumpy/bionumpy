import bionumpy as bnp
from bionumpy.strops import str_equal, split
from bionumpy.gtf import get_transcripts


data = bnp.open("example_data/small.gtf").read()

transcripts = get_transcripts(data)
print(transcripts)
# 
#  data[str_equal(data.feature_type, "transcript")]
# print(transcripts)
# all_features = split(transcripts.atributes.ravel(), " ")
# print(all_features)
# transcript_id = all_features[1::2][str_equal(all_features[:-1:2], "transcript_id")]
# transcript_id = all_features[1::2][str_equal(all_features[:-1:2], "gene_id")]
# 
# 
# print(transcript_id)
