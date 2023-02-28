import re
from ..encoded_array import as_encoded_array

def match_regexp(encoded_array, regexp: str):
    encoded_array = as_encoded_array(encoded_array)
    return as_encoded_array(re.findall(regexp, encoded_array.to_string()))
