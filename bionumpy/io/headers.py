from typing import Dict


class BaseHeader:
    def __init__(self):
        self._header = ""

    def read_header(self):
        pass


class SamHeader:
    def __init__(self, header_text: str, contig_dict: Dict[str, int]):
        self._header_text = header_text
        self._contig_dict = contig_dict

    @property
    def contig_dict(self):
        return self._contig_dict

    @classmethod
    def from_text(self, text):
        contig_lines = (line for line in text.split('\n') if line.startswith('@SQ'))
        contig_dict = dict(self._get_name_and_length(line) for line in contig_lines)
        return SamHeader(text, contig_dict)

    @classmethod
    def _get_name_and_length(cls, line):
        fields = dict(part.split(':', maxsplit=1) for part in line.split()[1:])
        return fields['SN'], int(fields['LN'])
