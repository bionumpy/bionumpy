from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Callable, Dict, Set, Tuple, List, Mapping

import re
@dataclass
class VCFHeader:
    fileformat: str
    source: str
    fileDate: str
    reference: str
    FILTER: Mapping[str, Any]
    FORMAT: Mapping[str, Any]
    INFO: Mapping[str, Any]
    contig: Mapping[str, Any]
    optional: List[Any]

def preprocess_number(x):
    regex = re.match('(\d)', x)
    if regex:
        return int(regex.group(1))
    else:
        return None

def preprocess_type(x):
    mapping = {'Float': float, 'Integer': int, 'Flag': bool, 'String': str}
    return mapping[x]


def extract_identifier(line: str) -> Tuple[str]:
    regex = re.search('^##(\S+?)=(.*)$', line)
    if regex:
        return regex.group(1), regex.group(2)
    else:
        return 'optional', line[2:]

def parse_mapping_config(
        content: str,
        preprocessors: Dict[str, Callable],
        identifier_regex: str) -> Dict[str, Any]:
        
    results = dict()
    for field, preprocess_fn in preprocessors.items():
        regex = re.search(f'{field}{identifier_regex[field]}', content)
        if regex:
            results.setdefault(field, preprocess_fn(regex.group(1)))
        else:
            results.setdefault(field, None)

    return results

class IdentifiersReflection:
    string_identifiers: Set[str] = {'fileformat', 'fileDate', 'source', 'reference'}
    mapping_identifiers: Set[str] = {'FILTER', 'FORMAT', 'INFO', 'contig'}
    #keys: List[str] = ['ID', 'Number', 'Type', 'Description']

    preprocessor: Dict[str, Callable] = {
        'ID': lambda x: x,
        'Number': preprocess_number,
        'Type': preprocess_type,
        'Description': lambda x: x
    }
    regex: Dict[str, str] = {
        'ID': '=(.+?)[,>]',
        'Number': '=(.+?)[,>]',
        'Type': '=(.+?)[,>]',
        'Description': '="(.+?)"'
    }

def parse_header(lines):
    headers = dict()
    for line in lines.split('\n'):
        # check if line is optional
        identifier, content = extract_identifier(line)
        if identifier in IdentifiersReflection.string_identifiers:
            headers[identifier] = content
        elif identifier in IdentifiersReflection.mapping_identifiers:
            if identifier not in headers:
                headers[identifier] = OrderedDict()
            result = parse_mapping_config(
                    content,
                    IdentifiersReflection.preprocessor,
                    IdentifiersReflection.regex)
            if result.get('ID', None):
                headers[identifier][result.get('ID')] = result
            else:
                if 'Without ID' not in headers[identifier]:
                    headers['Without ID'] = list()
                headers[identifier]['Without ID'].append(result)
        else:
            if 'optional' not in headers:
                headers['optional'] = list()
            headers['optional'].append(line[2:])
    
    return VCFHeader(**headers)


#def parse_info_line(line):
#    line = line.strip()
#    assert line.startswith("##INFO=<")
#    assert line.endswith(">")
#    line = line[8:-1]
#    parts = line.split(",")[:3]
#    return dict(tuple(part.split("=")) for part in parts)

#def read_header(file_obj):
#    lines = []    
