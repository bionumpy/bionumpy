from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Set, Tuple, List, Mapping, Optional

import re
import warnings


@dataclass
class VCFHeader:
    """VCFHeader"""
    fileformat: str = None
    source: str = None
    fileDate: str = None
    reference: str = None
    FILTER: Mapping[str, Any] = None
    FORMAT: Mapping[str, Any] = None
    INFO: Mapping[str, Any] = field(default_factory=dict)
    contig: Mapping[str, Any] = None
    optional: Mapping[str, List[Any]] = None


def preprocess_number(x: str) -> Optional[int]:
    """Preprocess the number in field such as Number. Return None if x
    cannot be formated as integer.

    Parameters
    ----------
    x: str
        the value for the field. e.g. Number={$x}

    Returns
    -------
    Optional[int]
        the integer representation of x. None if x cannot be formated 
        as integer.
    """
    regex = re.match('(\d)', x)
    if regex:
        return int(regex.group(1))
    else:
        return None


def preprocess_type(x: str) -> type:
    """Preprocess the type in field such as Type. Reuturn the class 
    type.

    Parameters
    ----------
    x: str
        the value for the field. e.g. Type={$x}

    Returns
    -------
    type
        the corresponding class type.
    """
    mapping = {'Float': Optional[float], 'Integer': Optional[int], 'Flag': bool, 'String': str}
    return mapping[x]


def extract_identifier_and_content(line: str) -> Tuple[Optional[str]]:
    """Separate the identifier and the content. Return (None, None) if
    not in the correct VCF header format (##{$identifier}={$content}).

    Parameters
    ----------
    line: str
        the header line. e.g. ##{$identifier}={$content}

    Returns
    -------
    Tuple[str]
        the first element is the identifier, the second element is the 
        content.  Return (None, None) if not in the correct VCF header 
        format (##{$identifier}={$content}).

    """
    regex = re.search('^##(\S+?)=(.*)$', line)
    if regex:
        return regex.group(1), regex.group(2)
    else:
        message = "".join((
            f"Header is not in the correct format. Expect ",
            f'^##(\S+?)=(.*)$, get {line}. Ignore this line.'))
        warnings.warn(message, RuntimeWarning)
        return None, None

def parse_mapping_header(
        content: str,
        preprocessors: Dict[str, Callable],
        field_regex: str) -> Dict[str, Any]:
    """Parse the mapping header.

    Parameters
    ----------
    content : str
        the content of the identifier. e.g. 
        '##($identifier)=<$content>$'

    preprocessors : Dict[str, Callable]
        preprocessor for the value, the keys should be the field, the
        values shoul be a function that takes the string content as 
        inputs. e.g. {'Type':  preprocess_type}, where 
        preprocess_type('Float')

    identifier_regex : str
        regular expression to extract the string content. e.g. 
        '=(.+?)[,>]' would extract 'Float' from '..Type=Float,..'

    Returns
    -------
    Dict[str, Any]
        the mapping for the identifier.

    """
    results = dict()
    for field, preprocess_fn in preprocessors.items():
        regex = re.search(f'{field}{field_regex[field]}', content)
        if regex:
            results.setdefault(field, preprocess_fn(regex.group(1)))

    return results

class HeaderReflection:
    """HeaderReflection that deals with the mapping of preprocessor, 
    regular expression, string identifiers, mapping identifiers.

    Attributes
    ----------
    string_identifiers: Set[str]
        the set for identifiers that has a string content.

    mapping_identifiers: Set[str]
        the set for identifiers that has a mapping content.

    preprocessor: Dict[str, Callable]
        the preprocessor for different fields in the mapping headers.
    
    regex: Dict[str, str]
        the regular expression used to extract the values of mapping 
        headers.
    """
    string_identifiers: Set[str] = {'fileformat', 'fileDate', 'source', 'reference'}
    mapping_identifiers: Set[str] = {'FILTER', 'FORMAT', 'INFO', 'contig'}

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

    @classmethod
    def is_string_identifier(cls, field: str) -> bool:
        return field in cls.string_identifiers

    @classmethod
    def is_mapping_identifier(cls, field: str) -> bool:
        return field in cls.mapping_identifiers


def parse_header(lines: str) -> VCFHeader:
    """Parse the VCF header.

    Parameters
    ----------
    lines: str
        string of the VCF file.

    Returns
    -------
    VCFHeader
        The parsed VCFHeader
    """
    headers = dict()
    for line in lines.split('\n'):
        if not line.startswith('##'):
            continue
        identifier, content = extract_identifier_and_content(line)
        if identifier is None:
            continue

        if HeaderReflection.is_string_identifier(identifier):
            # for line '##FIELD=STRING
            headers[identifier] = content
        elif HeaderReflection.is_mapping_identifier(identifier):
            # for line '##FIELD=<KEY1=VALUE1,KEY2=VALUE2...>'
            if identifier not in headers:
                headers[identifier] = OrderedDict()

            mapping = parse_mapping_header(
                    content,
                    HeaderReflection.preprocessor,
                    HeaderReflection.regex)

            if mapping.get('ID', None):
                headers[identifier][mapping.get('ID')] = mapping
            else:
                absent_id = 'Without ID'
                if absent_id not in headers[identifier]:
                    headers[identifier][absent_id] = list()
                headers[identifier][absent_id].append(mapping)

        else:
            # for line ##FIELD=STRING where FIELD not in string identifier or mapping identifier
            optional_identifier = 'optional'
            if optional_identifier not in headers:
                headers[optional_identifier] = dict()
            if identifier not in headers[optional_identifier]:
                headers[optional_identifier][identifier] = list()
            headers[optional_identifier][identifier].append(content)
    
    optional_identifier = 'optional'
    # for identifier, contents in headers[optional_identifier].items():
    #     if len(contents) > 1:
    #         message = f"Duplicated values for {identifier} ({len(contents)} values)"
    #         warnings.warn(message, RuntimeWarning)

    return VCFHeader(**headers)