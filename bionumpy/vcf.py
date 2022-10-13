
def parse_info_line(line):
    line = line.strip()
    assert line.startswith("##INFO=<")
    assert line.endswith(">")
    line = line[8:-1]
    parts = line.split(",")[:3]
    return dict(tuple(part.split("=")) for part in parts)

def parse_header(lines):
    info= [parse_info_line(line) for line in lines if line.startswith("##INFO")]
    return info

def read_header(file_obj):
    lines = []
    
    
