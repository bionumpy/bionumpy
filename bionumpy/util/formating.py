
def table(nested_list, headers=None, column_width=25):
    lines = []
    col_length = column_width
    if headers is not None:
        header_names = []
        for name in headers:
            header_names.append(f"{name:>{col_length}}")
        lines.append("".join(header_names))

    for _, entry in zip(range(10), nested_list):
        cols = [e for e in entry]
        lines.append("".join(f"{str(col)[:col_length - 2]:>{col_length}}" for col in cols))

    return "\n".join(lines)
