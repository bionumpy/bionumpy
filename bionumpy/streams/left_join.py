from typing import Iterable, Dict, Tuple, Any
def left_join(grouped_left: Iterable[Tuple[str, Any]], grouped_right: Iterable[Tuple[str, Any]]) -> Iterable[Tuple[str, Any, Any]]:
    default_value = None
    name_right, data_right = next(grouped_right, (None, None))
    for name_left, data_left in grouped_left:
        if name_left != name_right:
            yield (name_left, data_left, default_value)
            continue
        yield (name_left, data_left, data_right)
        name_right, data_right = next(grouped_right, (None, None))
    assert name_right is None and data_right is None
    try:
        next_group = next(grouped_right)
        print(next_group)
    except StopIteration:
        pass
    else:
        raise Exception(f'Data left  in right group: {next_group}')

            
        
              
