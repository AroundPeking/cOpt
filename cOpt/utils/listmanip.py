def flatten(x, return_index=False):
    '''
    Flattens a nested list and optionally attaches the original indices.

    '''
    def _flatten(x):
        for i, elem in enumerate(x):
            if isinstance(elem, list):
                for sub_elem in _flatten(elem):
                    yield (sub_elem[0], (i,) + sub_elem[1]) if return_index else sub_elem
            else:
                yield (elem, (i,)) if return_index else elem

    assert isinstance(x, list)
    return list(_flatten(x))