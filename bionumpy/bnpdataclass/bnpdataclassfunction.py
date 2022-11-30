class bnpdataclassfunction:
    def __init__(self, *args):
        arg_names = args

    def __call__(self, func):
        def new_func(data_object, *args, **kwargs):
            pass
