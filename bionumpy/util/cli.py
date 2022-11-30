import sys
def run_as_commandline(function):
    args = sys.argv[1:]
    args = (_type(arg) for _type, arg in zip(function.__annotations__.values(), args))
    function(*args)
