import os


def get_file_name(filename: str):
    current_directory = os.getcwd()
    last_folder = os.path.basename(current_directory)
    if last_folder == 'tests' or last_folder == 'scripts':
        return os.path.join('..', filename)
    elif last_folder == 'bionumpy':
        return filename

    assert False, "Unknown folder structure, {}".format(current_directory)