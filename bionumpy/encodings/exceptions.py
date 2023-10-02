class EncodingError(Exception):
    def __init__(self, message, offset=0):
        self.message = message
        self.offset = offset
