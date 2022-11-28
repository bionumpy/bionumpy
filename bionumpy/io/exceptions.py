class FormatException(Exception):
    def __init__(self, message, byte_position=None, line_number=None, offending_text=None):
        super().__init__(message)
        self.byte_position = byte_position
        self.line_number = line_number
        self.offending_text = offending_text
