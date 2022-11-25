class IdentityEncoding:
    """Used internally as an empty temp encoding when encodings stuff."""
    @classmethod
    def is_base_encoding(cls):
        return True

    def is_one_to_one_encoding(self):
        return True

    @classmethod
    def encode(cls, s):
        return s

    @classmethod
    def decode(cls, s):
        return s
