class IdentityEncoding:
    """Used internally as an empty temp encoding when encodings stuff."""
    @classmethod
    def is_base_encoding(cls):
        return True

    def is_one_to_one_encoding(self):
        return True

    def encode(self, s):
        return s

    def decode(self, s):
        return s
