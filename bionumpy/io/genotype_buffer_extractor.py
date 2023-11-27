from bionumpy.io.file_buffers import TextBufferExtractor


class GenotypeDataExtractor(TextBufferExtractor):
    def __init__(self, data, field_starts, field_lens):
        super().__init__(data, field_starts, field_lens=field_lens)
