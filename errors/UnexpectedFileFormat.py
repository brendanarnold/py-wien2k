class UnexpectedFileFormat(Exception):
    def __init__(self, text=None):
        self.text = text
    def __str__(self):
        if self.text is None:
            return repr('Unknown WIEN2k file format error')
        else:
            return repr(self.text)
