class BandedSymmetricMatrix:
    def __init__(self, size, half_bandwidth):
        self.n = size
        self.bw = half_bandwidth
        self.data = [[0.0 for _ in range(self.bw + 1)] for _ in range(self.n)]

    def assemble(self, i, j, value):
        if i > j: i, j = j, i
        offset = j - i
        if 0 <= offset <= self.bw:
            self.data[i][offset] += value

    def set_element(self, i, j, value):
        """Reset and set a specific element to a new value."""
        if i > j: i, j = j, i
        offset = j - i
        if 0 <= offset <= self.bw:
            self.data[i][offset] = value

    def get_element(self, i, j):
        if i > j: i, j = j, i
        offset = j - i
        if 0 <= offset <= self.bw:
            return self.data[i][offset]
        return 0.0
