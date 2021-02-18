import sys
from struct import unpack


def remove_ff00(data: bytes):
    ret = []
    i = 0
    while True:
        b, bnext = unpack("BB", data[i:i + 2])
        if b == 0xff:
            if bnext != 0:
                break
            ret.append(data[i])
            i += 2
        else:
            ret.append(data[i])
            i += 1
    return ret, i


class JPEGBaselineDecoder:
    def __init__(self, filename: str):
        self.f = open(filename, 'rb')
        self.q_table = [[], [], [], []]
        self.xyp = (0, 0, 0)
        self.huffman_ac_tables = [{}, {}, {}, {}]
        self.huffman_dc_tables = [{}, {}, {}, {}]
        self.component = {}
        self.n_components = 0
        self.dc = []

    def read_word(self) -> int:
        """ Read a 16 bit word"""
        out = ord(self.f.read(1)) << 8
        out |= ord(self.f.read(1))
        return out

    def read_byte(self) -> int:
        """Read a byte"""
        out = ord(self.f.read(1))
        return out

    def read_app(self, app_type):
        """Read APP marker"""
        Lp = self.read_word()
        Lp -= 2
        # skip app header
        self.f.seek(Lp, 1)

    def read_dqt(self):
        """Read the quantization table. The table is in zigzag order"""
        Lq = self.read_word()
        Lq -= 2
        while Lq > 0:
            table = []
            Tq = self.read_byte()
            Pq = Tq >> 4
            Tq &= 0xF
            Lq -= 1

            if Pq == 0:
                for i in range(64):
                    table.append(self.read_byte())
                    Lq -= 1
            else:
                for i in range(64):
                    val = self.read_word()
                    table.append(val)
                    Lq -= 2

            self.q_table[Tq] = table

    def read_dnl(self):
        """Read the DNL marker. Changes the number of lines"""
        Ld = self.read_word()
        Ld -= 2
        NL = self.read_word()
        Ld -= 2

        X, Y, P = self.xyp

        if Y == 0:
            self.xyp = X, NL, P

    def huffman_codes(self, huff_lens: list):
        """Calculate the huffman code of each length"""
        huffcode = []
        code = 0

        # magic
        n = len(huff_lens)
        for i in range(n):
            si = huff_lens[i]
            for _ in range(si):
                huffcode.append((i + 1, code))
                code += 1

            code <<= 1

        return huffcode

    def map_codes_to_values(self, codes: list, values: list):
        """Map the huffman code to the right value"""
        out = {}

        n = len(codes)
        for i in range(n):
            out[codes[i]] = values[i]

        return out

    def read_dht(self):
        """Read and compute the huffman tables"""
        # Read the marker length
        Lh = self.read_word()
        Lh -= 2  # minus the size of DHT marker itself, which is 16 bits
        while Lh > 0:
            huff_lens = []
            huff_vals = []
            print(f"DHT Table length: {Lh}")

            T = self.read_byte()
            Th = T & 0x0F

            print(f"DHT Table destination: {Th}")
            Tc = (T >> 4) & 0x0F

            print(f"DHT Table class: {Tc}")
            Lh = Lh - 1

            # read how many symbols of each length up to 16 bits
            for i in range(16):
                huff_lens.append(self.read_byte())
                Lh -= 1

            # generate the huffman codes
            huffcode = self.huffman_codes(huff_lens)
            print("Huffman codes:", huffcode)

            # read the values that should be mapped to huffman codes
            for _ in huffcode:
                huff_vals.append(self.read_byte())
                Lh -= 1

            # generate lookup tables
            if Tc == 0:
                self.huffman_dc_tables[Th] = self.map_codes_to_values(huffcode, huff_vals)
            else:
                self.huffman_ac_tables[Th] = self.map_codes_to_values(huffcode, huff_vals)

    def read_sof(self, type_):
        """Read the start of frame marker"""
        # read the marker length
        Lf = self.read_word()
        Lf -= 2
        # read the sample precision
        P = self.read_byte()
        Lf -= 1
        # read number of lines
        Y = self.read_word()
        Lf -= 2
        # read the number of sample per line
        X = self.read_word()
        Lf -= 2
        # Read number of components
        Nf = self.read_byte()
        Lf -= 1

        self.xyp = X, Y, P
        print(self.xyp)

        while Lf > 0:
            # Read component identifier
            C = self.read_byte()
            # Read sampling factors
            V = self.read_byte()
            Tq = self.read_byte()
            Lf -= 3
            H = V >> 4
            V &= 0xF
            self.component[C] = {}
            # assign horizontal sampling factor
            self.component[C]['H'] = H
            # assign vertical sampling factor
            self.component[C]['V'] = V
            # assign quantization table
            self.component[C]['Tq'] = Tq

    def read_sos(self):
        """Read the start of scan marker"""
        Ls = self.read_word()
        Ls -= 2

        # read number of components in scan
        Ns = self.read_byte()
        print("Number of components", Ns)
        Ls -= 1

        for i in range(Ns):
            # read the scan component selector
            Cs = self.read_byte()
            print("Scan component selector:", Cs)
            Ls -= 1

            # read the huffman table selectors
            Ta = self.read_byte()
            Ls -= 1
            Td = Ta >> 4
            Ta &= 0xF
            print("AC entropy coding table destination selector:", Ta)
            print("DC entropy coding table destination selector:", Td)

            # assign the DC and AC huffman table
            self.component[Cs]['Td'] = Td
            self.component[Cs]['Ta'] = Ta

        # should be zero if baseline DCT
        Ss = self.read_byte()
        Ls -= 1
        # should be 63 if baseline DCT
        Se = self.read_byte()
        Ls -= 1
        # should be zero if baseline DCT
        A = self.read_byte()
        Ls -= 1

        assert Ss == 0
        assert Se == 63
        assert A == 0

        self.n_components = Ns
        self.dc = [0 for _ in range(self.n_components + 1)]

    def decode(self):
        ch = self.read_byte()

        while ch:
            if ch == 0xFF:  # prefix of a marker
                ch = self.read_byte()
                if ch == 0xD8:
                    print("Start of image")
                elif 0xE0 <= ch <= 0xEF:
                    print("APP%x" % (ch - 0xe0))
                    self.read_app(ch - 0XE0)
                elif 0xd0 <= ch <= 0xD7:
                    print("RST%x" % (ch - 0XD0))
                elif ch == 0xDB:
                    print("DQT")
                    self.read_dqt()
                elif ch == 0xDC:
                    print("DNL")
                    self.read_dnl()
                elif ch == 0xC4:
                    print("DHT")
                    self.read_dht()
                elif ch == 0xC8:
                    print("JPG")
                elif ch == 0xCC:
                    print("DAC")
                elif 0XC0 <= ch <= 0XCF:
                    print("Start of frame %x" % (ch - 0xC0))
                    self.read_sof(ch - 0xC0)
                elif ch == 0xDA:
                    print("Start of scan")
                    self.read_sos()
                    # bit_stream = bit_read(input_file)
                    # while not EOI:
                    #     data.append(read_mcu())
                elif ch == 0xD9:
                    print("End of image")
                # print("FF%02X" % in_num)
            ch = self.read_byte()


def main():
    filename = sys.argv[1]
    decoder = JPEGBaselineDecoder(filename)
    decoder.decode()


if __name__ == "__main__":
    main()
