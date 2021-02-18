import sys
from struct import unpack
from math import *
from memoize import memoize

mcus_read = 0
inline_dc = 0
idct_precision = 8
EOI = False
data = []

@memoize
def calc_add_bits(len, val):
    """ Calculate the value from the "additional" bits in the huffman data. """
    if (val & (1 << len - 1)):
        pass
    else:
        val -= (1 << len) - 1

    return val


def bit_read(file):
    """Read one bit from file and handle markers and byte stuffing"""
    global EOI
    global dc
    global inline_dc

    in_ = file.read(1)
    while in_ and not EOI:
        if in_ == chr(0xFF):
            cmd = file.read(1)
            if cmd:
                # byte stuffing
                if cmd == chr(0x00):
                    in_ = chr(0xFF)
                # end of image marker
                elif cmd == chr(0xD9):
                    EOI = True
                # restart markers
                elif 0xD0 <= ord(cmd) <= 0xD7 and inline_dc:
                    # reset dc value
                    dc = [0 for i in range(num_components + 1)]
                    in_ = file.read(1)
                else:
                    in_ = file.read(1)
                    print("CMD: %x" % ord(cmd))

        if not EOI:
            for i in range(7, -1, -1):
                # output next bit
                yield (ord(in_) >> i) & 0x01

            in_ = file.read(1)

    while True:
        yield []


def get_bits(num, gen):
    """Get "num" bits from gen"""
    out = 0
    for _ in range(num):
        out <<= 1
        val = gen.__next__()
        if val != []:
            out += val & 0x01
        else:
            return []

    return out


def print_and_pause(fn):
    def new(*args):
        x = fn(*args)
        print(x)
        input()
        return x

    return new


def read_data_unit(comp_num):
    """Read one DU with component id comp_num"""
    global bit_stream
    global component
    global dc

    data = []

    comp = component[comp_num]
    huff_tbl = huffman_dc_tables[comp['Td']]

    # fill data with 64 coefficients
    while len(data) < 64:
        key = 0

        for bits in range(1, 17):
            key_len = []
            key <<= 1
            # get one bit from bit_stream
            val = get_bits(1, bit_stream)
            if val == []:
                break
            key |= val
            # if huffman code exists
            key_len = huff_tbl.get[(bits, key), None]
            if key_len is not None:
                break

        # after getting the DC value switch to the AC table
        huff_tbl = huffman_ac_tables[comp['Ta']]

        if key_len == []:
            print((bits, key, bin(key)), "key not found")
            break
        # If ZRL fill with 16 zero coefficients
        elif key_len == 0xF0:
            for i in range(16):
                data.append(0)
            continue

        # if not DC coefficient
        if len(data) != 0:
            # If End of block
            if key_len == 0x00:
                # Fill the rest of the DU with zeros
                while len(data) < 64:
                    data.append(0)
                break

            # the first part of the AC key_len is the number of leading zeros
            for i in range(key_len >> 4):
                if len(data) < 64:
                    data.append(0)
            key_len &= 0x0F

        if len(data) >= 64:
            break

        if key_len != 0:
            # The rest of key_len is the amount of "additional" bits
            val = get_bits(key_len, bit_stream)
            if val == []:
                break
            # Decode the additional bits
            num = calc_add_bits(key_len, val)

            # experimental, doesn't work right
            if len(data) == 0 and inline_dc:
                # The DC coefficient value is added to the DC value from the corresponding DU in the previous MCU
                num += dc[comp_num]
                dc[comp_num] = num

            data.append(num)
        else:
            data.append(0)

    if len(data) != 64:
        print("Wrong size", len(data))

    return data


def read_mcu():
    """Read an MCU"""
    global component
    global num_components
    global mcus_read

    mcu = [[] for _ in range(num_components)]
    for i in range(num_components):
        comp = component[i + 1]
        # for each DU
        for j in range(comp['H'] * comp['V']):
            if not EOI:
                mcu[i].append(read_data_unit(i + 1))

    mcus_read += 1

    return mcu


def restore_dc(data):
    """
    Restore the DC values as the DC values are coded as the difference from the previous DC value of the same
    component
    """
    dc_prev = [0 for x in range(len(data[0]))]
    out = []

    # For each MCU
    for mcu in data:
        # For each component
        for comp_num in range(len(mcu)):
            # For each DU
            for du in range(len(mcu[comp_num])):
                if mcu[comp_num][du]:
                    mcu[comp_num][du][0] += dc_prev[comp_num]
                    dc_prev[comp_num] = mcu[comp_num][du][0]

        out.append(mcu)

    return out


def dequantify(mcu):
    """ Dequantify MCU """
    global component

    out = mcu

    # For each component
    for c in range(len(out)):
        # For each DU
        for du in range(len(out[c])):
            # For each coefficient
            for i in range(len(out[c][du])):
                # Multiply by the the corresponding
                # value in the quantization table
                out[c][du][i] *= q_table[component[c + 1]['Tq']][i]

    return out


def zagzig(du):
    """ Put the coefficients in the right order """
    map = [[0, 1, 5, 6, 14, 15, 27, 28],
           [2, 4, 7, 13, 16, 26, 29, 42],
           [3, 8, 12, 17, 25, 30, 41, 43],
           [9, 11, 18, 24, 31, 40, 44, 53],
           [10, 19, 23, 32, 39, 45, 52, 54],
           [20, 22, 33, 38, 46, 51, 55, 60],
           [21, 34, 37, 47, 50, 56, 59, 61],
           [35, 36, 48, 49, 57, 58, 62, 63]]

    # Iterate over 8x8
    for x in range(8):
        for y in range(8):
            if map[x][y] < len(du):
                map[x][y] = du[map[x][y]]
            else:
                # If DU is too short
                # This shouldn't happen.
                map[x][y] = 0

    return map


def for_each_du_in_mcu(mcu, func):
    """ Helper function that iterates over all DU's in an MCU and runs "func" on it """
    return [map(func, comp) for comp in mcu]


# @memoize
def C(x):
    if x == 0:
        return 1.0 / sqrt(2.0)
    else:
        return 1.0


# Lookup table to speed up IDCT somewhat
idct_table = [[(C(u) * cos(((2.0 * x + 1.0) * u * pi) / 16.0)) for x in range(idct_precision)] for u in
              range(idct_precision)]
range8 = range(8)
rangeIDCT = range(idct_precision)


def idct(matrix):
    """ Converts from frequency domain ordinary(?) """
    out = [range(8) for i in range(8)]

    # Iterate over every pixel in the block
    for x in range8:
        for y in range8:
            sum = 0

            # Iterate over every coefficient
            # in the DU
            for u in rangeIDCT:
                for v in rangeIDCT:
                    sum += matrix[v][u] * idct_table[u][x] * idct_table[v][y]

            out[y][x] = sum // 4

    return out


def expand(mcu, H, V):
    """ Reverse subsampling """
    Hout = max(H)
    Vout = max(V)
    out = [[[] for x in range(8 * Hout)] for y in range(8 * Vout)]

    for i in range(len(mcu)):
        Hs = Hout // H[i]
        Vs = Vout // V[i]
        Hin = H[i]
        Vin = V[i]
        comp = mcu[i]

        if len(comp) != (Hin * Vin):
            return []

        for v in range(Vout):
            for h in range(Hout):
                for y in range(8):
                    for x in range(8):
                        out[y + v * 8][x + h * 8].append(comp[(h // Hs) + Hin * (v // Vs)][y // Vs][x // Hs])

    return out


def combine_mcu(mcu):
    global num_components

    H = []
    V = []

    for i in range(num_components):
        H.append(component[i + 1]['H'])
        V.append(component[i + 1]['V'])

    return expand(mcu, H, V)


def combine_blocks(data):
    global XYP

    X, Y, P = XYP

    out = [[(0, 0, 0) for x in range(X + 32)] for y in range(Y + 64)]
    offsetx = 0
    offsety = 0

    for block in data:
        ybmax = len(block)
        for yb in range(ybmax):
            xbmax = len(block[yb])
            for xb in range(xbmax):
                out[yb + offsety][xb + offsetx] = tuple(block[yb][xb])
        offsetx += xbmax
        if offsetx > X:
            offsetx = 0
            offsety += ybmax

    return out


def crop_image(data):
    global XYP
    global Xdensity
    global Ydensity

    X, Y, P = XYP
    return [[data[y][x] for x in range(X)] for y in range(Y)]


def clip(x):
    if x > 255:
        return 255
    elif x < 0:
        return 0
    else:
        return int(x)


def clamp(x):
    x = (abs(x) + x) // 2
    if x > 255:
        return 255
    else:
        return x


@memoize
def YCbCr2RGB(Y, Cb, Cr):
    Cred = 0.299
    Cgreen = 0.587
    Cblue = 0.114

    R = Cr * (2 - 2 * Cred) + Y
    B = Cb * (2 - 2 * Cblue) + Y
    G = (Y - Cblue * B - Cred * R) / Cgreen

    return clamp(R + 128), clamp(G + 128), clamp(B + 128)


def YCbCr2Y(Y, Cb, Cr):
    return Y, Y, Y


def for_each_pixel(data, func):
    out = [[0 for pixel in range(len(data[0]))] for line in range(len(data))]

    for line in range(len(data)):
        for pixel in range(len(data[0])):
            out[line][pixel] = func(*data[line][pixel])

    return out


def tuplify(data):
    out = []

    for line in data:
        out.append(tuple(line))

    return tuple(out)


@memoize
def prepare(x, y, z):
    return "#%02x%02x%02x" % (x, y, z)


if __name__ == "__main__":
    print("AC Huffman tables:", huffman_ac_tables)
    print("DC Huffman tables:", huffman_dc_tables)
    print("Quantization tables:", q_table)
    # print "Component table:", component

    if not inline_dc:
        print("restore dc")
        data = restore_dc(data)

    print("dequantify")
    data = map(dequantify, data)

    print("deserialize")
    data = [for_each_du_in_mcu(mcu, zagzig) for mcu in data]

    print("inverse discrete cosine transform")
    data = [for_each_du_in_mcu(mcu, idct) for mcu in data]

    print("combine mcu")
    data = map(combine_mcu, data)

    print("combine blocks")
    data = combine_blocks(data)

    print("crop image")
    data = crop_image(data)

    print("color conversion")
    data = for_each_pixel(data, YCbCr2RGB)
    # data= for_each_pixel(data, YCbCr2Y)

    print("prepare")
    data = for_each_pixel(data, prepare)

    print("tuplify")
    data = tuplify(data)

    # display_image(data)
