#include "io/inflate.h"

#include <cstdint>
#include <cstring>

#include "util/common.h"

namespace pramana {

namespace {

struct BitReader {
  const uint8_t* data;
  size_t size;
  size_t pos = 0;
  uint32_t bitBuf = 0;
  int bitCount = 0;

  int bits(int need) {
    uint32_t val = bitBuf;
    while (bitCount < need) {
      if (pos >= size) throw PramanaError("gzip: unexpected end of data");
      val |= static_cast<uint32_t>(data[pos++]) << bitCount;
      bitCount += 8;
    }
    bitBuf = val >> need;
    bitCount -= need;
    return static_cast<int>(val & ((1u << need) - 1));
  }
  void alignByte() { bitBuf = 0; bitCount = 0; }
};

// Canonical Huffman decoding table: counts per length + symbols sorted by code.
struct Huffman {
  short count[16] = {};
  std::vector<short> symbol;

  // Returns 0 on complete code, >0 incomplete, <0 over-subscribed.
  int build(const short* lengths, int n) {
    std::memset(count, 0, sizeof count);
    symbol.assign(n, 0);
    for (int s = 0; s < n; ++s) count[lengths[s]]++;
    if (count[0] == n) return 0;
    int left = 1;
    for (int len = 1; len < 16; ++len) {
      left <<= 1;
      left -= count[len];
      if (left < 0) return left;
    }
    short offs[16];
    offs[1] = 0;
    for (int len = 1; len < 15; ++len) offs[len + 1] = static_cast<short>(offs[len] + count[len]);
    for (int s = 0; s < n; ++s)
      if (lengths[s] != 0) symbol[offs[lengths[s]]++] = static_cast<short>(s);
    return left;
  }

  int decode(BitReader& br) const {
    int code = 0, first = 0, index = 0;
    for (int len = 1; len < 16; ++len) {
      code |= br.bits(1);
      int c = count[len];
      if (code - c < first) return symbol[index + (code - first)];
      index += c;
      first += c;
      first <<= 1;
      code <<= 1;
    }
    throw PramanaError("gzip: invalid Huffman code");
  }
};

const short kLenBase[29] = {3,  4,  5,  6,  7,  8,  9,  10, 11,  13,  15,  17,  19,  23, 27,
                            31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
const short kLenExtra[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2,
                             2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
const short kDistBase[30] = {1,   2,   3,   4,   5,   7,    9,    13,   17,   25,
                             33,  49,  65,  97,  129, 193,  257,  385,  513,  769,
                             1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
const short kDistExtra[30] = {0, 0, 0, 0, 1, 1, 2, 2,  3,  3,  4,  4,  5,  5,  6,
                              6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

void inflateCodes(BitReader& br, std::string& out, const Huffman& lencode, const Huffman& distcode) {
  for (;;) {
    int sym = lencode.decode(br);
    if (sym < 256) {
      out.push_back(static_cast<char>(sym));
    } else if (sym == 256) {
      return;
    } else {
      sym -= 257;
      if (sym >= 29) throw PramanaError("gzip: bad length symbol");
      int len = kLenBase[sym] + br.bits(kLenExtra[sym]);
      int dsym = distcode.decode(br);
      if (dsym >= 30) throw PramanaError("gzip: bad distance symbol");
      size_t dist = static_cast<size_t>(kDistBase[dsym] + br.bits(kDistExtra[dsym]));
      if (dist > out.size()) throw PramanaError("gzip: distance too far back");
      size_t from = out.size() - dist;
      for (int k = 0; k < len; ++k) out.push_back(out[from + k]);
    }
  }
}

void inflateStored(BitReader& br, std::string& out) {
  br.alignByte();
  if (br.pos + 4 > br.size) throw PramanaError("gzip: truncated stored block");
  unsigned len = br.data[br.pos] | (br.data[br.pos + 1] << 8);
  unsigned nlen = br.data[br.pos + 2] | (br.data[br.pos + 3] << 8);
  br.pos += 4;
  if (len != (~nlen & 0xffff)) throw PramanaError("gzip: stored length mismatch");
  if (br.pos + len > br.size) throw PramanaError("gzip: truncated stored data");
  out.append(reinterpret_cast<const char*>(br.data + br.pos), len);
  br.pos += len;
}

void inflateFixed(BitReader& br, std::string& out) {
  static Huffman lencode, distcode;
  static bool built = false;
  if (!built) {
    short lengths[288];
    int s = 0;
    for (; s < 144; ++s) lengths[s] = 8;
    for (; s < 256; ++s) lengths[s] = 9;
    for (; s < 280; ++s) lengths[s] = 7;
    for (; s < 288; ++s) lengths[s] = 8;
    lencode.build(lengths, 288);
    for (s = 0; s < 30; ++s) lengths[s] = 5;
    distcode.build(lengths, 30);
    built = true;
  }
  inflateCodes(br, out, lencode, distcode);
}

void inflateDynamic(BitReader& br, std::string& out) {
  static const short order[19] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
  int nlen = br.bits(5) + 257;
  int ndist = br.bits(5) + 1;
  int ncode = br.bits(4) + 4;
  if (nlen > 286 || ndist > 30) throw PramanaError("gzip: bad dynamic header");
  short lengths[320] = {};
  for (int i = 0; i < ncode; ++i) lengths[order[i]] = static_cast<short>(br.bits(3));
  Huffman lencode, distcode;
  if (lencode.build(lengths, 19) != 0) throw PramanaError("gzip: incomplete code-length code");
  int index = 0;
  while (index < nlen + ndist) {
    int sym = lencode.decode(br);
    if (sym < 16) {
      lengths[index++] = static_cast<short>(sym);
    } else {
      short len = 0;
      int rep;
      if (sym == 16) {
        if (index == 0) throw PramanaError("gzip: repeat with no previous length");
        len = lengths[index - 1];
        rep = 3 + br.bits(2);
      } else if (sym == 17) {
        rep = 3 + br.bits(3);
      } else {
        rep = 11 + br.bits(7);
      }
      if (index + rep > nlen + ndist) throw PramanaError("gzip: too many lengths");
      while (rep--) lengths[index++] = len;
    }
  }
  if (lengths[256] == 0) throw PramanaError("gzip: no end-of-block code");
  int err = lencode.build(lengths, nlen);
  if (err < 0 || (err > 0 && nlen - lencode.count[0] != 1)) throw PramanaError("gzip: bad literal code");
  err = distcode.build(lengths + nlen, ndist);
  if (err < 0 || (err > 0 && ndist - distcode.count[0] != 1)) throw PramanaError("gzip: bad distance code");
  inflateCodes(br, out, lencode, distcode);
}

size_t inflateRaw(const uint8_t* data, size_t size, std::string& out) {
  BitReader br{data, size};
  int last;
  do {
    last = br.bits(1);
    int type = br.bits(2);
    if (type == 0) inflateStored(br, out);
    else if (type == 1) inflateFixed(br, out);
    else if (type == 2) inflateDynamic(br, out);
    else throw PramanaError("gzip: invalid block type");
  } while (!last);
  return br.pos;
}

}  // namespace

bool looksGzipped(const std::string& data) {
  return data.size() >= 2 && static_cast<uint8_t>(data[0]) == 0x1f && static_cast<uint8_t>(data[1]) == 0x8b;
}

std::string gunzip(const std::string& compressed) {
  const uint8_t* d = reinterpret_cast<const uint8_t*>(compressed.data());
  size_t n = compressed.size();
  size_t pos = 0;
  std::string out;
  out.reserve(n * 4);
  while (pos + 10 <= n && d[pos] == 0x1f && d[pos + 1] == 0x8b) {
    if (d[pos + 2] != 8) throw PramanaError("gzip: unsupported compression method");
    int flags = d[pos + 3];
    pos += 10;
    if (flags & 4) {  // FEXTRA
      size_t xlen = d[pos] | (d[pos + 1] << 8);
      pos += 2 + xlen;
    }
    if (flags & 8) { while (pos < n && d[pos]) ++pos; ++pos; }   // FNAME
    if (flags & 16) { while (pos < n && d[pos]) ++pos; ++pos; }  // FCOMMENT
    if (flags & 2) pos += 2;                                     // FHCRC
    if (pos > n) throw PramanaError("gzip: truncated header");
    pos += inflateRaw(d + pos, n - pos, out);
    pos += 8;  // CRC32 + ISIZE
  }
  return out;
}

}  // namespace pramana
