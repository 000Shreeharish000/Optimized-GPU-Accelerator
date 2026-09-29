// Self-contained DEFLATE (RFC 1951) / gzip (RFC 1952) decoder so the solver can
// read MIPLIB-style .mps.gz files without linking zlib.
#pragma once

#include <string>
#include <vector>

namespace pramana {

// Decompresses a gzip stream (possibly multi-member). Throws PramanaError.
std::string gunzip(const std::string& compressed);
bool looksGzipped(const std::string& data);

}  // namespace pramana
