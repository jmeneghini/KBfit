#include "byte_handler.h"
#include <cstdlib>
#include <iostream>
#include <stdexcept>

/*
Copyright (C) 1995-1996 Jean-loup Gailly and Mark Adler

This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
claim that you wrote the original software. If you use this software
in a product, an acknowledgment in the product documentation would be
appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.

Jean-loup Gailly        Mark Adler
gzip@prep.ai.mit.edu    madler@alumni.caltech.edu


The data format used by the zlib library is described by RFCs (Request for
Comments) 1950 to 1952 in the files ftp://ds.internic.net/rfc/rfc1950.txt
(zlib format), rfc1951.txt (deflate format) and rfc1952.txt (gzip format).
*/

/*
  Generate a table for a byte-wise 32-bit CRC calculation on the polynomial:
  x^32+x^26+x^23+x^22+x^16+x^12+x^11+x^10+x^8+x^7+x^5+x^4+x^2+x+1.

  Polynomials over GF(2) are represented in binary, one bit per coefficient,
  with the lowest powers in the most significant bit.  Then adding polynomials
  is just exclusive-or, and multiplying a polynomial by x is a right shift by
  one.  If we call the above polynomial p, and represent a byte as the
  polynomial q, also with the lowest power in the most significant bit (so the
  byte 0xb1 is the polynomial x^7+x^3+x+1), then the CRC is (q*x^32) mod p,
  where a mod b means the remainder after dividing a by b.

  This calculation is done using the shift-register method of multiplying and
  taking the remainder.  The register is initialized to zero, and for each
  incoming bit, x^32 is added mod p to the register if the bit is a one (where
  x^32 mod p is p+x^32 = x^26+...+1), and the register is multiplied mod p by
  x (which is shifting right by one and adding x^32 mod p if the bit shifted
  out is a one).  We start with the highest power (least significant bit) of
  q and repeat for all eight bits of q.

  The table is simply the CRC of all possible eight bit values.  This is all
  the information needed to generate CRC's on data a byte at a time for all
  combinations of CRC register values and incoming bytes.
*/

using namespace std;

// **************************************************************

ByteHandler::ByteHandler()
    : crc_table(
          {0x00000000L, 0x77073096L, 0xee0e612cL, 0x990951baL, 0x076dc419L,
           0x706af48fL, 0xe963a535L, 0x9e6495a3L, 0x0edb8832L, 0x79dcb8a4L,
           0xe0d5e91eL, 0x97d2d988L, 0x09b64c2bL, 0x7eb17cbdL, 0xe7b82d07L,
           0x90bf1d91L, 0x1db71064L, 0x6ab020f2L, 0xf3b97148L, 0x84be41deL,
           0x1adad47dL, 0x6ddde4ebL, 0xf4d4b551L, 0x83d385c7L, 0x136c9856L,
           0x646ba8c0L, 0xfd62f97aL, 0x8a65c9ecL, 0x14015c4fL, 0x63066cd9L,
           0xfa0f3d63L, 0x8d080df5L, 0x3b6e20c8L, 0x4c69105eL, 0xd56041e4L,
           0xa2677172L, 0x3c03e4d1L, 0x4b04d447L, 0xd20d85fdL, 0xa50ab56bL,
           0x35b5a8faL, 0x42b2986cL, 0xdbbbc9d6L, 0xacbcf940L, 0x32d86ce3L,
           0x45df5c75L, 0xdcd60dcfL, 0xabd13d59L, 0x26d930acL, 0x51de003aL,
           0xc8d75180L, 0xbfd06116L, 0x21b4f4b5L, 0x56b3c423L, 0xcfba9599L,
           0xb8bda50fL, 0x2802b89eL, 0x5f058808L, 0xc60cd9b2L, 0xb10be924L,
           0x2f6f7c87L, 0x58684c11L, 0xc1611dabL, 0xb6662d3dL, 0x76dc4190L,
           0x01db7106L, 0x98d220bcL, 0xefd5102aL, 0x71b18589L, 0x06b6b51fL,
           0x9fbfe4a5L, 0xe8b8d433L, 0x7807c9a2L, 0x0f00f934L, 0x9609a88eL,
           0xe10e9818L, 0x7f6a0dbbL, 0x086d3d2dL, 0x91646c97L, 0xe6635c01L,
           0x6b6b51f4L, 0x1c6c6162L, 0x856530d8L, 0xf262004eL, 0x6c0695edL,
           0x1b01a57bL, 0x8208f4c1L, 0xf50fc457L, 0x65b0d9c6L, 0x12b7e950L,
           0x8bbeb8eaL, 0xfcb9887cL, 0x62dd1ddfL, 0x15da2d49L, 0x8cd37cf3L,
           0xfbd44c65L, 0x4db26158L, 0x3ab551ceL, 0xa3bc0074L, 0xd4bb30e2L,
           0x4adfa541L, 0x3dd895d7L, 0xa4d1c46dL, 0xd3d6f4fbL, 0x4369e96aL,
           0x346ed9fcL, 0xad678846L, 0xda60b8d0L, 0x44042d73L, 0x33031de5L,
           0xaa0a4c5fL, 0xdd0d7cc9L, 0x5005713cL, 0x270241aaL, 0xbe0b1010L,
           0xc90c2086L, 0x5768b525L, 0x206f85b3L, 0xb966d409L, 0xce61e49fL,
           0x5edef90eL, 0x29d9c998L, 0xb0d09822L, 0xc7d7a8b4L, 0x59b33d17L,
           0x2eb40d81L, 0xb7bd5c3bL, 0xc0ba6cadL, 0xedb88320L, 0x9abfb3b6L,
           0x03b6e20cL, 0x74b1d29aL, 0xead54739L, 0x9dd277afL, 0x04db2615L,
           0x73dc1683L, 0xe3630b12L, 0x94643b84L, 0x0d6d6a3eL, 0x7a6a5aa8L,
           0xe40ecf0bL, 0x9309ff9dL, 0x0a00ae27L, 0x7d079eb1L, 0xf00f9344L,
           0x8708a3d2L, 0x1e01f268L, 0x6906c2feL, 0xf762575dL, 0x806567cbL,
           0x196c3671L, 0x6e6b06e7L, 0xfed41b76L, 0x89d32be0L, 0x10da7a5aL,
           0x67dd4accL, 0xf9b9df6fL, 0x8ebeeff9L, 0x17b7be43L, 0x60b08ed5L,
           0xd6d6a3e8L, 0xa1d1937eL, 0x38d8c2c4L, 0x4fdff252L, 0xd1bb67f1L,
           0xa6bc5767L, 0x3fb506ddL, 0x48b2364bL, 0xd80d2bdaL, 0xaf0a1b4cL,
           0x36034af6L, 0x41047a60L, 0xdf60efc3L, 0xa867df55L, 0x316e8eefL,
           0x4669be79L, 0xcb61b38cL, 0xbc66831aL, 0x256fd2a0L, 0x5268e236L,
           0xcc0c7795L, 0xbb0b4703L, 0x220216b9L, 0x5505262fL, 0xc5ba3bbeL,
           0xb2bd0b28L, 0x2bb45a92L, 0x5cb36a04L, 0xc2d7ffa7L, 0xb5d0cf31L,
           0x2cd99e8bL, 0x5bdeae1dL, 0x9b64c2b0L, 0xec63f226L, 0x756aa39cL,
           0x026d930aL, 0x9c0906a9L, 0xeb0e363fL, 0x72076785L, 0x05005713L,
           0x95bf4a82L, 0xe2b87a14L, 0x7bb12baeL, 0x0cb61b38L, 0x92d28e9bL,
           0xe5d5be0dL, 0x7cdcefb7L, 0x0bdbdf21L, 0x86d3d2d4L, 0xf1d4e242L,
           0x68ddb3f8l, 0x1fda836eL, 0x81be16cdL, 0xf6b9265bL, 0x6fb077e1L,
           0x18b74777L, 0x88085ae6L, 0xff0f6a70L, 0x66063bcaL, 0x11010b5cL,
           0x8f659effL, 0xf862ae69L, 0x616bffd3L, 0x166ccf45L, 0xa00ae278L,
           0xd70dd2eeL, 0x4e048354L, 0x3903b3c2L, 0xa7672661L, 0xd06016f7L,
           0x4969474dL, 0x3e6e77dbL, 0xaed16a4aL, 0xd9d65adcL, 0x40df0b66L,
           0x37d83bf0L, 0xa9bcae53L, 0xdebb9ec5L, 0x47b2cf7fL, 0x30b5ffe9L,
           0xbdbdf21cL, 0xcabac28aL, 0x53b39330L, 0x24b4a3a6L, 0xbad03605L,
           0xcdd70693L, 0x54de5729L, 0x23d967bfL, 0xb3667a2eL, 0xc4614ab8L,
           0x5d681b02L, 0x2a6f2b94L, 0xb40bbe37L, 0xc30c8ea1L, 0x5a05df1bL,
           0x2d02ef8dL}) {}

// Is the native byte order big endian?

bool ByteHandler::big_endian() {
  union {
    int l;
    char c[sizeof(int)];
  } u;
  u.l = 1;
  return (u.c[sizeof(int) - 1] == 1);
}

// Byte-swap an array of data each of size nmemb

void ByteHandler::byte_swap(void* ptr, size_t size, size_t nmemb) {
  unsigned int j;

  char char_in[16]; /* characters used in byte swapping */

  char* in_ptr;
  double* double_ptr; /* Pointer used in the double routines */

  switch (size) {
  case 4: /* n_uint32_t */
  {
    n_uint32_t* w = (n_uint32_t*)ptr;
    n_uint32_t old, recent;

    for (j = 0; j < nmemb; j++) {
      old = w[j];
      recent = old >> 24 & 0x000000ff;
      recent |= old >> 8 & 0x0000ff00;
      recent |= old << 8 & 0x00ff0000;
      recent |= old << 24 & 0xff000000;
      w[j] = recent;
    }
  } break;

  case 1: /* n_uint8_t: byte - do nothing */
    break;

  case 8: /* n_uint64_t */
  {
    for (j = 0, double_ptr = (double*)ptr; j < nmemb; j++, double_ptr++) {
      in_ptr = (char*)double_ptr; /* Set the character pointer to
                                     point to the start of the double */

      /*
       *  Assign all the byte variables to a character
       */
      char_in[0] = in_ptr[0];
      char_in[1] = in_ptr[1];
      char_in[2] = in_ptr[2];
      char_in[3] = in_ptr[3];
      char_in[4] = in_ptr[4];
      char_in[5] = in_ptr[5];
      char_in[6] = in_ptr[6];
      char_in[7] = in_ptr[7];

      /*
       *  Now just swap the order
       */
      in_ptr[0] = char_in[7];
      in_ptr[1] = char_in[6];
      in_ptr[2] = char_in[5];
      in_ptr[3] = char_in[4];
      in_ptr[4] = char_in[3];
      in_ptr[5] = char_in[2];
      in_ptr[6] = char_in[1];
      in_ptr[7] = char_in[0];
    }
  } break;

  case 16: /* Long Long */
  {
    for (j = 0, double_ptr = (double*)ptr; j < nmemb; j++, double_ptr += 2) {

      in_ptr = (char*)double_ptr; /* Set the character pointer to
                                     point to the start of the double */

      /*
       *  Assign all the byte variables to a character
       */
      char_in[0] = in_ptr[0];
      char_in[1] = in_ptr[1];
      char_in[2] = in_ptr[2];
      char_in[3] = in_ptr[3];
      char_in[4] = in_ptr[4];
      char_in[5] = in_ptr[5];
      char_in[6] = in_ptr[6];
      char_in[7] = in_ptr[7];
      char_in[8] = in_ptr[8];
      char_in[9] = in_ptr[9];
      char_in[10] = in_ptr[10];
      char_in[11] = in_ptr[11];
      char_in[12] = in_ptr[12];
      char_in[13] = in_ptr[13];
      char_in[14] = in_ptr[14];
      char_in[15] = in_ptr[15];

      /*
       *  Now just swap the order
       */
      in_ptr[0] = char_in[15];
      in_ptr[1] = char_in[14];
      in_ptr[2] = char_in[13];
      in_ptr[3] = char_in[12];
      in_ptr[4] = char_in[11];
      in_ptr[5] = char_in[10];
      in_ptr[6] = char_in[9];
      in_ptr[7] = char_in[8];

      in_ptr[8] = char_in[7];
      in_ptr[9] = char_in[6];
      in_ptr[10] = char_in[5];
      in_ptr[11] = char_in[4];
      in_ptr[12] = char_in[3];
      in_ptr[13] = char_in[2];
      in_ptr[14] = char_in[1];
      in_ptr[15] = char_in[0];
    }
  } break;

  case 2: /* n_uint16_t */
  {
    n_uint16_t* w = (n_uint16_t*)ptr;
    n_uint16_t old, recent;

    for (j = 0; j < nmemb; j++) {
      old = w[j];
      recent = old >> 8 & 0x00ff;
      recent |= old << 8 & 0xff00;
      w[j] = recent;
    }
  } break;

  default:
    throw(std::invalid_argument("unsupported word size"));
  }
}

ByteHandler::n_uint32_t ByteHandler::get_checksum(n_uint32_t crc,
                                                  const unsigned char* buf,
                                                  size_t len) {
  crc = crc ^ 0xffffffffL;
  while (len >= 8) {
    for (int k = 0; k < 8; ++k)
      crc = crc_table[((int)crc ^ (*buf++)) & 0xff] ^ (crc >> 8);
    len -= 8;
  }
  if (len) {
    do {
      crc = crc_table[((int)crc ^ (*buf++)) & 0xff] ^ (crc >> 8);
    } while (--len);
  }
  return crc ^ 0xffffffffL;
}

ByteHandler::n_uint32_t ByteHandler::get_checksum(n_uint32_t crc,
                                                  const char* buf, size_t len) {
  return get_checksum(crc, (const unsigned char*)(buf), len);
}

// ****************************************************************************************
