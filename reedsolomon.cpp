//  This file is part of par2cmdline (a PAR 2.0 compatible file verification and
//  repair tool). See http://parchive.sourceforge.net for details of PAR 2.0.
//
//  Copyright (c) 2003 Peter Brian Clements
//
//  par2cmdline is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  par2cmdline is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "par2cmdline.h"

#ifdef _MSC_VER
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
#endif

u32 gcd(u32 a, u32 b)
{
  if (a && b)
  {
    while (a && b)
    {
      if (a>b)
      {
        a = a%b;
      }
      else
      {
        b = b%a;
      }
    }

    return a+b;
  }
  else
  {
    return 0;
  }
}

template<> bool ReedSolomon<Galois8>::SetInput(const vector<bool> &present)
{
  inputcount = (u32)present.size();

  datapresentindex = new u32[inputcount];
  datamissingindex = new u32[inputcount];
  database         = new G::ValueType[inputcount];

  G::ValueType base = 1;

  for (unsigned int index=0; index<inputcount; index++)
  {
    // Record the index of the file in the datapresentindex array 
    // or the datamissingindex array
    if (present[index])
    {
      datapresentindex[datapresent++] = index;
    }
    else
    {
      datamissingindex[datamissing++] = index;
    }

    database[index] = base++;
  }

  return true;
}

template<> bool ReedSolomon<Galois8>::SetInput(u32 count)
{
  inputcount = count;

  datapresentindex = new u32[inputcount];
  datamissingindex = new u32[inputcount];
  database         = new G::ValueType[inputcount];

  G::ValueType base = 1;

  for (unsigned int index=0; index<count; index++)
  {
    // Record that the file is present
    datapresentindex[datapresent++] = index;

    database[index] = base++;
  }

  return true;
}

template<> bool ReedSolomon<Galois8>::InternalProcess(const Galois8 factor, size_t size, const void *inputbuffer, void *outputbuffer)
{
#ifdef LONGMULTIPLY  // Based on the Galois16 version, stripped down.  just lookups in a fully pre-computed multiplication table
  // The 8-bit long multiplication tables
  Galois8 *table = glmt->tables;

  // Split the factor into Low and High bytes
  unsigned int fl = (factor.Value() >> 0) & 0xff;

  Galois8 *LL = &table[(0*256 + fl) * 256 + 0]; // factor.low  * source.low

  unsigned int L[256];
  unsigned int *pL = &L[0];
  for (unsigned int i=0; i<256; i++)
    *pL++ = (*LL++).Value();	// expand the bytes to 32bits, for no reason


  // Treat the buffers as arrays of 32-bit unsigned ints.
  u32 *src4 = (u32 *)inputbuffer;
  u32 *end4 = (u32 *)&((u8*)inputbuffer)[size & ~3];
  u32 *dst4 = (u32 *)outputbuffer;

  // Process the data
  while (src4 < end4)
  {
    u32 s = *src4++;

    // Use the lookup table
    *dst4++ ^= (L[(s >> 0) & 0xff]      )
            ^  (L[(s >> 8) & 0xff] << 8 )
            ^  (L[(s >> 16)& 0xff] << 16)
            ^  (L[(s >> 24)& 0xff] << 24);
  }

  // Process any left over bytes at the end of the buffer
  if (size & 3)
  {
    u8 *src1 = &((u8*)inputbuffer)[size & ~3];
    u8 *end1 = &((u8*)inputbuffer)[size];
    u8 *dst1 = &((u8*)outputbuffer)[size & ~3];

    // Process the data
    while (src1 < end1)
    {
      u8 s = *src1++;
      *dst1++ ^= L[s];
    }
  }
#else
  // Treat the buffers as arrays of 8-bit Galois values.

  Galois8 *src = (Galois8 *)inputbuffer;
  Galois8 *end = (Galois8 *)&((u8*)inputbuffer)[size];
  Galois8 *dst = (Galois8 *)outputbuffer;

  // Process the data
  while (src < end)
  {
    *dst++ += *src++ * factor;
  }
#endif

  return eSuccess;
}



////////////////////////////////////////////////////////////////////////////////////////////



// Set which of the source files are present and which are missing
// and compute the base values to use for the vandermonde matrix.
template<> bool ReedSolomon<Galois16>::SetInput(const vector<bool> &present)
{
  inputcount = (u32)present.size();

  datapresentindex = new u32[inputcount];
  datamissingindex = new u32[inputcount];
  database         = new G::ValueType[inputcount];

  unsigned int logbase = 0;

  for (unsigned int index=0; index<inputcount; index++)
  {
    // Record the index of the file in the datapresentindex array 
    // or the datamissingindex array
    if (present[index])
    {
      datapresentindex[datapresent++] = index;
    }
    else
    {
      datamissingindex[datamissing++] = index;
    }

    // Determine the next useable base value.
    // Its log must must be relatively prime to 65535
    while (gcd(G::Limit, logbase) != 1)
    {
      logbase++;
    }
    if (logbase >= G::Limit)
    {
      cerr << "Too many input blocks for Reed Solomon matrix." << endl;
      return false;
    }
    G::ValueType base = G(logbase++).ALog();

    database[index] = base;
  }

  return true;
}

// Record that the specified number of source files are all present
// and compute the base values to use for the vandermonde matrix.
template<> bool ReedSolomon<Galois16>::SetInput(u32 count)
{
  inputcount = count;

  datapresentindex = new u32[inputcount];
  datamissingindex = new u32[inputcount];
  database         = new G::ValueType[inputcount];

  unsigned int logbase = 0;

  for (unsigned int index=0; index<count; index++)
  {
    // Record that the file is present
    datapresentindex[datapresent++] = index;

    // Determine the next useable base value.
    // Its log must must be relatively prime to 65535
    while (gcd(G::Limit, logbase) != 1)
    {
      logbase++;
    }
    if (logbase >= G::Limit)
    {
      cerr << "Too many input blocks for Reed Solomon matrix." << endl;
      return false;
    }
    G::ValueType base = G(logbase++).ALog();

    database[index] = base;
  }

  return true;
}

#if __BYTE_ORDER == __LITTLE_ENDIAN
#define le16_to_cpu(x) (x)
#define cpu_to_le16(x) (x)
#elif __BYTE_ORDER == __BIG_ENDIAN
static inline u16 le16_to_cpu(u16 temp) { return (temp >> 8) & 0xff | (temp << 8) & 0xff00; }
#define cpu_to_le16 le16_to_cpu
#endif
// TODO: get the implementation from linux byteorder.h or whatever it is.
// Or not, telling the compiler about our byteswapping is probably optimal
// using an inline asm bswap or whatever would hide the data movement from the compiler.


template<> bool ReedSolomon<Galois16>::InternalProcess(const Galois16 factor, size_t size, const void *inputbuffer, void *outputbuffer)
{
//#undef LONGMULTIPLY
#ifdef LONGMULTIPLY
  // The 8-bit long multiplication tables
  const Galois16 *table = glmt->tables;

  // Split the factor into Low and High bytes
  unsigned int fl = (factor.Value() >> 0) & 0xff;
  unsigned int fh = (factor.Value() >> 8) & 0xff;

  // Get the four separate multiplication tables
  const Galois16 *LL = &table[(0*256 + fl) * 256 + 0]; // factor.low  * source.low
  const Galois16 *LH = &table[(1*256 + fl) * 256 + 0]; // factor.low  * source.high
  const Galois16 *HL = &table[(1*256 + 0) * 256 + fh]; // factor.high * source.low
  const Galois16 *HH = &table[(2*256 + fh) * 256 + 0]; // factor.high * source.high

  // Combine the four multiplication tables into two
  unsigned int L[256];
  unsigned int H[256];

#if __BYTE_ORDER == __LITTLE_ENDIAN
  unsigned int *pL = &L[0];
  unsigned int *pH = &H[0];
#else
  unsigned int *pL = &H[0];
  unsigned int *pH = &L[0];
#endif

  for (unsigned int i=0; i<256; i++)
  {
    *pL = cpu_to_le16( (*LL + *HL).Value() );
    pL++;    LL++;    HL+=256;

    *pH = cpu_to_le16( (*LH + *HH).Value() );
    pH++;    LH++;    HH++;
  }

  // Treat the buffers as arrays of 32-bit unsigned ints.
  u32 *src = (u32 *)inputbuffer;
  u32 *end = (u32 *)&((u8*)inputbuffer)[size];
  u32 *dst = (u32 *)outputbuffer;
  
  // Process the data
  while (src < end)
  {
    u32 s = *src++;

    // Use the two lookup tables computed earlier
//#if __BYTE_ORDER == __LITTLE_ENDIAN	// in the BE case, things are swapped around so this is still right.
    u32 d = *dst ^ (L[(s >> 0) & 0xff]      )
                 ^ (H[(s >> 8) & 0xff]      )
                 ^ (L[(s >> 16)& 0xff] << 16)
                 ^ (H[(s >> 24)& 0xff] << 16);
    *dst++ = d;
//#else
//    *dst++ ^= (L[(s >> 8) & 0xff]      )
//           ^  (H[(s >> 0) & 0xff]      )
//           ^  (L[(s >> 24)& 0xff] << 16)
//           ^  (H[(s >> 16)& 0xff] << 16);
//#endif
  }
#else
  // Treat the buffers as arrays of 16-bit Galois values.

  const u16 *src = (const u16*)inputbuffer;
  const u16 *end = (const u16*)&((const u8*)inputbuffer)[size];
  u16 *dst = (u16*) outputbuffer;

  // Process the data
  while (src < end)
  {
      Galois16 stmp = Galois16(le16_to_cpu(*src));
      Galois16 dtmp = Galois16(le16_to_cpu(*dst));
      Galois16 newval = dtmp + stmp * factor;
      *dst = cpu_to_le16(newval.Value());
      ++dst; ++src;
  }
#endif

  return eSuccess;
}

