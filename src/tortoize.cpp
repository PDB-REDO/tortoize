/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <fstream>
#include <numeric>
#include <vector>

#include "tortoize.hpp"
#include "../dssp/src/DSSP.hpp"

#include "revision.hpp"

namespace fs = std::filesystem;
namespace ba = boost::algorithm;

using json = zeep::json::element;

// --------------------------------------------------------------------
// simple integer compression, based somewhat on MRS code

class OBitStream
{
  public:
	OBitStream(std::vector<uint8_t>& buffer)
		: m_buffer(buffer)
	{
		m_buffer.push_back(0);
	}

	OBitStream(const OBitStream&) = delete;
	OBitStream& operator=(const OBitStream&) = delete;

	void writebit(bool bit)
	{
		if (bit)
			m_buffer.back() |= 1 << m_bitOffset;
		
		if (--m_bitOffset < 0)
		{
			m_buffer.push_back(0);
			m_bitOffset = 7;
		}
	}

	// write fixed size
	void write(uint32_t value, int bits)
	{
		while (bits-- > 0)
		{
			if (value & (1UL << bits))
				m_buffer.back() |= 1 << m_bitOffset;
			
			if (--m_bitOffset < 0)
			{
				m_buffer.push_back(0);
				m_bitOffset = 7;
			}
		}
	}

	void sync()
	{
		writebit(0);

		while (m_bitOffset != 7)
			writebit(1);
	}

	const uint8_t* data() const					{ return m_buffer.data(); }
	size_t size() const							{ return m_buffer.size(); }

	friend void WriteArray(OBitStream& bs, const std::vector<uint32_t>& data);

  private:
	std::vector<uint8_t>& m_buffer;
	int m_bitOffset = 7;
};

class IBitStream
{
  public:
	IBitStream(const uint8_t* data)
		: m_data(data), m_byte(*m_data++), m_bitOffset(7) {}

	IBitStream(const OBitStream& bits)
		: IBitStream(bits.data()) {}

	IBitStream(const IBitStream&) = delete;
	IBitStream& operator=(const IBitStream&) = delete;

	uint32_t read(int bc)
	{
		uint32_t result = 0;

		while (bc > 0)
		{
			static const uint8_t kM[] = { 0x00, 0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3F, 0x7F, 0xFF };

			int bw = m_bitOffset + 1;
			if (bw > bc)
				bw = bc;

			m_bitOffset -= bw;
			result = result << bw | (kM[bw] & (m_byte >> (m_bitOffset + 1)));

			if (m_bitOffset < 0)
			{
				m_byte = *m_data++;
				m_bitOffset = 7;
			}

			bc -= bw;
		}

		return result;
	}

	friend std::vector<uint32_t> ReadArray(IBitStream& bs);

  private:
	const uint8_t* m_data;
	uint8_t m_byte;
	int m_bitOffset;
};

// --------------------------------------------------------------------
//    Arrays
//    This is a simplified version of the array compression routines in MRS
//    Only supported datatype is uint32_t and only supported width it 24 bit.

struct Selector
{
	int32_t databits;
	uint32_t span;
} const kSelectors[16] = {
	{  0, 1 },
	{ -4, 1 },
	{ -2, 1 }, { -2, 2 },
	{ -1, 1 }, { -1, 2 }, { -1, 4 },
	{  0, 1 }, {  0, 2 }, {  0, 4 },
	{  1, 1 }, {  1, 2 }, {  1, 4 },
	{  2, 1 }, {  2, 2 },
	{  4, 1 }
};
 
// store ints of at most 24 bits, should be enough.
const uint32_t kStartWidth = 8, kMaxWidth = 24;
 
inline uint32_t bitWidth(uint32_t v)
{
	uint32_t result = 0;
	while (v > 0)
	{
		v >>= 1;
		++result;
	}
	return result;
}

void CompressSimpleArraySelector(OBitStream& inBits, const std::vector<uint32_t>& inArray)
{
	int32_t width = kStartWidth;

	int32_t bn[4];
	uint32_t dv[4];
	uint32_t bc = 0;
	auto a = inArray.begin(), e = inArray.end();

	while (a != e or bc > 0)
	{
		while (bc < 4 and a != e)
		{
			dv[bc] = *a++;
			bn[bc] = bitWidth(dv[bc]);
			++bc;
		}

		uint32_t s = 0;
		int32_t c = bn[0] - kMaxWidth;

		for (uint32_t i = 1; i < 16; ++i)
		{
			if (kSelectors[i].span > bc)
				continue;

			int32_t w = width + kSelectors[i].databits;

			if (static_cast<uint32_t>(w) > kMaxWidth)
				continue;

			bool fits = true;
			int32_t waste = 0;

			switch (kSelectors[i].span)
			{
				case 4:
					fits = fits and bn[3] <= w;
					waste += w - bn[3];
					[[fallthrough]];
				case 3:
					fits = fits and bn[2] <= w;
					waste += w - bn[2];
					[[fallthrough]];
				case 2:
					fits = fits and bn[1] <= w;
					waste += w - bn[1];
					[[fallthrough]];
				case 1:
					fits = fits and bn[0] <= w;
					waste += w - bn[0];
			}

			if (fits == false)
				continue;

			int32_t n = (kSelectors[i].span - 1) * 4 - waste;

			if (n > c)
			{
				s = i;
				c = n;
			}
		}

		if (s == 0)
			width = kMaxWidth;
		else
			width += kSelectors[s].databits;

		uint32_t n = kSelectors[s].span;

		inBits.write(s, 4);

		if (width > 0)
		{
			for (uint32_t i = 0; i < n; ++i)
				inBits.write(dv[i], width);
		}

		bc -= n;

		if (bc > 0)
		{
			for (uint32_t i = 0; i < (4 - n); ++i)
			{
				bn[i] = bn[i + n];
				dv[i] = dv[i + n];
			}
		}
	}
}

void DecompressSimpleArraySelector(IBitStream& inBits, std::vector<uint32_t>& outArray)
{
	uint32_t width = kStartWidth;
	uint32_t span = 0;
	
	// The array should be initilialized to the expected size!
	auto size = outArray.size();
	auto a = outArray.begin();

	while (size-- > 0)
	{
		if (span == 0)
		{
			uint32_t selector = inBits.read(4);
			span = kSelectors[selector].span;

			if (selector == 0)
				width = kMaxWidth;
			else
				width += kSelectors[selector].databits;
		}

		if (width > 0)
			*a++ = inBits.read(width);
		else
			*a++ = 0;

		--span;
	}
}

// --------------------------------------------------------------------

enum class SecStrType : char
{
	helix = 'H',
	strand = 'E',
	other = '.',
	cis = 'c',
	prepro = 'p'
};

std::ostream& operator<<(std::ostream& os, SecStrType ss)
{
	switch (ss)
	{
		case SecStrType::helix: os << "helix"; break;
		case SecStrType::strand: os << "strand"; break;
		case SecStrType::other: os << "other"; break;
		case SecStrType::cis: os << "cis"; break;
		case SecStrType::prepro: os << "prepro"; break;
	}

	return os;
}

// --------------------------------------------------------------------
// The header for the data blocks as written in de resource

struct StoredData
{
	char		aa[3];
	SecStrType	ss;
	float		mean, mean_vs_random, sd, sd_vs_random, binSpacing;
	uint32_t	offset;			// offset into compressed data area
};

class Data
{
	friend class DataTable;

  public:
	Data(Data&& d)
		: aa(d.aa), ss(d.ss), mean(d.mean), sd(d.sd)
		, mean_vs_random(d.mean_vs_random), sd_vs_random(d.sd_vs_random)
		, binSpacing(d.binSpacing), counts(move(d.counts))
		, dim(d.dim), d2(d.d2)
	{
	}

	Data(const Data&) = delete;
	Data& operator=(const Data&) = delete;

	Data(const char* type, const std::string& aa, SecStrType ss, std::istream& is);
	Data(bool torsion, const StoredData& data, const uint8_t* bits);

	void store(StoredData& data, std::vector<uint8_t>& databits);

	float interpolatedCount(float phi, float a2) const;
	float zscore(float a1, float a2) const
	{
		return (interpolatedCount(a1, a2) - mean) / sd;
	}

	void dump() const
	{
		for (size_t i = 0; i < counts.size(); ++i)
		{
			float a1, a2;
			std::tie(a1, a2) = angles(i);
			std::cout << a1 << ' ' << a2 << ' ' << counts[i] << std::endl;
		}
	}

  private:

	std::string aa;
	SecStrType ss;
	float mean, sd, mean_vs_random, sd_vs_random;
	float binSpacing;
	std::vector<uint32_t> counts;

	// calculated
	size_t dim;
	bool d2;
	
	float count(size_t a1Ix, size_t a2Ix) const
	{
		a1Ix %= dim;
		a2Ix %= dim;
		return d2 ? counts.at(a1Ix * dim + a2Ix) : counts.at(a1Ix);
	}

	size_t index(float a1, float a2 = 0) const
	{
		size_t x = 0, y = 0;

		if (d2)
		{
			x = static_cast<size_t>((a1 + 180) / binSpacing);
			y = static_cast<size_t>((a2 + 180) / binSpacing);
		}
		else
			y = static_cast<size_t>((a1 + 180) / binSpacing);

		return x * (360 / binSpacing) + y;
	}

	std::tuple<float,float> angles(size_t index) const
	{
		size_t x = index / dim;
		size_t y = index % dim;

		return std::make_tuple(x * binSpacing - 180, y * binSpacing - 180);
	}
};

Data::Data(const char* type, const std::string& aa, SecStrType ss, std::istream& is)
	: aa(aa), ss(ss)
{
	// example:
	// 14400 bins, aver 19.2878, sd 15.4453, binspacing 3
	// torsion vs random: 2.0553 2.8287
	static const std::regex
		kRX1(R"((\d+) bins, aver ([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?), sd ([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?), binspacing ([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?))"),
		kRX2(R"((torsion|rama) vs random: ([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?) ([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?))");

	std::string line;
	getline(is, line);

	d2 = strcmp(type, "torsion") != 0 or std::set<std::string>{"CYS", "SER", "THR", "VAL"}.count(aa) == 0;

	std::smatch m;
	if (not std::regex_match(line, m, kRX1))
		throw std::runtime_error("Invalid file");

	size_t nBins = stoi(m[1]);
	mean = stof(m[2]);
	sd = stof(m[3]);
	binSpacing = stof(m[4]);

	dim = static_cast<size_t>(360 / binSpacing);
	if ((d2 and nBins != dim * dim) or (not d2 and nBins != dim))
		throw std::runtime_error("Unexpected number of bins");

	counts.resize(nBins);

	getline(is, line);

	if (not std::regex_match(line, m, kRX2) or m[1] != type)
		throw std::runtime_error("Invalid file");
	
	mean_vs_random = stof(m[2]);
	sd_vs_random = stof(m[3]);

	for (size_t i = 0; i < nBins; ++i)
	{
		float a1, a2;
		uint32_t count;

		if (d2)
			is >> a1 >> a2 >> count;
		else
			is >> a1 >> count;

		if (is.eof())
			throw std::runtime_error("truncated file?");
		
		counts.at(index(a1, a2)) = count;
	}
}

Data::Data(bool torsion, const StoredData& data, const uint8_t* databits)
{
	aa.assign(data.aa, data.aa + 3);
	ss = data.ss;
	mean = data.mean;
	mean_vs_random = data.mean_vs_random;
	sd = data.sd;
	sd_vs_random = data.sd_vs_random;
	binSpacing = data.binSpacing;
	
	d2 = not torsion or std::set<std::string>{"CYS", "SER", "THR", "VAL"}.count(aa) == 0;

	size_t nBins = static_cast<size_t>(360 / binSpacing);
	dim = nBins;

	if (d2)
		nBins *= nBins;

	counts.insert(counts.begin(), nBins, 0);

	IBitStream bits(databits + data.offset);
	DecompressSimpleArraySelector(bits, counts);
}

void Data::store(StoredData& data, std::vector<uint8_t>& databits)
{
	assert(aa.length() == 3);
	copy(aa.begin(), aa.end(), data.aa);
	data.ss = ss;
	data.mean = mean;
	data.sd = sd;
	data.mean_vs_random = mean_vs_random;
	data.sd_vs_random = sd_vs_random;
	data.offset = databits.size();
	data.binSpacing = binSpacing;
	
	OBitStream bits(databits);
	CompressSimpleArraySelector(bits, counts);
	bits.sync();
}

float Data::interpolatedCount(float a1, float a2) const
{
	const size_t N = dim;

	float result;

	if (d2)
	{
		size_t a1FloorIx = static_cast<size_t>(N * (a1 + 180) / 360);
		size_t a2FloorIx = static_cast<size_t>(N * (a2 + 180) / 360);
		
		size_t a1CeilIx = (a1FloorIx + 1);
		size_t a2CeilIx = (a2FloorIx + 1);

		float a1FloorAngle = (a1FloorIx * 360.0f) / N - 180;
		float a2FloorAngle = (a2FloorIx * 360.0f) / N - 180;

		float a1CeilAngle = (a1CeilIx * 360.0f) / N - 180;
		float a2CeilAngle = (a2CeilIx * 360.0f) / N - 180;

		float a1Factor = a1CeilIx > a1FloorIx ? (a1 - a1FloorAngle) / (a1CeilAngle - a1FloorAngle) : 1;
		float a2Factor = a2CeilIx > a2FloorIx ? (a2 - a2FloorAngle) / (a2CeilAngle - a2FloorAngle) : 1;

		float c1 = count(a1FloorIx, a2FloorIx) + (count(a1CeilIx, a2FloorIx) - count(a1FloorIx, a2FloorIx)) * a1Factor;
		float c2 = count(a1FloorIx, a2CeilIx) + (count(a1CeilIx, a2CeilIx) - count(a1FloorIx, a2CeilIx)) * a1Factor;
		
		result = c1 + (c2 - c1) * a2Factor;
	}
	else
	{
		size_t a1FloorIx = static_cast<size_t>(N * (a1 + 180) / 360);
		size_t a1CeilIx = (a1FloorIx + 1);

		float a1FloorAngle = (a1FloorIx * 360.0f) / N - 180;
		float a1CeilAngle = (a1CeilIx * 360.0f) / N - 180;

		float a1Factor = a1CeilIx > a1FloorIx ? (a1 - a1FloorAngle) / (a1CeilAngle - a1FloorAngle) : 1;

		result = count(a1FloorIx, 0) + (count(a1CeilIx, 0) - count(a1FloorIx, 0)) * a1Factor;
	}

	return result;
}

void buildDataFile(const fs::path &dir)
{
	using namespace std::literals;

	// first read the global mean and sd

	float mean_torsion, sd_torsion, mean_ramachandran, sd_ramachandran;

	std::ifstream in(dir / "zscores_proteins.txt");
	std::string line;
	const std::regex krx(R"((Rama|Rota): average ([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?), sd ([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?))");

	while (getline(in, line))
	{
		std::smatch m;
		if (not std::regex_match(line, m, krx))
			continue;

		if (m[1] == "Rama")
		{
			mean_ramachandran = stof(m[2]);
			sd_ramachandran = stof(m[3]);
		}
		else
		{
			mean_torsion = stof(m[2]);
			sd_torsion = stof(m[3]);
		}
	}

	std::vector<StoredData> data;
	std::vector<uint8_t> bits;

	// first ramachandran counts
	for (const char *aa: {
		 "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
	})
	{
		for (std::pair<SecStrType, const char*> ss: {
				std::make_pair(SecStrType::helix, "helix"),
				std::make_pair(SecStrType::strand, "strand"),
				std::make_pair(SecStrType::other, "other") })
		{
			auto p = dir / ("rama_count_"s + ss.second + '_' + aa + ".txt");
			if (not fs::exists(p))
				continue;
			
			std::ifstream f(p);
			Data d("rama", aa, ss.first, f);

			StoredData sd = { };
			d.store(sd, bits);
			data.push_back(sd);
		}
	}

	for (std::tuple<SecStrType, const char*, const char*> ss: {
			std::make_tuple(SecStrType::cis, "PRO", "cis_PRO"),
			std::make_tuple(SecStrType::prepro, "***", "prepro_all_noGIV"),
			std::make_tuple(SecStrType::prepro, "GLY", "prepro_GLY"),
			std::make_tuple(SecStrType::prepro, "IV_", "prepro_ILEVAL") })
	{
		auto p = dir / ("rama_count_"s + std::get<2>(ss) + ".txt");
		if (not fs::exists(p))
			continue;
		
		std::ifstream f(p);
		Data d("rama", std::get<1>(ss), std::get<0>(ss), f);

		StoredData sd = { };
		d.store(sd, bits);
		data.push_back(sd);
	}

	data.push_back({});

	if (fs::exists("rama-data.bin"))
		fs::remove("rama-data.bin");
	std::ofstream out("rama-data.bin", std::ios::binary);
	if (not out.is_open())
		throw std::runtime_error("Could not create rama-data.bin file");
	out.write(reinterpret_cast<char*>(&mean_ramachandran), sizeof(mean_ramachandran));
	out.write(reinterpret_cast<char*>(&sd_ramachandran), sizeof(sd_ramachandran));
	out.write(reinterpret_cast<char*>(data.data()), data.size() * sizeof(StoredData));
	out.write(reinterpret_cast<char*>(bits.data()), bits.size());
	out.close();

	data.clear();
	bits.clear();

	// next torsion counts
	for (const char *aa: {
		 "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
	})
	{
		for (std::pair<SecStrType, const char*> ss: {
				std::make_pair(SecStrType::helix, "helix"),
				std::make_pair(SecStrType::strand, "strand"),
				std::make_pair(SecStrType::other, "other") })
		{
			auto p = dir / ("torsion_count_"s + ss.second + '_' + aa + ".txt");
			if (not fs::exists(p))
				continue;
			
			std::ifstream f(p);
			Data d("torsion", aa, ss.first, f);

			StoredData sd = { };
			d.store(sd, bits);
			data.push_back(sd);
		}
	}

	data.push_back({});

	if (fs::exists("torsion-data.bin"))
		fs::remove("torsion-data.bin");
	out.open("torsion-data.bin", std::ios::binary);
	if (not out.is_open())
		throw std::runtime_error("Could not create torsion-data.bin file");
	out.write(reinterpret_cast<char*>(&mean_torsion), sizeof(mean_torsion));
	out.write(reinterpret_cast<char*>(&sd_torsion), sizeof(sd_torsion));
	out.write(reinterpret_cast<char*>(data.data()), data.size() * sizeof(StoredData));
	out.write(reinterpret_cast<char*>(bits.data()), bits.size());
	out.close();
}

// --------------------------------------------------------------------

class DataTable
{
  public:

	static DataTable& instance()
	{
		static DataTable sInstance;
		return sInstance;
	}

	const Data& loadTorsionData(const std::string& aa, SecStrType ss) const;
	const Data& loadRamachandranData(const std::string& aa, SecStrType ss) const;

	float mean_torsion() const		{ return m_mean_torsion; }
	float sd_torsion() const		{ return m_sd_torsion; }
	float mean_ramachandran() const	{ return m_mean_ramachandran; }
	float sd_ramachandran() const	{ return m_sd_ramachandran; }

  private:

	DataTable(const DataTable& ) = delete;
	DataTable& operator=(const DataTable&) = delete;

	DataTable();

	void load(const char* name, std::vector<Data>& table, float& mean, float& sd);

	std::vector<Data>	m_torsion, m_ramachandran;

	float m_mean_torsion, m_sd_torsion, m_mean_ramachandran, m_sd_ramachandran;
};

DataTable::DataTable()
{
	load("torsion-data.bin", m_torsion, m_mean_torsion, m_sd_torsion);
	load("rama-data.bin", m_ramachandran, m_mean_ramachandran, m_sd_ramachandran);
}

const Data& DataTable::loadTorsionData(const std::string& aa, SecStrType ss) const
{
	auto i = find_if(m_torsion.begin(), m_torsion.end(), [aa, ss](auto& d)
		{ return d.aa == aa and d.ss == ss; });
	if (i == m_torsion.end())
		throw std::runtime_error("Data missing for aa = " + aa + " and ss = '" + static_cast<char>(ss) + '\'');

	return *i;
}

const Data& DataTable::loadRamachandranData(const std::string& aa, SecStrType ss) const
{
	std::vector<Data>::const_iterator i;

	switch (ss)
	{
		case SecStrType::cis:
			i = find_if(m_ramachandran.begin(), m_ramachandran.end(), [](auto& d) { return d.ss == SecStrType::cis and d.aa == "PRO"; });
			break;

		case SecStrType::prepro:
			i = find_if(m_ramachandran.begin(), m_ramachandran.end(), [aa](auto& d) {
					bool result = false;
					if (d.ss == SecStrType::prepro)
					{
						if (aa == "GLY")
							result = d.aa == "GLY";
						else if (aa == "ILE" or aa == "VAL")
							result = d.aa == "IV_";
						else
							result = d.aa == "***";
					}
					return result;
				});
			break;

		default:
			i = find_if(m_ramachandran.begin(), m_ramachandran.end(), [aa, ss](auto& d)
				{ return d.aa == aa and d.ss == ss; });
			break;
	}

	if (i == m_ramachandran.end())
		throw std::runtime_error("Data missing for aa= " + aa + " and ss = '" + static_cast<char>(ss) + '\'');
	
	return *i;
}

void DataTable::load(const char* name, std::vector<Data>& table, float& mean, float& sd)
{
	using namespace std::literals;

	auto rfd = cif::load_resource(name);

	if (not rfd)
		throw std::runtime_error("Missing resource "s + name);

	rfd->seekg(0, rfd->end);
	auto size = rfd->tellg();
	rfd->seekg(0, rfd->beg);

	const float* fv = new float[size / sizeof(float) + 1];
	rfd->read(reinterpret_cast<char*>(const_cast<float*>(fv)), size);

	mean = fv[0];
	sd = fv[1];

	const StoredData* data = reinterpret_cast<const StoredData*>(fv + 2);
	size_t ix = 0;
	while (data[ix].aa[0] != 0)
		++ix;
	
	size_t n = ix;
	const uint8_t* bits = reinterpret_cast<const uint8_t*>(fv + 2) + (n + 1) * sizeof(StoredData);

	for (ix = 0; ix < n; ++ix)
		table.emplace_back(strcmp(name, "torsion-data.bin") == 0, data[ix], bits);
}

// --------------------------------------------------------------------

float jackknife(const std::vector<float>& zScorePerResidue)
{
	// jackknife variance estimate, see: https://en.wikipedia.org/wiki/Jackknife_resampling

	const size_t N = zScorePerResidue.size();
	double zScoreSum = accumulate(zScorePerResidue.begin(), zScorePerResidue.end(), 0.0);
	std::vector<double> scores(N);

	DataTable& tbl = DataTable::instance();

	double scoreSum = 0;
	for (size_t i = 0; i < zScorePerResidue.size(); ++i)
	{
		double score = (zScoreSum - zScorePerResidue[i]) / (N - 1);
		score = (score - tbl.mean_ramachandran()) / tbl.sd_ramachandran();
		scores[i] = score;
		scoreSum += score;
	}

	double avg = scoreSum / N;
	double sumD = accumulate(scores.begin(), scores.end(), 0.0, [avg](double a, double z) { return a + (z - avg) * (z - avg); });

	return sqrt((N - 1) * sumD / N);
}

// --------------------------------------------------------------------

json calculateZScores(const cif::datablock &db, int model_nr)
{
	dssp::DSSP dssp(db, model_nr, 3, false);
	auto& tbl = DataTable::instance();

	double ramaZScoreSum = 0;
	size_t ramaZScoreCount = 0;
	double torsZScoreSum = 0;
	size_t torsZScoreCount = 0;

	json residues;
	std::vector<float> ramaZScorePerResidue, torsZScorePerResidue;

	for (auto res : dssp)
	{
		auto phi = res.phi();
		auto psi = res.psi();

		if (phi == 360 or psi == 360)
			continue;

		std::string aa = res.compound_id();

		json residue = {
			{ "asymID", res.asym_id() },
			{ "seqID", res.seq_id() },
			{ "compID", aa },
			{ "pdb", {
				{ "strandID", res.pdb_strand_id() },
				{ "seqNum", res.pdb_seq_num() },
				{ "compID", aa },
				{ "insCode", res.pdb_ins_code() }
			}}
		};

		// remap some common modified amino acids
		if (aa == "MSE")
		{
			if (cif::VERBOSE > 1)
				std::cerr << "Replacing MSE with MET" << std::endl;
			aa = "MET";
		}
		else if (aa == "HYP")
		{
			if (cif::VERBOSE > 1)
				std::cerr << "Replacing HYP with PRO" << std::endl;

			aa = "PRO";
		}
		else if (aa == "ASX")
		{
			if (cif::VERBOSE > 1)
				std::cerr << "Replacing ASX with ASP" << std::endl;

			aa = "ASP";
		}
		else if (aa == "GLX")
		{
			if (cif::VERBOSE > 1)
				std::cerr << "Replacing GLX with GLU" << std::endl;

			aa = "GLU";
		}
		else if (res.compound_letter() == 'X')
		{
			if (cif::VERBOSE)
				std::cerr << "Replacing " << aa << " with ALA" << std::endl;

			aa = "ALA";
		}
		
		SecStrType tors_ss, rama_ss;

		switch (res.type())
		{
			case dssp::structure_type::Alphahelix:	tors_ss = SecStrType::helix; break;
			case dssp::structure_type::Strand:		tors_ss = SecStrType::strand; break;
			default:					tors_ss = SecStrType::other; break;
		}

		if (res.is_pre_pro())
			rama_ss = SecStrType::prepro;
		else if (aa == "PRO" && res.is_cis())
			rama_ss = SecStrType::cis;
		else
			rama_ss = tors_ss;

// TODO: #pragma warning "todo" (but the question now is, what is here to do???)
		auto& rd = tbl.loadRamachandranData(aa, rama_ss);

		auto zr = rd.zscore(phi, psi);

		residue["ramachandran"] =
		{
			{ "ss-type", boost::lexical_cast<std::string>(rama_ss) },
			{ "z-score", zr }
		};

		ramaZScorePerResidue.push_back(zr);

		ramaZScoreSum += zr;
		++ramaZScoreCount;

		try
		{
			float zt = nan("1");

			auto chiCount = res.nr_of_chis();
			if (chiCount)
			{
				float chi1 = res.chi(0);
				float chi2 = chiCount > 1 ? res.chi(1) : 0;

				auto& td = tbl.loadTorsionData(aa, tors_ss);

				zt = td.zscore(chi1, chi2);

				torsZScoreSum += zt;
				++torsZScoreCount;

				torsZScorePerResidue.push_back(zt);

				residue["torsion"] =
				{
					{ "ss-type", boost::lexical_cast<std::string>(tors_ss) },
					{ "z-score", zt }
				};
			}
		}
		catch (const std::exception& e)
		{
			if (cif::VERBOSE)
				std::cerr << e.what() << '\n';
		}

		residues.push_back(residue);
	}

	float ramaVsRand = ramaZScoreSum / ramaZScoreCount;
	float torsVsRand = torsZScoreSum / torsZScoreCount;

	float jackknifeRama = jackknife(ramaZScorePerResidue);
	float jackknifeTors = jackknife(torsZScorePerResidue);

	return {
		{ "ramachandran-z", ((ramaVsRand - tbl.mean_ramachandran()) / tbl.sd_ramachandran()) },
		{ "ramachandran-jackknife-sd", jackknifeRama },
		{ "torsion-z", ((torsVsRand - tbl.mean_torsion()) / tbl.sd_torsion()) },
		{ "torsion-jackknife-sd", jackknifeTors },
		{ "residues", residues },
	};
}

// --------------------------------------------------------------------

json tortoize_calculate(const fs::path &xyzin)
{
	json data{
		{ "software",
			{
				{ "name", "tortoize" },
				{ "version", kVersionNumber },
				{ "reference", "Sobolev et al. A Global Ramachandran Score Identifies Protein Structures with Unlikely Stereochemistry, Structure (2020)" },
				{ "reference-doi", "https://doi.org/10.1016/j.str.2020.08.005" }
			}
		}
	};

	// --------------------------------------------------------------------

	cif::file f(xyzin);

	if (not f.is_valid())
		throw std::runtime_error("Invalid mmCIF file");

	std::set<uint32_t> models;
	for (auto r: f.front()["atom_site"])
	{
		if (not r["pdbx_PDB_model_num"].empty())
			models.insert(r["pdbx_PDB_model_num"].as<uint32_t>());
	}

	if (models.empty())
		models.insert(0);

	for (auto model: models)
		data["model"][std::to_string(model)] = calculateZScores(f.front(), model);

	return data;
}