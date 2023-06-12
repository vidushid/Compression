#include "./lz_diff.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include "./zstd/lib/zstd.h"

using namespace std;

// *******************************************************************************************
CLZDiffBase::CLZDiffBase(const uint32_t _min_match_len)
{
	min_match_len = _min_match_len;
	key_len = min_match_len - hashing_step + 1u;
	short_ht_ver = false;
	ht_mask = 0;
	ht_size = 0;
	index_ready = false;
}

// *******************************************************************************************
CLZDiffBase::~CLZDiffBase()
{
}

// *******************************************************************************************
bool CLZDiffBase::SetMinMatchLen(const uint32_t _min_match_len)
{
	if (!reference.empty() || index_ready)
		return false;

	min_match_len = _min_match_len;
	key_len = min_match_len - hashing_step + 1u;

	return true;
}

//#define USE_REVCOMP_REFERENCE
//#define USE_REV_REFERENCE
// *******************************************************************************************
void CLZDiffBase::prepare_gen(const contig_t& _reference)
{
	reference = _reference;

	reference.resize(reference.size() + key_len, invalid_symbol);

#ifdef USE_REVCOMP_REFERENCE
	reference.reserve(2 * reference.size());

	for (auto p = _reference.rbegin(); p != _reference.rend(); ++p)
		if (*p < 4)
			reference.emplace_back(3 - *p);
		else
			reference.emplace_back(*p);

	for (int i = 0; i < key_len; ++i)
		reference.emplace_back(invalid_symbol);
#endif

#ifdef USE_REV_REFERENCE
	reference.reserve(2 * reference.size());

	for (auto p = _reference.rbegin(); p != _reference.rend(); ++p)
		reference.emplace_back(*p);

	for (int i = 0; i < key_len; ++i)
		reference.emplace_back(invalid_symbol);
#endif

	reference.shrink_to_fit();
}

// *******************************************************************************************
void CLZDiffBase::prepare_index()
{
	ht_size = 0;

	uint32_t no_prev_valid = 0;

#ifdef USE_SPARSE_HT
	uint32_t cnt_mod = 0;
	uint32_t key_len_mod = key_len % hashing_step;

	for (auto c : reference)
	{
		if (c < 4)
			++no_prev_valid;
		else
			no_prev_valid = 0;

		if (++cnt_mod == hashing_step)
			cnt_mod = 0;

		if (cnt_mod == key_len_mod && no_prev_valid >= key_len)
			++ht_size;
	}
#else
	for (auto c : reference)
	{
		if (c < 4)
			++no_prev_valid;
		else
			no_prev_valid = 0;

		if (no_prev_valid >= key_len)
			++ht_size;
	}
#endif

	ht_size = (uint64_t)(ht_size / max_load_factor);

	while (ht_size & (ht_size - 1))
		ht_size &= ht_size - 1;

	ht_size <<= 1;
	ht_mask = ht_size - 1;

	if (short_ht_ver)
	{
		ht16.resize(ht_size, empty_key16);
		make_index16();
	}
	else
	{
		ht32.resize(ht_size, empty_key32);
		make_index32();
	}

	index_ready = true;
}

// *******************************************************************************************
void CLZDiffBase::Prepare(const contig_t& _reference)
{
	short_ht_ver = _reference.size() / hashing_step < 65535;

	prepare_gen(_reference);
}

// *******************************************************************************************
void CLZDiffBase::GetCodingCostVector(const contig_t& text, vector<uint32_t>& v_costs, const bool prefix_costs)
{
	if (!index_ready)
		prepare_index();

	v_costs.clear();
	v_costs.reserve(text.size());

	uint32_t text_size = (uint32_t)text.size();
	MurMur64Hash mmh;

	uint32_t i = 0;
	uint32_t pred_pos = 0;

	const uint8_t* text_ptr = text.data();

#ifdef USE_SPARSE_HT
	uint32_t no_prev_literals = 0;
#endif

	for (; i + key_len < text_size; )
	{
		uint64_t x = get_code(text_ptr);
		if (x == ~0ull)
		{
			uint32_t Nrun_len = get_Nrun_len(text_ptr, text_size - i);

			if (Nrun_len >= min_Nrun_len)
			{
				auto tc = coding_cost_Nrun(Nrun_len);
				if (prefix_costs)
				{
					v_costs.emplace_back(tc);
					v_costs.insert(v_costs.end(), Nrun_len - 1, 0);
				}
				else
				{
					v_costs.insert(v_costs.end(), Nrun_len - 1, 0);
					v_costs.emplace_back(tc);
				}

				text_ptr += Nrun_len;
				i += Nrun_len;
#ifdef USE_SPARSE_HT
				no_prev_literals = 0;
#endif
			}
			else
			{
				v_costs.emplace_back(1);
				++i;
				++pred_pos;
				++text_ptr;
#ifdef USE_SPARSE_HT
				++no_prev_literals;
#endif
			}

			continue;
		}

		uint64_t ht_pos = mmh(x) & ht_mask;

		uint32_t len_bck = 0;
		uint32_t len_fwd = 0;
		uint32_t match_pos = 0;
		uint32_t max_len = text_size - i;

		if (short_ht_ver ?
			!find_best_match16((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd) :
			!find_best_match32((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd))
		{
			v_costs.emplace_back(1);
			//			encode_literal(*text_ptr, encoded);
			++i;
			++text_ptr;
			++pred_pos;
#ifdef USE_SPARSE_HT
			++no_prev_literals;
#endif
			continue;
		}
		else
		{
#ifdef USE_SPARSE_HT
			if (len_bck)
			{
				for (uint32_t k = 0; k < len_bck; ++k)
					v_costs.pop_back();
				match_pos -= len_bck;
				pred_pos -= len_bck;
				text_ptr -= len_bck;
				i -= len_bck;
			}
#endif

			auto tc = coding_cost_match(match_pos, len_bck + len_fwd, pred_pos);

			if (prefix_costs)
			{
				v_costs.emplace_back(tc);
				v_costs.insert(v_costs.end(), len_bck + len_fwd - 1, 0);
			}
			else
			{
				v_costs.insert(v_costs.end(), len_bck + len_fwd - 1, 0);
				v_costs.emplace_back(tc);
			}

			pred_pos = match_pos + len_bck + len_fwd;
			i += len_bck + len_fwd;
			text_ptr += len_bck + len_fwd;

#ifdef USE_SPARSE_HT
			no_prev_literals = 0;
#endif
		}
	}

	for (; i < text_size; ++i)
		v_costs.emplace_back(1);
}

// *******************************************************************************************
bool CLZDiffBase::find_best_match16(uint32_t ht_pos, const uint8_t* s, const uint32_t max_len, const uint32_t no_prev_literals,
	uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd)
{
	len_fwd = 0;
	len_bck = 0;

	const uint8_t* ref_ptr = reference.data();

	for (uint32_t i = 0; i < max_no_tries; ++i)
	{
		if (ht16[ht_pos] == empty_key16)
			break;

		uint32_t f_len = 0;
		uint32_t h_pos = ht16[ht_pos] * hashing_step;

		const uint8_t* p = ref_ptr + h_pos;

		for (; f_len < max_len; ++f_len)
			if (s[f_len] != p[f_len])
				break;

		uint32_t b_len = 0;
		for (; b_len < min(no_prev_literals, h_pos); ++b_len)
			if (*(s - b_len - 1) != *(p - b_len - 1))
				break;

		if (b_len + f_len > len_bck + len_fwd)
		{
			len_bck = b_len;
			len_fwd = f_len;
			ref_pos = h_pos;
		}

		ht_pos = (uint32_t)((ht_pos + 1u) & ht_mask);
	}

	return len_bck + len_fwd >= min_match_len;
}

// *******************************************************************************************
bool CLZDiffBase::find_best_match32(uint32_t ht_pos, const uint8_t* s, const uint32_t max_len, const uint32_t no_prev_literals,
	uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd)
{
	len_fwd = 0;
	len_bck = 0;

	const uint8_t* ref_ptr = reference.data();

	for (uint32_t i = 0; i < max_no_tries; ++i)
	{
		if (ht32[ht_pos] == empty_key32)
			break;

		uint32_t f_len = 0;
		uint32_t h_pos = ht32[ht_pos] * hashing_step;
		const uint8_t* p = ref_ptr + h_pos;

		for (; f_len < max_len; ++f_len)
			if (s[f_len] != p[f_len])
				break;

		uint32_t b_len = 0;
		for (; b_len < min(no_prev_literals, h_pos); ++b_len)
			if (*(s - b_len - 1) != *(p - b_len - 1))
				break;

		if (b_len + f_len > len_bck + len_fwd)
		{
			len_bck = b_len;
			len_fwd = f_len;
			ref_pos = h_pos;
		}

		ht_pos = (uint32_t) ((ht_pos + 1u) & ht_mask);
	}

	return len_bck + len_fwd >= min_match_len;
}

// *******************************************************************************************
void CLZDiffBase::append_int(contig_t& text, int64_t x)
{
	if (x == 0)
	{
		text.emplace_back('0');
		return;
	}

	if (x < 0)
	{
		text.push_back('-');
		x = -x;
	}

	contig_t tmp;

	for (; x; x /= 10)
		tmp.emplace_back((uint8_t) ('0' + (x % 10)));

	text.insert(text.end(), tmp.rbegin(), tmp.rend());
}

// *******************************************************************************************
void CLZDiffBase::read_int(contig_t::const_iterator& p, int64_t& x)
{
	bool is_neg = false;
	x = 0;

	if (*p == '-')
	{
		is_neg = true;
		++p;
	}

	while(*p >= '0' && *p <= '9')
		x = x * 10 + (int64_t)(*p++ - '0');

	if (is_neg)
		x = -x;
}

// *******************************************************************************************
void CLZDiffBase::encode_literal(const uint8_t c, contig_t& encoded)
{
	encoded.push_back('A' + c);
}

// *******************************************************************************************
void CLZDiffBase::encode_literal_diff(const uint8_t c, const uint8_t r, contig_t& encoded)
{
	if (r == 0 || (r > 3 || c > 3))
		encoded.push_back(c);
	else
	{
		if (c < r)
			encoded.push_back(3 - c);
		else
			encoded.push_back(c - r);
	}

}

// *******************************************************************************************
void CLZDiffBase::encode_Nrun(const uint32_t len, contig_t &encoded)
{
	encoded.emplace_back(N_run_starter_code);		// N-run start marker
	append_int(encoded, len - min_Nrun_len);
	encoded.emplace_back(N_code);					// N-run stop marker
}

// *******************************************************************************************
uint32_t CLZDiffBase::coding_cost_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos) const
{
	uint32_t r;
	int dif_pos = (int)ref_pos - (int)pred_pos;

	if (dif_pos >= 0)
		r = int_len((uint32_t)dif_pos);
	else
		r = int_len((uint32_t)-dif_pos) + 1;

	r += int_len(len - min_match_len) + 2;

	return r;
}

// *******************************************************************************************
uint32_t CLZDiffBase::coding_cost_Nrun(const uint32_t len) const
{
	return 2 + int_len(len - min_Nrun_len);
}

// *******************************************************************************************
uint32_t CLZDiffBase::get_Nrun_len(const uint8_t* s, const uint32_t max_len) const
{
	if (*s != N_code || *(s + 1) != N_code || *(s + 2) != N_code)
		return 0;

	uint32_t len;
	for (len = 3; len < max_len && *(s + len) == N_code; ++len)
		;

	return len;
}

// *******************************************************************************************
bool CLZDiffBase::is_literal(const contig_t::const_iterator& p) const
{
	return (*p >= 'A' && *p <= 'A' + 20) || (*p == '!');
}

// *******************************************************************************************
bool CLZDiffBase::is_Nrun(const contig_t::const_iterator& p) const
{
	return *p == N_run_starter_code;
}

// *******************************************************************************************
void CLZDiffBase::decode_literal(contig_t::const_iterator& p, uint8_t &c)
{
	if (*p == '!')
	{
		c = '!';
		++p;
	}
	else
		c = *p++ - 'A';
}

// *******************************************************************************************
void CLZDiffBase::decode_Nrun(contig_t::const_iterator& p, uint32_t& len)
{
	int64_t raw_len;

	++p;		// prefix
	read_int(p, raw_len);
	++p;		// suffix

	len = (uint32_t) (raw_len + min_Nrun_len);
}

// *******************************************************************************************
uint64_t CLZDiffBase::get_code(const uint8_t* s) const
{
	uint64_t x = 0;

	for (uint32_t i = 0; i < key_len; ++i, ++s)
	{
		if (*s > 3)
			return ~0ull;
		x = (x << 2) + (uint64_t) *s;
	}

	return x;
}

// *******************************************************************************************
void CLZDiffBase::make_index16()
{
	uint32_t ref_size = (uint32_t)reference.size();
	MurMur64Hash mmh;

	const uint8_t* ptr = reference.data();

#ifdef USE_SPARSE_HT
	for (uint32_t i = 0; i + key_len < ref_size; i += hashing_step, ptr += hashing_step)
#else
	for (uint32_t i = 0; i + key_len < ref_size; ++i, ++ptr)
#endif
	{
		uint64_t x = get_code(ptr);
		if (x == ~0ull)
			continue;
		uint64_t pos = mmh(x) & ht_mask;

		for (uint32_t j = 0; j < max_no_tries; ++j)
			if (ht16[(pos + j) & ht_mask] == empty_key16)
			{
				ht16[(pos + j) & ht_mask] = i / hashing_step;
				break;
			}
	}
}

// *******************************************************************************************
void CLZDiffBase::make_index32()
{
	uint32_t ref_size = (uint32_t)reference.size();
	MurMur64Hash mmh;

	const uint8_t* ptr = reference.data();

#ifdef USE_SPARSE_HT
	for (uint32_t i = 0; i + key_len < ref_size; i += hashing_step, ptr += hashing_step)
#else
		for (uint32_t i = 0; i + key_len < ref_size; ++i, ++ptr)
#endif
	{
		uint64_t x = get_code(ptr);
		if (x == ~0ull)
			continue;
		uint64_t pos = mmh(x) & ht_mask;

		for (uint32_t j = 0; j < max_no_tries; ++j)
			if (ht32[(pos + j) & ht_mask] == empty_key32)
			{
				ht32[(pos + j) & ht_mask] = i / hashing_step;
				break;
			}
	}
}

// *******************************************************************************************
void CLZDiffBase::GetReference(contig_t& s)
{
	if (reference.empty())
		s.clear();
	else
		s.assign(reference.begin(), reference.begin() + (reference.size() - key_len));
}


// *******************************************************************************************
//
// *******************************************************************************************
void CLZDiff_V1::encode_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos, contig_t& encoded)
{
	int dif_pos = (int)ref_pos - (int)pred_pos;

	append_int(encoded, dif_pos);
	encoded.emplace_back(',');
	append_int(encoded, len - min_match_len);

	encoded.emplace_back('.');
}

// *******************************************************************************************
void CLZDiff_V1::decode_match(contig_t::const_iterator& p, uint32_t& ref_pos, uint32_t& len, uint32_t& pred_pos)
{
	int64_t raw_pos;
	int64_t raw_len;

	read_int(p, raw_pos);
	++p;		// ','

	ref_pos = (uint32_t)(raw_pos + (int64_t)pred_pos);

	if (*p == '.')
		len = ~0u;
	else
	{
		read_int(p, raw_len);
		len = (uint32_t)(raw_len + min_match_len);
	}

	++p;		// '.'
}

// *******************************************************************************************
void CLZDiff_V1::Encode(const contig_t& text, contig_t& encoded)
{
	if (!index_ready)
		prepare_index();

	uint32_t text_size = (uint32_t)text.size();

	encoded.clear();

#ifdef IMPROVED_LZ_ENCODING
	if (text_size == reference.size() - key_len)
		if (equal(text.begin(), text.end(), reference.begin()))
			return;														// equal sequences
#endif

	MurMur64Hash mmh;

	uint32_t i = 0;
	uint32_t pred_pos = 0;


	const uint8_t* text_ptr = text.data();

#ifdef USE_SPARSE_HT
//	const uint32_t thr_literal_len = 1000;

	uint32_t no_prev_literals = 0;
#endif

	for (; i + key_len < text_size; )
	{
		uint64_t x = get_code(text_ptr);

		if (x == ~0ull)
		{
			uint32_t Nrun_len = get_Nrun_len(text_ptr, text_size - i);

			if (Nrun_len >= min_Nrun_len)
			{
				encode_Nrun(Nrun_len, encoded);
				text_ptr += Nrun_len;
				i += Nrun_len;
#ifdef USE_SPARSE_HT
				/*				if (no_prev_literals > thr_literal_len)
									cerr << "literal len: " + to_string(no_prev_literals) + "    \n";*/
				no_prev_literals = 0;
#endif
			}
			else
			{
				encode_literal(*text_ptr, encoded);

				++i;
				++pred_pos;
				++text_ptr;
#ifdef USE_SPARSE_HT
				++no_prev_literals;
#endif
			}

			continue;
		}

		uint64_t ht_pos = mmh(x) & ht_mask;

		uint32_t len_bck = 0;
		uint32_t len_fwd = 0;
		uint32_t match_pos = 0;
		uint32_t max_len = text_size - i;

		if (short_ht_ver ?
			!find_best_match16((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd) :
			!find_best_match32((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd))
		{
			encode_literal(*text_ptr, encoded);

			++i;
			++text_ptr;
			++pred_pos;
#ifdef USE_SPARSE_HT
			++no_prev_literals;
#endif
			continue;
		}
		else
		{
#ifdef USE_SPARSE_HT
			if (len_bck)
			{
				for (uint32_t k = 0; k < len_bck; ++k)
					encoded.pop_back();
				match_pos -= len_bck;
				pred_pos -= len_bck;
				text_ptr -= len_bck;
				i -= len_bck;
			}
#endif

			encode_match(match_pos, len_bck + len_fwd, pred_pos, encoded);

			pred_pos = match_pos + len_bck + len_fwd;
			i += len_bck + len_fwd;
			text_ptr += len_bck + len_fwd;

#ifdef USE_SPARSE_HT
			/*			if (no_prev_literals > thr_literal_len)
							cerr << "litera len: " + to_string(no_prev_literals) + "     \n";*/

			no_prev_literals = 0;
#endif
		}
	}

	for (; i < text_size; ++i)
		encode_literal(text[i], encoded);
}

// *******************************************************************************************
void CLZDiff_V1::Decode(const contig_t& reference, const contig_t& encoded, contig_t& decoded)
{
	uint8_t c;
	uint32_t ref_pos, len;
	uint32_t pred_pos = 0;

	decoded.clear();

	for (auto p = encoded.begin(); p != encoded.end(); )
	{
		if (is_literal(p))
		{
			decode_literal(p, c);
			decoded.emplace_back(c);
			++pred_pos;
		}
		else if (is_Nrun(p))
		{
			decode_Nrun(p, len);
			decoded.insert(decoded.end(), len, N_code);
		}
		else
		{
			decode_match(p, ref_pos, len, pred_pos);
			decoded.insert(decoded.end(), reference.begin() + ref_pos, reference.begin() + ref_pos + len);
			pred_pos = ref_pos + len;
		}
	}
}


// *******************************************************************************************
//
// *******************************************************************************************
void CLZDiff_V2::encode_match(const uint32_t ref_pos, const uint32_t len, const uint32_t pred_pos, contig_t& encoded)
{
	int dif_pos = (int)ref_pos - (int)pred_pos;

	append_int(encoded, dif_pos);
	if (len != ~0u)
	{
		encoded.emplace_back(',');
		append_int(encoded, len - min_match_len);
	}

	encoded.emplace_back('.');
}

// *******************************************************************************************
void CLZDiff_V2::decode_match(contig_t::const_iterator& p, uint32_t& ref_pos, uint32_t& len, uint32_t& pred_pos)
{
	int64_t raw_pos;
	int64_t raw_len;

	read_int(p, raw_pos);
	ref_pos = (uint32_t)(raw_pos + (int64_t)pred_pos);

	if(*p == ',')
	{
		++p;
		read_int(p, raw_len);
		len = (uint32_t)(raw_len + min_match_len);
		++p;		// '.'
	}
	else
	{
		len = ~0u;
		++p;		// '.'
	}
}

// *******************************************************************************************
void CLZDiff_V2::Encode(const contig_t& text, contig_t& encoded)
{
	if (!index_ready)
		prepare_index();

	uint32_t text_size = (uint32_t)text.size();

	encoded.clear();

	if (text_size == reference.size() - key_len)
		if (equal(text.begin(), text.end(), reference.begin()))
			return;														// equal sequences

	MurMur64Hash mmh;

	uint32_t i = 0;
	uint32_t pred_pos = 0;

	const uint8_t* text_ptr = text.data();

#ifdef USE_SPARSE_HT
//	const uint32_t thr_literal_len = 1000;

	uint32_t no_prev_literals = 0;
#endif

	for (; i + key_len < text_size; )
	{
		uint64_t x = get_code(text_ptr);

		if (x == ~0ull)
		{
			uint32_t Nrun_len = get_Nrun_len(text_ptr, text_size - i);

			if (Nrun_len >= min_Nrun_len)
			{
				encode_Nrun(Nrun_len, encoded);
				text_ptr += Nrun_len;
				i += Nrun_len;
#ifdef USE_SPARSE_HT
				/*				if (no_prev_literals > thr_literal_len)
									cerr << "literal len: " + to_string(no_prev_literals) + "    \n";*/
				no_prev_literals = 0;
#endif
			}
			else
			{
				encode_literal(*text_ptr, encoded);

				++i;
				++pred_pos;
				++text_ptr;
#ifdef USE_SPARSE_HT
				++no_prev_literals;
#endif
			}

			continue;
		}

		uint64_t ht_pos = mmh(x) & ht_mask;

		uint32_t len_bck = 0;
		uint32_t len_fwd = 0;
		uint32_t match_pos = 0;
		uint32_t max_len = text_size - i;

		if (short_ht_ver ?
			!find_best_match16((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd) :
			!find_best_match32((uint32_t)ht_pos, text_ptr, max_len, no_prev_literals, match_pos, len_bck, len_fwd))
		{
			encode_literal(*text_ptr, encoded);

			++i;
			++text_ptr;
			++pred_pos;
#ifdef USE_SPARSE_HT
			++no_prev_literals;
#endif
			continue;
		}
		else
		{
#ifdef USE_SPARSE_HT
			if (len_bck)
			{
				for (uint32_t k = 0; k < len_bck; ++k)
					encoded.pop_back();
				match_pos -= len_bck;
				pred_pos -= len_bck;
				text_ptr -= len_bck;
				i -= len_bck;
			}
#endif

			if (match_pos == pred_pos)
			{
				uint32_t e_size = encoded.size();
				for (uint32_t i = 1; i < e_size && i < match_pos; ++i)
				{
					if (encoded[e_size - i] < 'A' || encoded[e_size - i] > 'Z')
						break;
					if (encoded[e_size - i] - 'A' == reference[match_pos - i])
						encoded[e_size - i] = '!';
				}
			}

			if (i + len_bck + len_fwd == text_size && match_pos + len_bck + len_fwd == reference.size() - key_len)			// is match to end of sequence?
				encode_match(match_pos, ~0u, pred_pos, encoded);
			else
				encode_match(match_pos, len_bck + len_fwd, pred_pos, encoded);

			pred_pos = match_pos + len_bck + len_fwd;
			i += len_bck + len_fwd;
			text_ptr += len_bck + len_fwd;

#ifdef USE_SPARSE_HT
			/*			if (no_prev_literals > thr_literal_len)
							cerr << "litera len: " + to_string(no_prev_literals) + "     \n";*/

			no_prev_literals = 0;
#endif
		}
	}

	for (; i < text_size; ++i)
		encode_literal(text[i], encoded);
}

// *******************************************************************************************
void CLZDiff_V2::Decode(const contig_t& reference, const contig_t& encoded, contig_t& decoded)
{
	uint8_t c;
	uint32_t ref_pos, len;
	uint32_t pred_pos = 0;

	decoded.clear();

	for (auto p = encoded.begin(); p != encoded.end(); )
	{
		if (is_literal(p))
		{
			decode_literal(p, c);

			if (c == '!')
				c = reference[pred_pos];
			decoded.emplace_back(c);
			++pred_pos;
		}
		else if (is_Nrun(p))
		{
			decode_Nrun(p, len);
			decoded.insert(decoded.end(), len, N_code);
		}
		else
		{
			decode_match(p, ref_pos, len, pred_pos);
			
			if (len == ~0u)
				len = reference.size() - ref_pos;

			decoded.insert(decoded.end(), reference.begin() + ref_pos, reference.begin() + ref_pos + len);
			pred_pos = ref_pos + len;
		}
	}
}

//*********************************************************************************************

pair<vector<string>, vector<string>> read(string file_name) {
    vector<string> head;
    vector<string> data;
	file_name="/global/homes/v/vidushi/sharedp/projects/compression/agc/decompressed/"+file_name;
    ifstream f(file_name);
    string s = "";
    for (string line; getline(f, line); ) {
        vector<string> temp;
        temp.push_back(line);
        if (line[0] == '>') {
            head.push_back(temp[0].substr(1));
            if (s.length() != 0) {
                data.push_back(s);
                s = "";
            }
        } else {
            s += temp[0];
        }
    }
    data.push_back(s);
    f.close();
    return make_pair(head, data);
}

//*********************************************************************************************

vector<uint8_t> write_in_vector(string data){
    vector<uint8_t> data2; // create a vector to store the data
    char c;
    for(int i=0; i<data.size(); i++) { // read the file character by character
        data2.push_back(static_cast<uint8_t>(dna_code(data[i]))); // store each character in the vector
    }

    return data2;
}

//*********************************************************************************************

void write_to_file(vector<uint8_t> inp, string fn){
    ofstream outFile(fn);
    int temp;
    for (int i = 0; i < inp.size(); ++i) {
        temp = static_cast<int>(inp[i]);
        if (temp == 0){
            outFile.put('A');
        }
        else if (temp == 1){
            outFile.put('C');
        }
        else if (temp == 2){
            outFile.put('G');
        }
        else if (temp == 3){
            outFile.put('T');
        }
    }
    outFile.close();
}

//*******************************************************************************************

pair<uint8_t*,size_t> zstd_compress(vector<uint8_t> NumToCompress){
    uint8_t* com_ptr = NULL;
    size_t NumSize = NumToCompress.size();
    size_t Boundsize = ZSTD_compressBound(NumSize);
    com_ptr =(uint8_t*) malloc(Boundsize);
    size_t ComSize;
    ComSize = ZSTD_compress(com_ptr, Boundsize, NumToCompress.data(), 
    NumToCompress.size()*sizeof(uint8_t), 3);
	//com_ptr.resize(ComSize);
    return make_pair(com_ptr,ComSize);
}

//*******************************************************************************************

vector<uint8_t> zstd_decompress(uint8_t* com_ptr, size_t ComSize){
    uint8_t* decom_ptr = NULL;
    unsigned long long decom_Boundsize;
    decom_Boundsize = ZSTD_getFrameContentSize(com_ptr, ComSize);
    decom_ptr = (uint8_t*)malloc(decom_Boundsize);
    size_t  DecomSize;
    DecomSize = ZSTD_decompress(decom_ptr, decom_Boundsize, com_ptr, ComSize);
    vector<uint8_t> NumAfterDecompress(decom_ptr, decom_ptr+DecomSize/sizeof(uint8_t));
    return NumAfterDecompress;
}

//*********************************************************************************************

unordered_map<string, vector<string>> create_bed(pair<vector<string>, vector<string>> refer){
	unordered_map<string, vector<string>> lifted;
	vector<vector<int>> interval;
    vector<int> temp;
	int count = 0;
    ofstream file("bed.bed");
    for (int i = 0; i < refer.first.size(); i++) {
        string str1 = refer.first[i];
		int a = 0;
        int b = refer.second[i].size() - 1;
        int l = 50000;
		if(b<l){
			temp.push_back(a);
			temp.push_back(b);
			interval.push_back(temp);
			temp.clear();
			file << str1 << "\t" << to_string(a) << "\t" << to_string(b) <<"\n";
			lifted[str1+"_"+to_string(a)+"_"+to_string(b)] = {} ;
			count++;
		}
		else{
			int j=0;
			temp.push_back(0);
			temp.push_back(l);
			interval.push_back(temp);
			temp.clear();
			file << str1 << "\t" << to_string(0) << "\t" << to_string(l) <<"\n";
			lifted[str1+"_"+to_string(0)+"_"+to_string(l)] = {} ;
			count++;
			for(j=1; j < b/l; j++){
				temp.push_back(j*l+1);
				temp.push_back((j+1)*l);
				interval.push_back(temp);
				temp.clear();
				file << str1 << "\t" << to_string(j*l+1) << "\t" << to_string((j+1)*l) <<"\n";
				lifted[str1+"_"+to_string(j*l+1)+"_"+to_string((j+1)*l)] = {} ;
				count++;
			}
			temp.push_back(j*l+1);
			temp.push_back(b);
			interval.push_back(temp);
			temp.clear();
			file << str1 << "\t" << to_string(j*l+1) << "\t" << to_string(b) <<"\n";
			lifted[str1+"_"+to_string(j*l+1)+"_"+to_string(b)] = {} ;
			count++;
    	}
		interval.clear();
	}
        
    file.close();
	return lifted;

}

//*********************************************************************************************

void use_liftover(vector<string> paf, vector<string> out){

	for (int i=0; i<paf.size(); i++){ 
		string s = "/global/homes/v/vidushi/k8 /global/homes/v/vidushi/sharedp/projects/compression/cal/paf.js liftover -l50000 ";
        string s1 = "/global/homes/v/vidushi/sharedp/projects/compression/cal/aligns/";
        system((s+s1+ paf[i] +" bed.bed >"+ out[i]).c_str());
	}

}

//*********************************************************************************************

unordered_map<string, int> mapping_contig(vector<string> recieved){
	unordered_map<string, int> abc;
    for (int i=0; i<recieved.size(); i++) {
		abc[recieved[i]] = i;
    }
	return abc;
}

//*********************************************************************************************

// unordered_map<string, vector<int>> create_lifted_dict(){
	
// }



//*********************************************************************************************

vector<string> split(string s, string d){
    vector<string> ret;
    size_t pos = 0;
    string token;
    while ((pos = s.find(d)) != string::npos) {
        token = s.substr(0, pos);
        ret.push_back(token);
        s.erase(0, pos + d.length());
    }
    ret.push_back(s);
    
    return ret;
}

//************************************************************************************************

unordered_map<string, vector<string>> create_lift(string file_name,unordered_map<string, int> seq1_index,pair<vector<string>, vector<string>> sequence1){
	unordered_map<string, vector<string>> lift1;
	ifstream file(file_name);
	if (file.is_open()) {
		string line;
		while (getline(file, line)) {
			string s = line;
			vector<string> str1 = split(s,"\t");
			vector<string> str2 = split(str1[3],"_");
			lift1[str2[0]+"_"+str2[1]+"_"+str2[2]].clear();
			if(stoi(str1[2])==sequence1.second[seq1_index[str1[0]]].size()){
				lift1[str2[0]+"_"+str2[1]+"_"+str2[2]]={str1[0]+"_"+str1[1]+"_"+to_string(stoi(str1[2])-1)};
			}
			else{
				lift1[str2[0]+"_"+str2[1]+"_"+str2[2]]={str1[0]+"_"+str1[1]+"_"+str1[2]};
			}
		}
		file.close();
	}

	return lift1;
}


//************************************************************************************************

int main() {

    
    vector<string> fa_files={"CHM13Y.fa","HG02148.1.fa","HG02148.2.fa","HG02257.1.fa","HG02257.2.fa","HG02486.1.fa",
    "HG02486.2.fa","HG02559.1.fa","HG02559.2.fa","HG02572.1.fa","HG02572.2.fa","HG02622.1.fa","HG02622.2.fa",
    "HG02630.1.fa","HG02630.2.fa","HG02717.1.fa","HG02717.2.fa","HG02723.1.fa","HG02723.2.fa","HG02818.1.fa",
    "HG02818.2.fa","HG02886.1.fa","HG02886.2.fa","HG03098.1.fa","HG03098.2.fa","HG03453.1.fa","HG03453.2.fa",
    "HG03486.1.fa","HG03486.2.fa","HG03492.1.fa","HG03492.2.fa","HG03516.1.fa","HG03516.2.fa","HG03540.1.fa",
    "HG03540.2.fa","HG03579.1.fa","HG03579.2.fa","NA18906.1.fa","NA18906.2.fa","NA19240.1.fa","NA19240.2.fa",
    "NA20129.1.fa","NA20129.2.fa","NA21309.1.fa","NA21309.2.fa"};

    vector<string> paf_files=  {"halign51.paf","halign52.paf","halign53.paf","halign54.paf","halign55.paf",
                                "halign56.paf","halign57.paf","halign58.paf","halign59.paf","halign60.paf",
                                "halign61.paf","halign62.paf","halign63.paf","halign64.paf","halign65.paf",
                                "halign66.paf","halign67.paf","halign68.paf","halign69.paf","halign70.paf",
                                "halign71.paf","halign72.paf","halign73.paf","halign74.paf","halign75.paf",
                                "halign76.paf","halign77.paf","halign78.paf","halign79.paf","halign80.paf",
                                "halign81.paf","halign82.paf","halign83.paf","halign84.paf","halign85.paf",
                                "halign86.paf","halign87.paf","halign88.paf","halign89.paf","halign90.paf",
                                "halign91.paf","halign92.paf","halign93.paf","halign94.paf"};

	vector<string> output_lift={"lift51.txt","lift52.txt","lift53.txt","lift54.txt","lift55.txt",
                                "lift56.txt","lift57.txt","lift58.txt","lift59.txt","lift60.txt",
                                "lift61.txt","lift62.txt","lift63.txt","lift64.txt","lift65.txt",
                                "lift66.txt","lift67.txt","lift68.txt","lift69.txt","lift70.txt",
                                "lift71.txt","lift72.txt","lift73.txt","lift74.txt","lift75.txt",
                                "lift76.txt","lift77.txt","lift78.txt","lift79.txt","lift80.txt",
                                "lift81.txt","lift82.txt","lift83.txt","lift84.txt","lift85.txt",
                                "lift86.txt","lift87.txt","lift88.txt","lift89.txt","lift90.txt",
                                "lift91.txt","lift92.txt","lift93.txt","lift94.txt"};

    pair<vector<string>, vector<string>> reference = read(fa_files[0]);
    pair<vector<string>, vector<string>> sequence51 = read(fa_files[1]);
    // cout<<"1ran"<<endl;
    pair<vector<string>, vector<string>> sequence52 = read(fa_files[2]);
    // cout<<"2ran"<<endl;
    pair<vector<string>, vector<string>> sequence53 = read(fa_files[3]);
    // cout<<"3ran"<<endl;
    pair<vector<string>, vector<string>> sequence54 = read(fa_files[4]);
    // cout<<"4ran"<<endl;
    pair<vector<string>, vector<string>> sequence55 = read(fa_files[5]);
    // cout<<"5ran"<<endl;
    pair<vector<string>, vector<string>> sequence56 = read(fa_files[6]);
    // cout<<"6ran"<<endl;
    pair<vector<string>, vector<string>> sequence57 = read(fa_files[7]);
    // cout<<"7ran"<<endl;
    pair<vector<string>, vector<string>> sequence58 = read(fa_files[8]);
    // cout<<"8ran"<<endl;
    pair<vector<string>, vector<string>> sequence59 = read(fa_files[9]);
    // cout<<"9ran"<<endl;
    pair<vector<string>, vector<string>> sequence60 = read(fa_files[10]);
    // cout<<"10ran"<<endl;
    pair<vector<string>, vector<string>> sequence61 = read(fa_files[11]);
    // cout<<"11ran"<<endl;
    pair<vector<string>, vector<string>> sequence62 = read(fa_files[12]);
    // cout<<"12ran"<<endl;
    pair<vector<string>, vector<string>> sequence63 = read(fa_files[13]);
    // cout<<"13ran"<<endl;
    pair<vector<string>, vector<string>> sequence64 = read(fa_files[14]);
    // cout<<"14ran"<<endl;
    pair<vector<string>, vector<string>> sequence65 = read(fa_files[15]);
    // cout<<"15ran"<<endl;
    pair<vector<string>, vector<string>> sequence66 = read(fa_files[16]);
    // cout<<"16ran"<<endl;
    pair<vector<string>, vector<string>> sequence67 = read(fa_files[17]);
    // cout<<"17ran"<<endl;
    pair<vector<string>, vector<string>> sequence68 = read(fa_files[18]);
    // cout<<"18ran"<<endl;
    pair<vector<string>, vector<string>> sequence69 = read(fa_files[19]);
    // cout<<"19ran"<<endl;
    pair<vector<string>, vector<string>> sequence70 = read(fa_files[20]);
    // cout<<"20ran"<<endl;
    pair<vector<string>, vector<string>> sequence71 = read(fa_files[21]);
    // cout<<"21ran"<<endl;
    pair<vector<string>, vector<string>> sequence72 = read(fa_files[22]);
    // cout<<"22ran"<<endl;
    pair<vector<string>, vector<string>> sequence73 = read(fa_files[23]);
    // cout<<"23ran"<<endl;
    pair<vector<string>, vector<string>> sequence74 = read(fa_files[24]);
    // cout<<"24ran"<<endl;
    pair<vector<string>, vector<string>> sequence75 = read(fa_files[25]);
    // cout<<"25ran"<<endl;
    pair<vector<string>, vector<string>> sequence76 = read(fa_files[26]);
    // cout<<"26ran"<<endl;
    pair<vector<string>, vector<string>> sequence77 = read(fa_files[27]);
    // cout<<"27ran"<<endl;
    pair<vector<string>, vector<string>> sequence78 = read(fa_files[28]);
    // cout<<"28ran"<<endl;
    pair<vector<string>, vector<string>> sequence79 = read(fa_files[29]);
    // cout<<"29ran"<<endl;
    pair<vector<string>, vector<string>> sequence80 = read(fa_files[30]);
    // cout<<"30ran"<<endl;
    pair<vector<string>, vector<string>> sequence81 = read(fa_files[31]);
    // cout<<"31ran"<<endl;
    pair<vector<string>, vector<string>> sequence82 = read(fa_files[32]);
    // cout<<"32ran"<<endl;
    pair<vector<string>, vector<string>> sequence83 = read(fa_files[33]);
    // cout<<"33ran"<<endl;
    pair<vector<string>, vector<string>> sequence84 = read(fa_files[34]);
    // cout<<"34ran"<<endl;
    pair<vector<string>, vector<string>> sequence85 = read(fa_files[35]);
    // cout<<"35ran"<<endl;
    pair<vector<string>, vector<string>> sequence86 = read(fa_files[36]);
    // cout<<"36ran"<<endl;
    pair<vector<string>, vector<string>> sequence87 = read(fa_files[37]);
    // cout<<"37ran"<<endl;
    pair<vector<string>, vector<string>> sequence88 = read(fa_files[38]);
    // cout<<"38ran"<<endl;
    pair<vector<string>, vector<string>> sequence89 = read(fa_files[39]);
    // cout<<"39ran"<<endl;
    pair<vector<string>, vector<string>> sequence90 = read(fa_files[40]);
    // cout<<"40ran"<<endl;
    pair<vector<string>, vector<string>> sequence91 = read(fa_files[41]);
    // cout<<"41ran"<<endl;
    pair<vector<string>, vector<string>> sequence92 = read(fa_files[42]);
    // cout<<"42ran"<<endl;
    pair<vector<string>, vector<string>> sequence93 = read(fa_files[43]);
    // cout<<"43ran"<<endl;
    pair<vector<string>, vector<string>> sequence94 = read(fa_files[44]);
    // cout<<"44ran"<<endl;
    
	unordered_map<int,vector<string>> list;
    
    unordered_map<string, vector<string>>   lifted = create_bed(reference);
	unordered_map<string, vector<string>>   lift51, lift52, lift53, lift54, lift55, 
                                            lift56, lift57, lift58, lift59, lift60, 
                                            lift61, lift62, lift63, lift64, lift65, 
                                            lift66, lift67, lift68, lift69, lift70, 
                                            lift71, lift72, lift73, lift74, lift75, 
                                            lift76, lift77, lift78, lift79, lift80, 
                                            lift81, lift82, lift83, lift84, lift85, 
											lift86, lift87, lift88, lift89, lift90, 
											lift91, lift92, lift93, lift94;


	use_liftover(paf_files,output_lift);
	//create_lifted_dict();
    unordered_map<string, int> ref_index,   seq51_index, seq52_index, seq53_index, seq54_index, seq55_index, 
                                seq56_index, seq57_index, seq58_index, seq59_index, seq60_index, 
                                seq61_index, seq62_index, seq63_index, seq64_index, seq65_index, 
                                seq66_index, seq67_index, seq68_index, seq69_index, seq70_index, 
                                seq71_index, seq72_index, seq73_index, seq74_index, seq75_index, 
                                seq76_index, seq77_index, seq78_index, seq79_index, seq80_index, 
                                seq81_index, seq82_index, seq83_index, seq84_index, seq85_index, 
								seq86_index, seq87_index, seq88_index, seq89_index, seq90_index, 
								seq91_index, seq92_index, seq93_index, seq94_index;

	ref_index = mapping_contig(reference.first);
	seq51_index = mapping_contig(sequence51.first);
    seq52_index = mapping_contig(sequence52.first);
    seq53_index = mapping_contig(sequence53.first);
    seq54_index = mapping_contig(sequence54.first);
    seq55_index = mapping_contig(sequence55.first);
    seq56_index = mapping_contig(sequence56.first);
    seq57_index = mapping_contig(sequence57.first);
    seq58_index = mapping_contig(sequence58.first);
    seq59_index = mapping_contig(sequence59.first);
    seq60_index = mapping_contig(sequence60.first);
    seq61_index = mapping_contig(sequence61.first);
    seq62_index = mapping_contig(sequence62.first);
    seq63_index = mapping_contig(sequence63.first);
    seq64_index = mapping_contig(sequence64.first);
    seq65_index = mapping_contig(sequence65.first);
    seq66_index = mapping_contig(sequence66.first);
    seq67_index = mapping_contig(sequence67.first);
    seq68_index = mapping_contig(sequence68.first);
    seq69_index = mapping_contig(sequence69.first);
    seq70_index = mapping_contig(sequence70.first);
    seq71_index = mapping_contig(sequence71.first);
    seq72_index = mapping_contig(sequence72.first);
    seq73_index = mapping_contig(sequence73.first);
    seq74_index = mapping_contig(sequence74.first);
    seq75_index = mapping_contig(sequence75.first);
    seq76_index = mapping_contig(sequence76.first);
    seq77_index = mapping_contig(sequence77.first);
    seq78_index = mapping_contig(sequence78.first);
    seq79_index = mapping_contig(sequence79.first);
    seq80_index = mapping_contig(sequence80.first);
    seq81_index = mapping_contig(sequence81.first);
    seq82_index = mapping_contig(sequence82.first);
    seq83_index = mapping_contig(sequence83.first);
    seq84_index = mapping_contig(sequence84.first);
    seq85_index = mapping_contig(sequence85.first);
    seq86_index = mapping_contig(sequence86.first);
    seq87_index = mapping_contig(sequence87.first);
    seq88_index = mapping_contig(sequence88.first);
    seq89_index = mapping_contig(sequence89.first);
    seq90_index = mapping_contig(sequence90.first);
    seq91_index = mapping_contig(sequence91.first);
    seq92_index = mapping_contig(sequence92.first);
    seq93_index = mapping_contig(sequence93.first);
    seq94_index = mapping_contig(sequence94.first);



	// ofstream f("seq1c.bin",ios::out | ios :: binary);
	lift51=create_lift("lift51.txt",seq51_index,sequence51);
    lift52=create_lift("lift52.txt",seq52_index,sequence52);
    lift53=create_lift("lift53.txt",seq53_index,sequence53);
    lift54=create_lift("lift54.txt",seq54_index,sequence54);
    lift55=create_lift("lift55.txt",seq55_index,sequence55);
    lift56=create_lift("lift56.txt",seq56_index,sequence56);
    lift57=create_lift("lift57.txt",seq57_index,sequence57);
    lift58=create_lift("lift58.txt",seq58_index,sequence58);
    lift59=create_lift("lift59.txt",seq59_index,sequence59);
    lift60=create_lift("lift60.txt",seq60_index,sequence60);
    lift61=create_lift("lift61.txt",seq61_index,sequence61);
    lift62=create_lift("lift62.txt",seq62_index,sequence62);
    lift63=create_lift("lift63.txt",seq63_index,sequence63);
    lift64=create_lift("lift64.txt",seq64_index,sequence64);
    lift65=create_lift("lift65.txt",seq65_index,sequence65);
    lift66=create_lift("lift66.txt",seq66_index,sequence66);
    lift67=create_lift("lift67.txt",seq67_index,sequence67);
    lift68=create_lift("lift68.txt",seq68_index,sequence68);
    lift69=create_lift("lift69.txt",seq69_index,sequence69);
    lift70=create_lift("lift70.txt",seq70_index,sequence70);
    lift71=create_lift("lift71.txt",seq71_index,sequence71);
    lift72=create_lift("lift72.txt",seq72_index,sequence72);
    lift73=create_lift("lift73.txt",seq73_index,sequence73);
    lift74=create_lift("lift74.txt",seq74_index,sequence74);
    lift75=create_lift("lift75.txt",seq75_index,sequence75);
    lift76=create_lift("lift76.txt",seq76_index,sequence76);
    lift77=create_lift("lift77.txt",seq77_index,sequence77);
    lift78=create_lift("lift78.txt",seq78_index,sequence78);
    lift79=create_lift("lift79.txt",seq79_index,sequence79);
    lift80=create_lift("lift80.txt",seq80_index,sequence80);
    lift81=create_lift("lift81.txt",seq81_index,sequence81);
    lift82=create_lift("lift82.txt",seq82_index,sequence82);
    lift83=create_lift("lift83.txt",seq83_index,sequence83);
    lift84=create_lift("lift84.txt",seq84_index,sequence84);
    lift85=create_lift("lift85.txt",seq85_index,sequence85);
    lift86=create_lift("lift86.txt",seq86_index,sequence86);
    lift87=create_lift("lift87.txt",seq87_index,sequence87);
    lift88=create_lift("lift88.txt",seq88_index,sequence88);
    lift89=create_lift("lift89.txt",seq89_index,sequence89);
    lift90=create_lift("lift90.txt",seq90_index,sequence90);
    lift91=create_lift("lift91.txt",seq91_index,sequence91);
    lift92=create_lift("lift92.txt",seq92_index,sequence92);
    lift93=create_lift("lift93.txt",seq93_index,sequence93);
    lift94=create_lift("lift94.txt",seq94_index,sequence94);


	for(auto i:lift51){ 
    lifted[i.first].push_back(i.second[0]+"\t"+"51");
    }
    
    for(auto i:lift52){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"52");
        }
        
    for(auto i:lift53){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"53");
        }
        
    for(auto i:lift54){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"54");
        }
        
    for(auto i:lift55){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"55");
        }
        
    for(auto i:lift56){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"56");
        }
        
    for(auto i:lift57){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"57");
        }
        
    for(auto i:lift58){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"58");
        }
        
    for(auto i:lift59){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"59");
        }
        
    for(auto i:lift60){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"60");
        }
        
    for(auto i:lift61){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"61");
        }
        
    for(auto i:lift62){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"62");
        }
        
    for(auto i:lift63){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"63");
        }
        
    for(auto i:lift64){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"64");
        }
        
    for(auto i:lift65){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"65");
        }
        
    for(auto i:lift66){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"66");
        }
        
    for(auto i:lift67){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"67");
        }
        
    for(auto i:lift68){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"68");
        }
        
    for(auto i:lift69){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"69");
        }
        
    for(auto i:lift70){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"70");
        }
        
    for(auto i:lift71){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"71");
        }
        
    for(auto i:lift72){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"72");
        }
        
    for(auto i:lift73){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"73");
        }
        
    for(auto i:lift74){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"74");
        }
        
    for(auto i:lift75){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"75");
        }
        
    for(auto i:lift76){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"76");
        }
        
    for(auto i:lift77){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"77");
        }
        
    for(auto i:lift78){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"78");
        }
        
    for(auto i:lift79){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"79");
        }
        
    for(auto i:lift80){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"80");
        }
        
    for(auto i:lift81){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"81");
        }
        
    for(auto i:lift82){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"82");
        }
        
    for(auto i:lift83){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"83");
        }
        
    for(auto i:lift84){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"84");
        }
        
    for(auto i:lift85){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"85");
        }
        
    for(auto i:lift86){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"86");
        }
        
    for(auto i:lift87){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"87");
        }
        
    for(auto i:lift88){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"88");
        }
        
    for(auto i:lift89){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"89");
        }
        
    for(auto i:lift90){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"90");
        }
        
    for(auto i:lift91){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"91");
        }
        
    for(auto i:lift92){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"92");
        }
        
    for(auto i:lift93){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"93");
        }
        
    for(auto i:lift94){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"94");
        }
	//start compressing
	ofstream f("seq35c.bin",ios::out | ios :: binary);
	for(auto i:lifted){
		if(i.second.size()>0){
            contig_t delta, decom,ref,seqf, compressed, decompressed;
			vector<string> str1 = split(i.first,"_");
            ref = write_in_vector(reference.second[ref_index[str1[0]]].substr(stoi(str1[1]),stoi(str1[2])-stoi(str1[1])));
            for(int j=0; j<i.second.size();j++){
                contig_t seq;
                vector<string> determine = split(i.second[j],"\t");
                vector<string> str2 = split(determine[0],"_");
                if(determine[1]=="51"){
                    seq = write_in_vector(sequence51.second[seq51_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="52"){
                    seq = write_in_vector(sequence52.second[seq52_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="53"){
                    seq = write_in_vector(sequence53.second[seq53_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="54"){
                    seq = write_in_vector(sequence54.second[seq54_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="55"){
                    seq = write_in_vector(sequence55.second[seq55_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="56"){
                    seq = write_in_vector(sequence56.second[seq56_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="57"){
                    seq = write_in_vector(sequence57.second[seq57_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="58"){
                    seq = write_in_vector(sequence58.second[seq58_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="59"){
                    seq = write_in_vector(sequence59.second[seq59_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="60"){
                    seq = write_in_vector(sequence60.second[seq60_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="61"){
                    seq = write_in_vector(sequence61.second[seq61_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="62"){
                    seq = write_in_vector(sequence62.second[seq62_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="63"){
                    seq = write_in_vector(sequence63.second[seq63_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="64"){
                    seq = write_in_vector(sequence64.second[seq64_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="65"){
                    seq = write_in_vector(sequence65.second[seq65_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="66"){
                    seq = write_in_vector(sequence66.second[seq66_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="67"){
                    seq = write_in_vector(sequence67.second[seq67_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="68"){
                    seq = write_in_vector(sequence68.second[seq68_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="69"){
                    seq = write_in_vector(sequence69.second[seq69_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="70"){
                    seq = write_in_vector(sequence70.second[seq70_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="71"){
                    seq = write_in_vector(sequence71.second[seq71_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="72"){
                    seq = write_in_vector(sequence72.second[seq72_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="73"){
                    seq = write_in_vector(sequence73.second[seq73_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="74"){
                    seq = write_in_vector(sequence74.second[seq74_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="75"){
                    seq = write_in_vector(sequence75.second[seq75_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="76"){
                    seq = write_in_vector(sequence76.second[seq76_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="77"){
                    seq = write_in_vector(sequence77.second[seq77_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="78"){
                    seq = write_in_vector(sequence78.second[seq78_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="79"){
                    seq = write_in_vector(sequence79.second[seq79_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="80"){
                    seq = write_in_vector(sequence80.second[seq80_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="81"){
                    seq = write_in_vector(sequence81.second[seq81_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="82"){
                    seq = write_in_vector(sequence82.second[seq82_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="83"){
                    seq = write_in_vector(sequence83.second[seq83_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="84"){
                    seq = write_in_vector(sequence84.second[seq84_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="85"){
                    seq = write_in_vector(sequence85.second[seq85_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="86"){
                    seq = write_in_vector(sequence86.second[seq86_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="87"){
                    seq = write_in_vector(sequence87.second[seq87_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="88"){
                    seq = write_in_vector(sequence88.second[seq88_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="89"){
                    seq = write_in_vector(sequence89.second[seq89_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="90"){
                    seq = write_in_vector(sequence90.second[seq90_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="91"){
                    seq = write_in_vector(sequence91.second[seq91_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="92"){
                    seq = write_in_vector(sequence92.second[seq92_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="93"){
                    seq = write_in_vector(sequence93.second[seq93_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="94"){
                    seq = write_in_vector(sequence94.second[seq94_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }

                seqf.insert(seqf.end(),seq.begin(),seq.end());
            }
            unique_ptr<CLZDiffBase> lz_diff;
			lz_diff = make_unique<CLZDiff_V2>();
			lz_diff->SetMinMatchLen(20);
			lz_diff->Prepare(ref);
			lz_diff->Encode(seqf,delta);
			pair<uint8_t*, size_t> comp = zstd_compress(delta);
			f.write((char*) comp.first, comp.second);
			f << "\t";
        }   

	}
	
	f.close();
	
    return 0;
}



