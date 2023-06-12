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
		system(("/home/karan/Downloads/Vidu/cal/k8-0.2.4/k8-Linux paf.js liftover -l50000 "+ paf[i] +" bed.bed >"+ out[i]).c_str());
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

    
    vector<string> fa_files={"CHM13Y.fa","HG002.1.fa","HG002.2.fa","HG00438.1.fa","HG00438.2.fa",
                            "HG005.1.fa","HG005.2.fa","HG00621.1.fa","HG00621.2.fa","HG00673.1.fa",
                            "HG00673.2.fa","HG00733.1.fa","HG00733.2.fa","HG00735.1.fa","HG00735.2.fa",
                            "HG00741.1.fa","HG00741.2.fa","HG01071.1.fa","HG01071.2.fa","HG01106.1.fa",
                            "HG01106.2.fa","HG01109.1.fa","HG01109.2.fa","HG01123.1.fa","HG01123.2.fa",
                            "HG01175.1.fa","HG01175.2.fa","HG01243.1.fa","HG01243.2.fa","HG01258.1.fa",
                            "HG01258.2.fa","HG01358.1.fa","HG01358.2.fa","HG01361.1.fa","HG01361.2.fa",
                            "HG01891.1.fa","HG01891.2.fa","HG01928.1.fa","HG01928.2.fa","HG01952.1.fa",
                            "HG01952.2.fa","HG01978.1.fa","HG01978.2.fa","HG02055.1.fa","HG02055.2.fa",
                            "HG02080.1.fa","HG02080.2.fa","HG02109.1.fa","HG02109.2.fa","HG02145.1.fa","HG02145.2.fa"};

    vector<string> paf_files=  {"halign1.paf","halign2.paf","halign3.paf","halign4.paf","halign5.paf",
                                "halign6.paf","halign7.paf","halign8.paf","halign9.paf","halign10.paf",
                                "halign11.paf","halign12.paf","halign13.paf","halign14.paf","halign15.paf",
                                "halign16.paf","halign17.paf","halign18.paf","halign19.paf","halign20.paf",
                                "halign21.paf","halign22.paf","halign23.paf","halign24.paf","halign25.paf",
                                "halign26.paf","halign27.paf","halign28.paf","halign29.paf","halign30.paf",
                                "halign31.paf","halign32.paf","halign33.paf","halign34.paf","halign35.paf",
                                "halign36.paf","halign37.paf","halign38.paf","halign39.paf","halign40.paf",
                                "halign41.paf","halign42.paf","halign43.paf","halign44.paf","halign45.paf",
                                "halign46.paf","halign47.paf","halign48.paf","halign49.paf","halign50.paf"};

	vector<string> output_lift={"lift1.txt","lift2.txt","lift3.txt","lift4.txt","lift5.txt",
                                "lift6.txt","lift7.txt","lift8.txt","lift9.txt","lift10.txt",
                                "lift11.txt","lift12.txt","lift13.txt","lift14.txt","lift15.txt",
                                "lift16.txt","lift17.txt","lift18.txt","lift19.txt","lift20.txt",
                                "lift21.txt","lift22.txt","lift23.txt","lift24.txt","lift25.txt",
                                "lift26.txt","lift27.txt","lift28.txt","lift29.txt","lift30.txt",
                                "lift31.txt","lift32.txt","lift33.txt","lift34.txt","lift35.txt",
                                "lift36.txt","lift37.txt","lift38.txt","lift39.txt","lift40.txt",
                                "lift41.txt","lift42.txt","lift43.txt","lift44.txt","lift45.txt",
                                "lift46.txt","lift47.txt","lift48.txt","lift49.txt","lift50.txt"};

    pair<vector<string>, vector<string>> reference = read(fa_files[0]);
    pair<vector<string>, vector<string>> sequence1 = read(fa_files[1]);
	pair<vector<string>, vector<string>> sequence2 = read(fa_files[2]);
	pair<vector<string>, vector<string>> sequence3 = read(fa_files[3]);
	pair<vector<string>, vector<string>> sequence4 = read(fa_files[4]);
	pair<vector<string>, vector<string>> sequence5 = read(fa_files[5]);
    pair<vector<string>, vector<string>> sequence6 = read(fa_files[6]);
	pair<vector<string>, vector<string>> sequence7 = read(fa_files[7]);
	pair<vector<string>, vector<string>> sequence8 = read(fa_files[8]);
	pair<vector<string>, vector<string>> sequence9 = read(fa_files[9]);
	pair<vector<string>, vector<string>> sequence10 = read(fa_files[10]);
	pair<vector<string>, vector<string>> sequence11 = read(fa_files[11]);
	pair<vector<string>, vector<string>> sequence12 = read(fa_files[12]);
	pair<vector<string>, vector<string>> sequence13 = read(fa_files[13]);
	pair<vector<string>, vector<string>> sequence14 = read(fa_files[14]);
	pair<vector<string>, vector<string>> sequence15 = read(fa_files[15]);
    pair<vector<string>, vector<string>> sequence16 = read(fa_files[16]);
	pair<vector<string>, vector<string>> sequence17 = read(fa_files[17]);
	pair<vector<string>, vector<string>> sequence18 = read(fa_files[18]);
	pair<vector<string>, vector<string>> sequence19 = read(fa_files[19]);
	pair<vector<string>, vector<string>> sequence20 = read(fa_files[20]);
	pair<vector<string>, vector<string>> sequence21 = read(fa_files[21]);
	pair<vector<string>, vector<string>> sequence22 = read(fa_files[22]);
	pair<vector<string>, vector<string>> sequence23 = read(fa_files[23]);
	pair<vector<string>, vector<string>> sequence24 = read(fa_files[24]);
	pair<vector<string>, vector<string>> sequence25 = read(fa_files[25]);
    pair<vector<string>, vector<string>> sequence26 = read(fa_files[26]);
    pair<vector<string>, vector<string>> sequence27 = read(fa_files[27]);
    pair<vector<string>, vector<string>> sequence28 = read(fa_files[28]);
    pair<vector<string>, vector<string>> sequence29 = read(fa_files[29]);
    pair<vector<string>, vector<string>> sequence30 = read(fa_files[30]);
    pair<vector<string>, vector<string>> sequence31 = read(fa_files[31]);
    pair<vector<string>, vector<string>> sequence32 = read(fa_files[32]);
    pair<vector<string>, vector<string>> sequence33 = read(fa_files[33]);
    pair<vector<string>, vector<string>> sequence34 = read(fa_files[34]);
    pair<vector<string>, vector<string>> sequence35 = read(fa_files[35]);
    pair<vector<string>, vector<string>> sequence36 = read(fa_files[36]);
    pair<vector<string>, vector<string>> sequence37 = read(fa_files[37]);
    pair<vector<string>, vector<string>> sequence38 = read(fa_files[38]);
    pair<vector<string>, vector<string>> sequence39 = read(fa_files[39]);
    pair<vector<string>, vector<string>> sequence40 = read(fa_files[40]);
    pair<vector<string>, vector<string>> sequence41 = read(fa_files[41]);
    pair<vector<string>, vector<string>> sequence42 = read(fa_files[42]);
    pair<vector<string>, vector<string>> sequence43 = read(fa_files[43]);
    pair<vector<string>, vector<string>> sequence44 = read(fa_files[44]);
    pair<vector<string>, vector<string>> sequence45 = read(fa_files[45]);
    pair<vector<string>, vector<string>> sequence46 = read(fa_files[46]);
    pair<vector<string>, vector<string>> sequence47 = read(fa_files[47]);
    pair<vector<string>, vector<string>> sequence48 = read(fa_files[48]);
    pair<vector<string>, vector<string>> sequence49 = read(fa_files[49]);
    pair<vector<string>, vector<string>> sequence50 = read(fa_files[50]);

	// return 0;
    
	unordered_map<int,vector<string>> list;
    
    unordered_map<string, vector<string>>   lifted = create_bed(reference);
	unordered_map<string, vector<string>>   lift1, lift2, lift3, lift4, lift5, 
                                            lift6, lift7, lift8, lift9, lift10, 
                                            lift11, lift12, lift13, lift14, lift15, 
                                            lift16, lift17, lift18, lift19, lift20, 
                                            lift21, lift22, lift23, lift24, lift25, 
                                            lift26, lift27, lift28, lift29, lift30, 
                                            lift31, lift32, lift33, lift34, lift35, 
                                            lift36, lift37, lift38, lift39, lift40, 
                                            lift41, lift42, lift43, lift44, lift45, 
                                            lift46, lift47, lift48, lift49, lift50;


	use_liftover(paf_files,output_lift);
	//create_lifted_dict();
    unordered_map<string, int> ref_index,   seq1_index, seq2_index, seq3_index, seq4_index, seq5_index,
                                            seq6_index , seq7_index, seq8_index, seq9_index, seq10_index,
                                            seq11_index, seq12_index, seq13_index, seq14_index, seq15_index,
                                            seq16_index, seq17_index, seq18_index, seq19_index, seq20_index,
                                            seq21_index, seq22_index, seq23_index, seq24_index, seq25_index,
                                            seq26_index, seq27_index, seq28_index, seq29_index, seq30_index, 
                                            seq31_index, seq32_index, seq33_index, seq34_index, seq35_index, 
                                            seq36_index, seq37_index, seq38_index, seq39_index, seq40_index, 
                                            seq41_index, seq42_index, seq43_index, seq44_index, seq45_index, 
                                            seq46_index, seq47_index, seq48_index, seq49_index, seq50_index;

	ref_index = mapping_contig(reference.first);
	seq1_index = mapping_contig(sequence1.first);
	seq2_index = mapping_contig(sequence2.first);
	seq3_index = mapping_contig(sequence3.first);
	seq4_index = mapping_contig(sequence4.first);
	seq5_index = mapping_contig(sequence5.first);
	seq6_index = mapping_contig(sequence6.first);
	seq7_index = mapping_contig(sequence7.first);
	seq8_index = mapping_contig(sequence8.first);
	seq9_index = mapping_contig(sequence9.first);
	seq10_index = mapping_contig(sequence10.first);
	seq11_index = mapping_contig(sequence11.first);
	seq12_index = mapping_contig(sequence12.first);
	seq13_index = mapping_contig(sequence13.first);
	seq14_index = mapping_contig(sequence14.first);
	seq15_index = mapping_contig(sequence15.first);
	seq16_index = mapping_contig(sequence16.first);
	seq17_index = mapping_contig(sequence17.first);
	seq18_index = mapping_contig(sequence18.first);
	seq19_index = mapping_contig(sequence19.first);
	seq20_index = mapping_contig(sequence20.first);
	seq21_index = mapping_contig(sequence21.first);
	seq22_index = mapping_contig(sequence22.first);
	seq23_index = mapping_contig(sequence23.first);
	seq24_index = mapping_contig(sequence24.first);
	seq25_index = mapping_contig(sequence25.first);
    seq26_index = mapping_contig(sequence26.first);
    seq27_index = mapping_contig(sequence27.first);
    seq28_index = mapping_contig(sequence28.first);
    seq29_index = mapping_contig(sequence29.first);
    seq30_index = mapping_contig(sequence30.first);
    seq31_index = mapping_contig(sequence31.first);
    seq32_index = mapping_contig(sequence32.first);
    seq33_index = mapping_contig(sequence33.first);
    seq34_index = mapping_contig(sequence34.first);
    seq35_index = mapping_contig(sequence35.first);
    seq36_index = mapping_contig(sequence36.first);
    seq37_index = mapping_contig(sequence37.first);
    seq38_index = mapping_contig(sequence38.first);
    seq39_index = mapping_contig(sequence39.first);
    seq40_index = mapping_contig(sequence40.first);
    seq41_index = mapping_contig(sequence41.first);
    seq42_index = mapping_contig(sequence42.first);
    seq43_index = mapping_contig(sequence43.first);
    seq44_index = mapping_contig(sequence44.first);
    seq45_index = mapping_contig(sequence45.first);
    seq46_index = mapping_contig(sequence46.first);
    seq47_index = mapping_contig(sequence47.first);
    seq48_index = mapping_contig(sequence48.first);
    seq49_index = mapping_contig(sequence49.first);
    seq50_index = mapping_contig(sequence50.first);


	// ofstream f("seq1c.bin",ios::out | ios :: binary);
	lift1=create_lift("lift1.txt",seq1_index,sequence1);
	lift2=create_lift("lift2.txt",seq2_index,sequence2);
	lift3=create_lift("lift3.txt",seq3_index,sequence3);
	lift4=create_lift("lift4.txt",seq4_index,sequence4);
	lift5=create_lift("lift5.txt",seq5_index,sequence5);
	lift6=create_lift("lift6.txt",seq6_index,sequence6);
	lift7=create_lift("lift7.txt",seq7_index,sequence7);
	lift8=create_lift("lift8.txt",seq8_index,sequence8);
	lift9=create_lift("lift9.txt",seq9_index,sequence9);
	lift10=create_lift("lift10.txt",seq10_index,sequence10);
	lift11=create_lift("lift11.txt",seq11_index,sequence11);
	lift12=create_lift("lift12.txt",seq12_index,sequence12);
	lift13=create_lift("lift13.txt",seq13_index,sequence13);
	lift14=create_lift("lift14.txt",seq14_index,sequence14);
	lift15=create_lift("lift15.txt",seq15_index,sequence15);
	lift16=create_lift("lift16.txt",seq16_index,sequence16);
	lift17=create_lift("lift17.txt",seq17_index,sequence17);
	lift18=create_lift("lift18.txt",seq18_index,sequence18);
	lift19=create_lift("lift19.txt",seq19_index,sequence19);
	lift20=create_lift("lift20.txt",seq20_index,sequence20);
	lift21=create_lift("lift21.txt",seq21_index,sequence21);
	lift22=create_lift("lift22.txt",seq22_index,sequence22);
	lift23=create_lift("lift23.txt",seq23_index,sequence23);
	lift24=create_lift("lift24.txt",seq24_index,sequence24);
	lift25=create_lift("lift25.txt",seq25_index,sequence25);
    lift26=create_lift("lift26.txt",seq26_index,sequence26);
    lift27=create_lift("lift27.txt",seq27_index,sequence27);
    lift28=create_lift("lift28.txt",seq28_index,sequence28);
    lift29=create_lift("lift29.txt",seq29_index,sequence29);
    lift30=create_lift("lift30.txt",seq30_index,sequence30);
    lift31=create_lift("lift31.txt",seq31_index,sequence31);
    lift32=create_lift("lift32.txt",seq32_index,sequence32);
    lift33=create_lift("lift33.txt",seq33_index,sequence33);
    lift34=create_lift("lift34.txt",seq34_index,sequence34);
    lift35=create_lift("lift35.txt",seq35_index,sequence35);
    lift36=create_lift("lift36.txt",seq36_index,sequence36);
    lift37=create_lift("lift37.txt",seq37_index,sequence37);
    lift38=create_lift("lift38.txt",seq38_index,sequence38);
    lift39=create_lift("lift39.txt",seq39_index,sequence39);
    lift40=create_lift("lift40.txt",seq40_index,sequence40);
    lift41=create_lift("lift41.txt",seq41_index,sequence41);
    lift42=create_lift("lift42.txt",seq42_index,sequence42);
    lift43=create_lift("lift43.txt",seq43_index,sequence43);
    lift44=create_lift("lift44.txt",seq44_index,sequence44);
    lift45=create_lift("lift45.txt",seq45_index,sequence45);
    lift46=create_lift("lift46.txt",seq46_index,sequence46);
    lift47=create_lift("lift47.txt",seq47_index,sequence47);
    lift48=create_lift("lift48.txt",seq48_index,sequence48);
    lift49=create_lift("lift49.txt",seq49_index,sequence49);
    lift50=create_lift("lift50.txt",seq50_index,sequence50);


	for(auto i:lift1){
		lifted[i.first].push_back(i.second[0]+"\t"+"1");
	}
	for(auto i:lift2){
		lifted[i.first].push_back(i.second[0]+"\t"+"2");
	}
	for(auto i:lift3){
		lifted[i.first].push_back(i.second[0]+"\t"+"3");
	}
	for(auto i:lift4){
		lifted[i.first].push_back(i.second[0]+"\t"+"4");
	}
	for(auto i:lift5){
		lifted[i.first].push_back(i.second[0]+"\t"+"5");
	}
	for(auto i:lift6){
		lifted[i.first].push_back(i.second[0]+"\t"+"6");
	}
	for(auto i:lift7){
		lifted[i.first].push_back(i.second[0]+"\t"+"7");
	}
	for(auto i:lift8){
		lifted[i.first].push_back(i.second[0]+"\t"+"8");
	}
	for(auto i:lift9){
		lifted[i.first].push_back(i.second[0]+"\t"+"9");
	}
	for(auto i:lift10){
		lifted[i.first].push_back(i.second[0]+"\t"+"10");
	}
	for(auto i:lift11){
		lifted[i.first].push_back(i.second[0]+"\t"+"11");
	}
	for(auto i:lift12){
		lifted[i.first].push_back(i.second[0]+"\t"+"12");
	}
	for(auto i:lift13){
		lifted[i.first].push_back(i.second[0]+"\t"+"13");
	}
	for(auto i:lift14){
		lifted[i.first].push_back(i.second[0]+"\t"+"14");
	}
	for(auto i:lift15){
		lifted[i.first].push_back(i.second[0]+"\t"+"15");
	}
	for(auto i:lift16){
		lifted[i.first].push_back(i.second[0]+"\t"+"16");
	}
	for(auto i:lift17){
		lifted[i.first].push_back(i.second[0]+"\t"+"17");
	}
	for(auto i:lift18){
		lifted[i.first].push_back(i.second[0]+"\t"+"18");
	}
	for(auto i:lift19){
		lifted[i.first].push_back(i.second[0]+"\t"+"19");
	}
	for(auto i:lift20){
		lifted[i.first].push_back(i.second[0]+"\t"+"20");
	}
	for(auto i:lift21){
		lifted[i.first].push_back(i.second[0]+"\t"+"21");
	}
	for(auto i:lift22){
		lifted[i.first].push_back(i.second[0]+"\t"+"22");
	}
	for(auto i:lift23){
		lifted[i.first].push_back(i.second[0]+"\t"+"23");
	}
	for(auto i:lift24){
		lifted[i.first].push_back(i.second[0]+"\t"+"24");
	}
	for(auto i:lift25){
		lifted[i.first].push_back(i.second[0]+"\t"+"25");
	}
    for(auto i:lift26){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"26");
    }
    
    for(auto i:lift27){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"27");
        }
        
    for(auto i:lift28){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"28");
        }
        
    for(auto i:lift29){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"29");
        }
        
    for(auto i:lift30){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"30");
        }
        
    for(auto i:lift31){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"31");
        }
        
    for(auto i:lift32){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"32");
        }
        
    for(auto i:lift33){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"33");
        }
        
    for(auto i:lift34){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"34");
        }
        
    for(auto i:lift35){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"35");
        }
        
    for(auto i:lift36){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"36");
        }
        
    for(auto i:lift37){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"37");
        }
        
    for(auto i:lift38){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"38");
        }
        
    for(auto i:lift39){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"39");
        }
        
    for(auto i:lift40){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"40");
        }
        
    for(auto i:lift41){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"41");
        }
        
    for(auto i:lift42){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"42");
        }
        
    for(auto i:lift43){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"43");
        }
        
    for(auto i:lift44){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"44");
        }
        
    for(auto i:lift45){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"45");
        }
        
    for(auto i:lift46){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"46");
        }
        
    for(auto i:lift47){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"47");
        }
        
    for(auto i:lift48){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"48");
        }
        
    for(auto i:lift49){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"49");
        }
        
    for(auto i:lift50){ 
        lifted[i.first].push_back(i.second[0]+"\t"+"50");
        }
	//start compressing
	ofstream f("seq50c.bin",ios::out | ios :: binary);
	for(auto i:lifted){
		if(i.second.size()>0){
            contig_t delta, decom,ref,seqf, compressed, decompressed;
			vector<string> str1 = split(i.first,"_");
            ref = write_in_vector(reference.second[ref_index[str1[0]]].substr(stoi(str1[1]),stoi(str1[2])-stoi(str1[1])));
            for(int j=0; j<i.second.size();j++){
                contig_t seq;
                vector<string> determine = split(i.second[j],"\t");
                vector<string> str2 = split(determine[0],"_");
                if(determine[1]=="1"){
                    seq = write_in_vector(sequence1.second[seq1_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                else if(determine[1]=="2"){
                    seq = write_in_vector(sequence2.second[seq2_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                else if(determine[1]=="3"){
                    seq = write_in_vector(sequence3.second[seq3_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                else if(determine[1]=="4"){
                    seq = write_in_vector(sequence4.second[seq4_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                else if(determine[1]=="5"){
                    seq = write_in_vector(sequence5.second[seq5_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="6"){
                    seq = write_in_vector(sequence6.second[seq6_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="7"){
                    seq = write_in_vector(sequence7.second[seq7_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="8"){
                    seq = write_in_vector(sequence8.second[seq8_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="9"){
                    seq = write_in_vector(sequence9.second[seq9_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="10"){
                    seq = write_in_vector(sequence10.second[seq10_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="11"){
                    seq = write_in_vector(sequence11.second[seq11_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                else if(determine[1]=="12"){
                    seq = write_in_vector(sequence12.second[seq12_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                else if(determine[1]=="13"){
                    seq = write_in_vector(sequence13.second[seq13_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                else if(determine[1]=="14"){
                    seq = write_in_vector(sequence14.second[seq14_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                else if(determine[1]=="15"){
                    seq = write_in_vector(sequence15.second[seq15_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="16"){
                    seq = write_in_vector(sequence16.second[seq16_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="17"){
                    seq = write_in_vector(sequence17.second[seq17_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="18"){
                    seq = write_in_vector(sequence18.second[seq18_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="19"){
                    seq = write_in_vector(sequence19.second[seq19_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="20"){
                    seq = write_in_vector(sequence20.second[seq20_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="21"){
                    seq = write_in_vector(sequence21.second[seq21_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="22"){
                    seq = write_in_vector(sequence22.second[seq22_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="23"){
                    seq = write_in_vector(sequence23.second[seq23_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="24"){
                    seq = write_in_vector(sequence24.second[seq24_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
				else if(determine[1]=="25"){
                    seq = write_in_vector(sequence25.second[seq25_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                else if(determine[1]=="26"){
                    seq = write_in_vector(sequence26.second[seq26_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="27"){
                    seq = write_in_vector(sequence27.second[seq27_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="28"){
                    seq = write_in_vector(sequence28.second[seq28_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="29"){
                    seq = write_in_vector(sequence29.second[seq29_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="30"){
                    seq = write_in_vector(sequence30.second[seq30_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="31"){
                    seq = write_in_vector(sequence31.second[seq31_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="32"){
                    seq = write_in_vector(sequence32.second[seq32_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="33"){
                    seq = write_in_vector(sequence33.second[seq33_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="34"){
                    seq = write_in_vector(sequence34.second[seq34_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="35"){
                    seq = write_in_vector(sequence35.second[seq35_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="36"){
                    seq = write_in_vector(sequence36.second[seq36_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="37"){
                    seq = write_in_vector(sequence37.second[seq37_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="38"){
                    seq = write_in_vector(sequence38.second[seq38_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="39"){
                    seq = write_in_vector(sequence39.second[seq39_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="40"){
                    seq = write_in_vector(sequence40.second[seq40_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="41"){
                    seq = write_in_vector(sequence41.second[seq41_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="42"){
                    seq = write_in_vector(sequence42.second[seq42_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="43"){
                    seq = write_in_vector(sequence43.second[seq43_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="44"){
                    seq = write_in_vector(sequence44.second[seq44_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="45"){
                    seq = write_in_vector(sequence45.second[seq45_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="46"){
                    seq = write_in_vector(sequence46.second[seq46_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="47"){
                    seq = write_in_vector(sequence47.second[seq47_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="48"){
                    seq = write_in_vector(sequence48.second[seq48_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="49"){
                    seq = write_in_vector(sequence49.second[seq49_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
                }
                    
                else if(determine[1]=="50"){
                    seq = write_in_vector(sequence50.second[seq50_index[str2[0]]].substr(stoi(str2[1]),stoi(str2[2])-stoi(str2[1])));
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



