#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <vector>
#include <istream>
#include <string>

struct translated_stream {
	std::vector<char> table;
	const unsigned asize;
	std::istream& is;
	size_t offset;
	translated_stream(const char* str, unsigned size, std::istream& i)
		: table(256, -1)
		, asize(size)
		, is(i)
		, offset(0)
		{ initialize_table(str, size); }

	operator bool() const { return (bool)is; }

	translated_stream& operator>>(char& c);

	translated_stream& header(); // Process fasta header if any

	const std::string seq_name() const { return m_seq_name; }

protected:
	void initialize_table(const char* str, unsigned size);
	std::string m_seq_name;
};


#endif // SEQUENCE_H_
