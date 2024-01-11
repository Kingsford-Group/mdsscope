#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <vector>
#include <istream>

struct translated_stream {
	std::vector<char> table;
	std::istream& is;
	size_t offset;
	translated_stream(const char* str, unsigned size, std::istream& i)
		: table(256, -1)
		, is(i)
		, offset(0)
		{ initialize_table(str, size); }

	operator bool() const { return (bool)is; }

	translated_stream& operator>>(char& c);


protected:
	void initialize_table(const char* str, unsigned size);
};


#endif // SEQUENCE_H_
