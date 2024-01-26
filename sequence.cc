#include "sequence.hpp"

#include <cstring>
#include <sstream>
#include <iostream>

void translated_stream::initialize_table(const char* str, unsigned size) {
	const auto str_len = strlen(str);
	if(str_len == 0) {
		for(unsigned i = 0; i < size; ++i)
			table['0' + i] = i;
	} else if(str_len == size) {
		for(unsigned i = 0; i < strlen(str); ++i) {
			if(!std::isprint(str[i]) || std::isspace(str[i]))
				throw std::runtime_error("Invalid character in alphabet");
			table[str[i]] = i;
		}
	} else {
		throw std::runtime_error("Invalid alphabet size");
	}
}

translated_stream& translated_stream::header() {
	m_seq_name.clear();
	if(is.peek() == '>') {
		is.get();
		std::getline(is, m_seq_name);
	}
	return *this;
}

translated_stream& translated_stream::operator>>(char& c) {
	char rc;
	while(true) {
		if(!(is >> rc)) break;
		if(rc == '>') {
			is.putback(rc);
			c = asize;
			break;
		}
		++offset;
		// std::cout << "read " << offset << ' ' << rc << ' ' << std::isspace(rc) << ' ' << (int)alphabet.table[rc] << '\n';
		if(std::isspace(rc)) continue;
		c = table[rc];
		if(c < 0) [[unlikely]] {
			std::ostringstream msg;
			msg << "Invalid character '" << rc << "' at position " << offset;
			throw std::runtime_error(msg.str());
		}
		break;
	}
	return *this;
}
