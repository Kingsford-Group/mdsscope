#ifndef MT_QUEUE_H_
#define MT_QUEUE_H_

#include <cstdint>
#include <atomic>
#include <optional>
#include <mutex>
#include <unistd.h>
#include <vector>
#include <cassert>
#include <ostream>
#include <iostream>

// /* base-2 logarithm, rounding down */
// static inline constexpr uint64_t lg_down(uint64_t x) {
//   return 63U - __builtin_clzl(x);
// }

// /* base-2 logarithm, rounding up */
// static inline constexpr uint64_t lg_up(uint64_t x) {
//   return lg_down(x - 1) + 1;
// }

// A multi-threaded queue for BFS. Two queues really for the current frontier
// and the next frontier. Pop pulls for the current, push appends to the next.
// Pull returns an empty std::optional element when empty. Use swap() to move to
// the next level.
//
// Objects are not deleted from the queue. Not sure that there is a trait that
// says that's OK.
template<typename T>
class mt_queue {
public:
	typedef T value_type;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	typedef value_type* pointer;
	typedef const value_type* const_pointer;
	typedef ssize_t difference_type;
	typedef size_t size_type;

protected:
	std::vector<T> _current;
	size_type _current_size; // Number of slots filled
	std::atomic<size_type> _index; // Location to pop

	std::vector<T> _next;// all pointers in _next.
	std::atomic<size_type> _next_size; // Number of slots filled
									   //

	const size_t chunk;

public:
	mt_queue(size_type size)
		: _current(size, 0)
		, _current_size(0)
		, _index(0)
		, _next(size, 0)
		, _next_size(0)
		, chunk(std::min(getpagesize() / sizeof(T), size))
		{ // std::cout << "chunk " << chunk << ' ' << sizeof(T) << ' ' << size << std::endl; 
		}

	// Is current frontier empty?
	inline bool current_empty() const { return _current_size == 0; }

	// Swap: make the next frontier be the current. Should be call after getting
	// an empty element from pop(), to start eploring the next level
	void swap(bool copy_next = true) {
		// std::cout << "swap " << _current_size << ' ' << _next_size << std::endl;
		std::swap(_current, _next);
		_current_size = _next_size;
		_index = 0;
		_next_size = 0;
	}

	// Pop an element from the current frontier. Returns a std::optional<T>
	auto pop() {
		const auto i = _index++;
		return i < _current_size ? std::optional<value_type>{_current[i]} : std::nullopt;
	}

	std::pair<T*,ssize_t> multi_pop() {
		const auto pos = _index.fetch_add(chunk);
		return std::make_pair(_current.data() + pos, std::min((ssize_t)chunk, (ssize_t)_current_size - (ssize_t)pos));
	}


	// Push an element to the next frontier.
	void push(const T& x) {
		const auto i = _next_size++;
		_next[i] = x;
	}

	std::pair<T*,ssize_t> multi_push() {
		const auto pos = _next_size.fetch_add(chunk);
		return std::make_pair(_next.data() + pos, chunk);
	}

	void clear() {
		_current_size = 0;
		_index = 0;
		_next_size = 0;
	}

	template<typename U>
	friend std::ostream& operator<<(std::ostream&, const mt_queue<U>&);
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const mt_queue<T>& q) {
	os << "current<" << q._index << ',' << q._current_size << ',' << (void*)&q._current[0] << ">[";
	for(size_t i = 0; i < q._current_size; ++i)
		os << q._current[i] << ' ';
	os << "] next<" << q._next_size << ',' << (void*)&q._next[0] << ">[";
	for(size_t i = 0; i < q._next_size; ++i)
			os << q._next[i] << ' ';
	return os << ']';
}

#endif // MT_QUEUE_H_
