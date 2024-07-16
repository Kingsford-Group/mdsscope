#ifndef SIMPLE_THREAD_POOL_H_
#define SIMPLE_THREAD_POOL_H_

#include <thread>
#include <barrier>
#include <functional>
#include <vector>
#include <iostream>

template<typename Fn>
class simple_thread_pool {
protected:
	std::barrier<std::function<void(void)>> _wait_barrier, _done_barrier;
	std::vector<std::thread> _ths;
	volatile bool _done;
	Fn _work;

	void run_thread(unsigned index) {
		while(true) {
			_wait_barrier.arrive_and_wait();
			if(_done) break;
			try {
				_work(index);
			} catch(...) {
				std::cerr << "Pool work ended with an exception!";
			}
			_done_barrier.arrive_and_wait();
		}
	}

public:
	simple_thread_pool(unsigned nb_threads)
		: _wait_barrier(nb_threads + 1, [](){}) // Barriers, no completion function
		, _done_barrier(nb_threads + 1, [](){})
		, _done(false)
		{
			for(unsigned i = 0; i < nb_threads; ++i)
				_ths.push_back(std::thread(&simple_thread_pool::run_thread, std::ref(*this), i));
		}

	void set_work(Fn fn) { _work = fn; }

	// Must be called when no thread is currently working. I.e., before calling
	// start(), or after start() has completed. Wait for the threads to
	// terminate.
	void stop() {
		_done = true;
		_wait_barrier.arrive_and_wait();
		for(auto& th : _ths)
			th.join();
	}

	// Start working by all the threads in the pool (i.e., execute the _work
	// function). Returns when the work is done (all the threads finished
	// executing _work). Concurrent call to start() are not supported.
	void start() {
		_wait_barrier.arrive_and_wait(); // Let all the thread start working
		_done_barrier.arrive_and_wait(); // Wait for the threads to be done working
	}
};


#endif // SIMPLE_THREAD_POOL_H_
