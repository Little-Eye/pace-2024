////////////////////////////////
/// usage : 1.	a fix-size loop queue.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_LOOP_QUEUE_H
#define CN_HUST_GOAL_COMMON_LOOP_QUEUE_H


#include <vector>


namespace goal {

template<typename T, typename IndexType = int>
class LoopQueue {
public:
	LoopQueue() {}
	LoopQueue(IndexType size) : len(size), head(0), tail(0), q(size) {}

	void clear() { head = tail = 0; }


	T& front() { return q[head]; }
	T& back() { return q[prevIndex(tail)]; }
	const T& front() const { return q[head]; }
	const T& back() const { return q[prevIndex(tail)]; }

	T& pushBack() { IndexType t = tail; increaseIndex(tail); return q[t]; }
	T& pushFront() { return q[decreaseIndex(head)]; }
	T& pushBack(const T& item) { return pushBack() = item; }
	T& pushFront(const T& item) { return pushFront() = item; }

	void popBack() { decreaseIndex(tail); }
	void popFront() { increaseIndex(head); }

	bool empty() const { return (head == tail); }
	bool full() const { return head == nextIndex(tail); }

protected:
	IndexType& increaseIndex(IndexType& index) const {
		if ((++index) >= len) { index = 0; }
		return index;
	}
	IndexType& decreaseIndex(IndexType& index) const {
		if ((--index) < 0) { index += len; }
		return index;
	}
	IndexType prevIndex(IndexType index) const { return decreaseIndex(index); }
	IndexType nextIndex(IndexType index) const { return increaseIndex(index); }


	IndexType len;
	IndexType head;
	IndexType tail;
	std::vector<T> q;
};

template<typename T, typename IndexType = int>
class SlidingWindow {
public:
	SlidingWindow(IndexType size, IndexType initIndex = 0) : len(size), cur(initIndex), q(size) {}


	T& now() { return q[cur]; }
	const T& now() const { return q[cur]; }

	T& next(IndexType offset) { return q[shiftIndex(cur, offset)]; }

	void slide1() {
		q[cur].clear();
		increaseIndex(cur);
	}

protected:
	IndexType shiftIndex(IndexType offset) { return (cur += offset) %= len; }

	IndexType increaseIndex() {
		if ((++cur) >= len) { cur -= len; }
		return cur;
	}


	IndexType len;
	IndexType cur;
	std::vector<T> q;
};


template<typename T, typename IndexType = int>
class RoundRobin {
public:
	RoundRobin(IndexType capacity = 64) : head(0), cur(0), tail(0) { q.reserve(capacity); }


	IndexType cursor() const { return cur; }
	void reset() { cur = 0; }

	T& now() { return q[cur]; }
	const T& now() const { return q[cur]; }
	T& operator[](IndexType index) { return q[index]; }
	const T& operator[](IndexType index) const { return q[index]; }

	T& next() { return q[increaseIndex()]; }

	void skip() { ++head; }
	void append(const T& obj) { ++tail; q.emplace_back(obj); }
	void resize(IndexType size) {
		tail = size;
		q.resize(size);
	}
	void resize(IndexType size, const T& obj) {
		tail = size;
		q.resize(size, obj);
	}

	IndexType span() const { return tail - head; }
	IndexType size() const { return tail; }
	bool empty() const { return head >= tail; }

protected:
	IndexType increaseIndex() {
		if ((++cur) >= tail) { cur = head; }
		return cur;
	}


	IndexType head;
	IndexType cur;
	IndexType tail;
	std::vector<T> q;
};

}


#endif // CN_HUST_GOAL_COMMON_LOOP_QUEUE_H
