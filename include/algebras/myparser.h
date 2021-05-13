/* 
** This header contains functions to convert between containers and strings or io streams
*/
#pragma once
#include <iostream>

/* Consume `pattern` from `sin`. Set badbit if not found. */
void consume(std::istream& sin, const char* pattern);
inline std::istream& operator>>(std::istream& sin, const char* pattern) { consume(sin, pattern); return sin; }

/* Load container from an istream */
template <typename Container>
void load_vector(std::istream& sin, Container& cont,
	const char* left, const char* sep, const char* right)
{
	cont.clear();
	sin >> std::ws; consume(sin, left);
	typename Container::value_type ele;
	while (sin.good()) {
		sin >> ele;
		cont.push_back(ele);
		sin >> std::ws; consume(sin, sep);
	}
	if (sin.bad()) {
		sin.clear();
		consume(sin, right);
	}
}

/* Load the container from an istream */
template <typename Container, typename _Fn_load>
void load_vector(std::istream& sin, Container& cont,
	const char* left, const char* sep, const char* right,
	_Fn_load load)
{
	cont.clear();
	sin >> std::ws; consume(sin, left);
	typename Container::value_type ele;
	while (sin.good()) {
		load(sin, ele);
		cont.push_back(std::move(ele));
		sin >> std::ws; consume(sin, sep);
	}
	if (sin.bad()) {
		sin.clear();
		consume(sin, right);
	}
}

/* Serialize the container to an ostream */
template <typename Container>
void dump_vector(std::ostream& sout, const Container& cont,
	const char* left, const char* sep, const char* right)
{
	sout << left;
	for (auto i = cont.begin(); i != cont.end(); ++i){
		if (i != cont.begin())
			sout << sep;
		sout << *i;
	}
	sout << right;
}

/* Serialize the container to an ostream */
template <typename Container, typename _Fn_dump>
void dump_vector(std::ostream& sout, const Container& cont,
	const char* left, const char* sep, const char* right,
	_Fn_dump dump)
{
	sout << left;
	for (auto i = cont.begin(); i != cont.end(); ++i) {
		if (i != cont.begin())
			sout << sep;
		dump(sout, *i);
	}
	sout << right;
}
