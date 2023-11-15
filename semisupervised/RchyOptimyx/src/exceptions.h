#ifndef _EXCEPTIONS_
#define _EXCEPTIONS_

#include <string>
using namespace std;

namespace tfl
{

	class Exception
	{
		private:
			string data;

		public:
			Exception() { }
			Exception(const string &desc) : data(desc) { }
			Exception(const char *desc) : data(desc) { }

			string getMsg() { return data; }
	};

	class NumberGeneratorOverflowException : public Exception
	{
	};
}
#endif
