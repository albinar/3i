#ifndef _NUMBER_GENERATOR_
#define _NUMBER_GENERATOR_

namespace tfl
{

	class NumberGenerator
	{
		private:
			char *_data;
			char *_text;
			int len, base, k;
			bool started;

			int firstIncreasableIndex();
			int firstZeroGreaterThanANonZeroIndex();
			int greatestNonZeroIndex();
			//reset all nonzero less significant digits than "index" to 1
			void resetLessSignificantNonZeroDigits(int index);
			void dataToText(const char *data, char *text);

			NumberGenerator() { }

		public:
			NumberGenerator(const NumberGenerator &other);
			NumberGenerator &operator=(const NumberGenerator &other);

			NumberGenerator(int len, int base, int k);
			NumberGenerator(int len, int base, const char *text);

			NumberGenerator &first();
			NumberGenerator &next();
			bool hasNext();
			int nonZeroCount();

			int firstDiff(const NumberGenerator &other);

			const char *text();
			const char *data();
			//zero based neighbor number
			NumberGenerator neighbor(int n_k);
			
			~NumberGenerator();
	};

}


#endif
