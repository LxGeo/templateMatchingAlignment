#pragma once
#include "defs.h"
#include "cli/base_parameters.h"
#include "CLI/CLI.hpp"


namespace LxGeo
{
	namespace templateMatchingAlignment
	{
		class Parameters: public baseParameters
		{
		public:
			Parameters(int argc, char* argv[]);

			~Parameters();

			bool initialized();

		protected:
			void init();

		public:
			bool printed_help;

			std::string template_shapefile;
			std::string template_image;
			std::string search_image;
			std::string imd1_path;
			std::string imd2_path;
			std::string output_shapefile;
			std::string temp_dir;

			bool overwrite_output;

		};

		extern Parameters* params;
	}
}
