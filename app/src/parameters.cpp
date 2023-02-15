#include "parameters.h"

namespace LxGeo
{
	namespace templateMatchingAlignment
	{
		Parameters::Parameters(int argc, char* argv[])
		{
			init();
			CLI::App app{ "Template Matching Alignment" };
			app.add_option("--tshp", template_shapefile, "Polygons shapefile of rooftops to align!")->required()->check(CLI::ExistingFile);
			app.add_option("--search_image", search_image, "Reference image used as search space for polygon patches!")->required()->check(CLI::ExistingFile);
			app.add_option("--imd1", imd1_path, "Metadata file respective to template image (template shapefile)!")->required()->check(CLI::ExistingFile);
			app.add_option("--imd2", imd2_path, "Metadata file respective to search image!")->required()->check(CLI::ExistingFile);
			app.add_option("-o, --output", output_shapefile, "Output path of matched polygons shapefile!");
			try {
				\
					(app).parse((argc), (argv));                                                                                   \
			}
			catch (const CLI::ParseError& e) {
				\
					(app).exit(e);                                                                                          \
			}
			//app.parse(argc, argv);
		}


		Parameters::~Parameters()
		{
		}


		bool Parameters::initialized()
		{
			return !template_shapefile.empty();
		}


		void Parameters::init()
		{
			printed_help = false;

			template_shapefile.clear();
			search_image.clear();
			imd1_path.clear(); imd2_path.clear();

			output_shapefile = "result.shp";
			temp_dir = "temp_dir";
			overwrite_output = false;

		}

		
		Parameters* params = nullptr;
	}
}