#include "parameters.h"
#include "cli/common_options/grid_options.h"
#include "cli/common_options/device.h"
#include "cli/common_options/single_value_enum.h"
#include "template_matcher_alignment.h"

namespace LxGeo
{
	namespace templateMatchingAlignment
	{
		Parameters::Parameters(int argc, char* argv[])
		{
			init();
			CLI::App app{ "Template Matching Alignment" };
			app.add_option("--tshp", template_shapefile, "Polygons shapefile of rooftops to align!")->required()->check(CLI::ExistingFile);
			app.add_option("--template_image", template_image, "Template image used to soak polygon patches!")->required()->check(CLI::ExistingFile);
			app.add_option("--search_image", search_image, "Reference image used as search space for polygon patches!")->required()->check(CLI::ExistingFile);
			app.add_option("--couple_path", couple_path, "Json file containing epipolar creation parameters!")->check(CLI::ExistingFile);
			app.add_option("--imd1", imd1_path, "Metadata file respective to template image (template shapefile)!")->check(CLI::ExistingFile);
			app.add_option("--imd2", imd2_path, "Metadata file respective to search image!")->check(CLI::ExistingFile);
			app.add_option("-o, --output", output_shapefile, "Output path of matched polygons shapefile!");
			add_grid_options(app, *this);
			add_device_option(app, *this);
			// Adding template method option
			std::map<std::string, TemplateMatchingMethod> tm_method_map{
				{"sq_diff", TemplateMatchingMethod::sq_diff},
				{"sq_diff_N", TemplateMatchingMethod::sq_diff_N},
				{"ccorr", TemplateMatchingMethod::ccorr},
				{"ccorr_N", TemplateMatchingMethod::ccorr_N},
				{"ccoef", TemplateMatchingMethod::ccoef},
				{"ccoef_N", TemplateMatchingMethod::ccoef_N},
				{"ssim", TemplateMatchingMethod::ssim},
				{"mu_inf_N", TemplateMatchingMethod::mu_inf_N}
			};
			add_single_value_enum_option(app, this->tm_method, tm_method_map, "-m,--tm_method", "Method used for template matching");
			app.add_option("--search_radius", search_radius_pixels, "Search radius in pixels!")->check(CLI::Range(1, 10000));

			try {
				\
					(app).parse((argc), (argv));\
			}
			catch (const CLI::ParseError& e) {
				\
					(app).exit(e);\
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
			imd1_path.clear(); imd2_path.clear(); couple_path.clear();
			output_shapefile = "result.shp";
			temp_dir = "temp_dir";
			overwrite_output = false;
			TemplateMatchingMethod tm_method = TemplateMatchingMethod::sq_diff;
			int search_radius_pixels = 300;
		}

		
		Parameters* params = nullptr;
	}
}