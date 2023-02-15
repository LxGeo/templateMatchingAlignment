#include "defs.h"
#include "io_shapefile.h"
#include "io_raster.h"
#include "parameters.h"
#include "template_matcher_alignment.h"
#include "satellites/imd.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace templateMatchingAlignment
	{

		bool polygonsMatcher::pre_check() {

			std::cout << "Pre check input parameters" << std::endl;

			PolygonsShapfileIO target_shape;
			RasterIO ref_raster;
			bool target_loaded = target_shape.load_shapefile(params->template_shapefile, true);
			if (!target_loaded) {
				std::cout << "Error loading shapefile to match at: " << params->template_shapefile << std::endl;
				return false;
			}
			bool ref_loaded = ref_raster.load_raster(params->search_image);
			if (!ref_loaded) {
				std::cout << "Error loading reference raster at: " << params->search_image << std::endl;
				return false;
			}

			if (!target_shape.spatial_refrence->IsSame(ref_raster.spatial_refrence)) {
				std::cout << "Input shapefiles have different spatial reference system!" << std::endl;
				return false;
			}

			IMetaData imd1(params->imd1_path);
			IMetaData imd2(params->imd2_path);

			//output dirs creation
			boost::filesystem::path output_path(params->output_shapefile);
			boost::filesystem::path output_parent_dirname = output_path.parent_path();
			boost::filesystem::path output_temp_path = output_parent_dirname / params->temp_dir;
			params->temp_dir = output_temp_path.string();
			boost::filesystem::create_directory(output_parent_dirname);
			boost::filesystem::create_directory(output_temp_path);

			if (boost::filesystem::exists(output_path) && !params->overwrite_output) {
				std::cout << fmt::format("output shapefile already exists: {}!", output_path.string()) << std::endl;
				std::cout << fmt::format("Add --overwrite_output !", output_path.string()) << std::endl;
				return false;
			}

			return true;

		}

		void polygonsMatcher::run() {
		}

	}
}