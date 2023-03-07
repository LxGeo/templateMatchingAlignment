#include "defs.h"
#include "io_shapefile.h"
#include "io_raster.h"
#include "parameters.h"
#include "template_matcher_alignment.h"
#include "satellites/imd.h"
#include "geometry_constructor/grid.h"
#include "geometry_rasterizer/individual_masks.h"
#include "soaked_geometries/def_soaked_geometries.h"
#include "tm_manager.h"
#include <fstream>
#include <nlohmann/json.hpp>
#include "satellites/formulas.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace templateMatchingAlignment
	{

		bool polygonsMatcher::pre_check() {

			std::cout << "Pre check input parameters" << std::endl;

			PolygonsShapfileIO target_shape;
			RasterIO template_raster, search_raster;
			bool target_loaded = target_shape.load_shapefile(params->template_shapefile, true);
			if (!target_loaded) {
				std::cout << "Error loading shapefile to match at: " << params->template_shapefile << std::endl;
				return false;
			}
			bool s_image_loaded = search_raster.load_raster(params->search_image);
			if (!s_image_loaded) {
				std::cout << "Error loading search image at: " << params->search_image << std::endl;
				return false;
			}
			bool t_image_loaded = template_raster.load_raster(params->template_image);
			if (!t_image_loaded) {
				std::cout << "Error loading search image at: " << params->template_image << std::endl;
				return false;
			}

			if (!target_shape.spatial_refrence->IsSame(search_raster.spatial_refrence) || !target_shape.spatial_refrence->IsSame(template_raster.spatial_refrence)) {
				std::cout << "Ensure that all input images and shapefiles have the same spatial reference system!" << std::endl;
				return false;
			}
			
			search_radius_pixels = params->search_radius_pixels;

			if (!params->couple_path.empty()) {
				std::ifstream f(params->couple_path);
				nlohmann::json data = nlohmann::json::parse(f);
				couple_rotation_angle = data["rotation_angle"];
				couple_v_displacement = data["v_disp"];
			}
			else {
				IMetaData imd1(params->imd1_path);
				IMetaData imd2(params->imd2_path);
				couple_rotation_angle = -DEGS(compute_rotation_angle(RADS(imd1.satAzimuth), RADS(imd1.satElevation), RADS(imd2.satAzimuth), RADS(imd2.satElevation)));
			}

			

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

			PolygonsShapfileIO target_shape;
			target_shape.load_shapefile(params->template_shapefile, true);			

			std::vector<OGREnvelope> grid_space = create_rectangular_grid(
				target_shape.bounding_box,
				params->optional_numerical_parameters["xstep"],
				params->optional_numerical_parameters["ystep"],
				params->optional_numerical_parameters["xsize"],
				params->optional_numerical_parameters["ysize"]
				);
			std::string grid_temp_path = params->temp_dir + "/grid.shp";
			auto grid_shp = PolygonsShapfileIO(grid_temp_path, target_shape.spatial_refrence);			
			grid_shp.write_shapefile(grid_to_geoms_with_attributes(grid_space));

			std::vector<Polygon_with_attributes> out;
			if (params->optional_str_parameters["device"] == "0")
				out = estimate_shift_one_block_1d<cv::Mat>(grid_space[0], params->tm_method, -couple_rotation_angle);
			else if (params->optional_str_parameters["device"] == "1")
				out = estimate_shift_one_block_1d<cv::cuda::GpuMat>(grid_space[0], params->tm_method, -couple_rotation_angle);

			PolygonsShapfileIO output_shapefile(params->output_shapefile, target_shape.spatial_refrence);
			output_shapefile.write_shapefile(out);
			std::cout << "end!" << std::endl;

		}

		

	}
}