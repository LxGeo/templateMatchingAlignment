#pragma once
#include "defs.h"
#include "geometries_with_attributes/polygon_with_attributes.h"
#include "tqdm/tqdm.h"

namespace LxGeo
{
	namespace templateMatchingAlignment
	{
		/**
		*  A vectorTransformerApp class to manage running required steps to generate final transformed vector.
		*/
		class polygonsMatcher
		{

		public:
			polygonsMatcher() {};

			~polygonsMatcher() {};

			/**
			*  A method used to run all steps of transformation.
			*/
			virtual void run();

			template <typename cv_mat_type>
			std::vector<Polygon_with_attributes> estimate_shift_one_block(OGREnvelope& c_block, const TemplateMatchingMethod& matching_strategy) {

				// Loading blocks of search and template image
				RasterIO search_raster(params->search_image, GA_ReadOnly, false), template_raster(params->template_image, GA_ReadOnly, false);
				GeoImage<cv_mat_type> template_geoimage(template_raster.raster_data, template_raster.geotransform);
				GeoImage<cv_mat_type> search_geoimage(search_raster.raster_data, search_raster.geotransform);
				
				TMStrategy<cv_mat_type> tms(matching_strategy);
				double search_radius_meters = search_raster.get_pixel_height() * search_radius_pixels;

				// Loading filtered view of shapefile
				PolygonsShapfileIO target_shape;
				OGRPolygon c_block_polygon = envelopeToPolygon(c_block);
				target_shape.load_shapefile(params->template_shapefile, false, &c_block_polygon);
				// Rasterizing masks for every geometry
				//auto geometries_masks = rasterize_geometries_masks<Boost_Polygon_2, cv_mat_type>(target_shape.geometries_container, search_raster);

				/************** Creating soaked geometries *************/
				std::vector<Soacked_Pixels_Boost_Polygon_2<cv_mat_type>> soaked_geometries; soaked_geometries.reserve(target_shape.geometries_container.size());
				for (size_t geometry_idx = 0; geometry_idx < target_shape.geometries_container.size(); geometry_idx++) {
					auto& c_geometry = target_shape.geometries_container[geometry_idx];
					Soacked_Pixels_Boost_Polygon_2<cv_mat_type> c_soaked_geom; c_soaked_geom.inners().resize(c_geometry.inners().size());
					boost::geometry::assign(c_soaked_geom, c_geometry);

					Boost_Box_2 c_geometry_bounds;
					boost::geometry::envelope(c_geometry.outer(), c_geometry_bounds);
					GeoImage<cv_mat_type> c_geom_image = template_geoimage.get_view_spatial<cv_mat_type>(
						c_geometry_bounds.min_corner().get<0>(), c_geometry_bounds.min_corner().get<1>(),
						c_geometry_bounds.max_corner().get<0>(), c_geometry_bounds.max_corner().get<1>()
						);

					c_soaked_geom.outer_ring_soak.img = c_geom_image.image;
					//c_soaked_geom.outer_ring_soak.mask = geometries_masks[geometry_idx].image;
					soaked_geometries.push_back(c_soaked_geom);
				}

				/******************* Estimating shift by template matching ******************/
				tqdm progress_bar;
				std::vector<Polygon_with_attributes> shift_estimated_polygons; shift_estimated_polygons.reserve(target_shape.geometries_container.size());
				for (size_t geometry_idx = 0; geometry_idx < target_shape.geometries_container.size(); geometry_idx++) {
					progress_bar.progress(geometry_idx, target_shape.geometries_container.size());
					Boost_Box_2 c_geometry_envelope; bg::envelope(target_shape.geometries_container[geometry_idx], c_geometry_envelope);					
					GeoImage<cv_mat_type> search_geoimage_crop = search_geoimage.get_view_spatial<cv_mat_type>(
						-search_radius_meters + c_geometry_envelope.min_corner().get<0>(),
						-search_radius_meters + c_geometry_envelope.min_corner().get<1>(),
						search_radius_meters + c_geometry_envelope.max_corner().get<0>(),
						search_radius_meters + c_geometry_envelope.max_corner().get<1>()
						);
					auto& template_image_crop = soaked_geometries[geometry_idx].outer_ring_soak.img;
					auto search_result = tms.compute(search_geoimage_crop.image, template_image_crop);
					auto best_val_pos = tms.get_best(search_result);
					Polygon_with_attributes out_polygon(target_shape.geometries_container[geometry_idx]);
					out_polygon.set_double_attribute("dx", (best_val_pos.second.x - search_radius_pixels) * search_raster.geotransform[1]);
					out_polygon.set_double_attribute("dy", (best_val_pos.second.y - search_radius_pixels) * search_raster.geotransform[5]);
					out_polygon.set_double_attribute("tm_val", best_val_pos.first);
					shift_estimated_polygons.push_back(out_polygon);
					/*** Temp*/
					cv::namedWindow("img_display", cv::WINDOW_AUTOSIZE); cv::Mat img_display;
					cv::namedWindow("result_display", cv::WINDOW_AUTOSIZE); cv::Mat result_display;
					cv::namedWindow("template_display", cv::WINDOW_AUTOSIZE); cv::Mat template_display;
					if constexpr (std::is_same_v<cv_mat_type, cv::Mat>) {
						img_display = search_geoimage_crop.image;
						result_display = search_result;
						template_display = template_image_crop;
					}
					else {
						search_geoimage_crop.image.download(img_display);
						search_result.download(result_display);
						template_image_crop.download(template_display);
					}
					cv::rectangle(img_display, best_val_pos.second, cv::Point(best_val_pos.second.x + template_image_crop.cols, best_val_pos.second.y + template_image_crop.rows), cv::Scalar(10,80,180), 2, 2, 0);
					cv::rectangle(result_display, best_val_pos.second, cv::Point(best_val_pos.second.x + template_image_crop.cols, best_val_pos.second.y + template_image_crop.rows), cv::Scalar(10, 80, 180), 2, 2, 0);

					cv::imshow("img_display", img_display);
					cv::imshow("template_display", template_display);
					cv::imshow("result_display", result_display);
					cv::waitKey(0);
					//***/
				}
				progress_bar.finish();
				return shift_estimated_polygons;
			}


			template <typename cv_mat_type>
			std::vector<Polygon_with_attributes> estimate_shift_one_block_1d(OGREnvelope& c_block, const TemplateMatchingMethod& matching_strategy, const double& rotation_angle) {

				RasterIO search_raster;
				search_raster.load_raster(params->search_image);

				auto template_geoimage = rotate(GeoImage<cv_mat_type>::from_file(params->template_image, c_block), rotation_angle);
				auto search_geoimage = rotate(GeoImage<cv_mat_type>::from_file(params->search_image, c_block), rotation_angle);

				//template_geoimage.to_file("C:/Users/geoimage/Music/t_copy.tif", search_raster.spatial_refrence);
				//search_geoimage.to_file("C:/Users/geoimage/Music/s_copy.tif", search_raster.spatial_refrence);

				TMStrategy<cv_mat_type> tms(matching_strategy);

				PolygonsShapfileIO target_shape;
				OGRPolygon c_block_polygon = envelopeToPolygon(c_block);
				target_shape.load_shapefile(params->template_shapefile, false, &c_block_polygon);

				bg::strategy::transform::matrix_transformer<double, 2, 2> trans_obj(
					template_geoimage.geotransform[1], template_geoimage.geotransform[2], template_geoimage.geotransform[0],
					template_geoimage.geotransform[4], template_geoimage.geotransform[5], template_geoimage.geotransform[3],
					0.0,0.0,1.0
				);
				bg::strategy::transform::inverse_transformer<double, 2, 2> coords_to_pix_transformer(trans_obj);
				/************** Creating soaked geometries *************/
				std::vector<Soacked_Pixels_Boost_Polygon_2<cv_mat_type>> soaked_geometries; soaked_geometries.reserve(target_shape.geometries_container.size());
				for (size_t geometry_idx = 0; geometry_idx < target_shape.geometries_container.size(); geometry_idx++) {
					auto& c_geometry = target_shape.geometries_container[geometry_idx];
					Soacked_Pixels_Boost_Polygon_2<cv_mat_type> c_soaked_geom; c_soaked_geom.inners().resize(c_geometry.inners().size());
					boost::geometry::transform(c_geometry, c_soaked_geom, coords_to_pix_transformer);
					//boost::geometry::assign(c_soaked_geom, c_geometry);

					Boost_Box_2 c_geometry_bounds;
					boost::geometry::envelope(c_soaked_geom.outer(), c_geometry_bounds);
					GeoImage<cv_mat_type> c_geom_image = template_geoimage.get_view_pixel<cv_mat_type>(
						c_geometry_bounds.min_corner().get<0>(), c_geometry_bounds.min_corner().get<1>(),
						c_geometry_bounds.max_corner().get<0>()- c_geometry_bounds.min_corner().get<0>(),
						c_geometry_bounds.max_corner().get<1>() - c_geometry_bounds.min_corner().get<1>()
						);

					c_soaked_geom.outer_ring_soak.img = c_geom_image.image;
					//c_soaked_geom.outer_ring_soak.mask = geometries_masks[geometry_idx].image;
					soaked_geometries.push_back(c_soaked_geom);
				}

				tqdm progress_bar;
				std::vector<Polygon_with_attributes> shift_estimated_polygons; shift_estimated_polygons.reserve(target_shape.geometries_container.size());
				for (size_t geometry_idx = 0; geometry_idx < target_shape.geometries_container.size(); geometry_idx++) {
					progress_bar.progress(geometry_idx, target_shape.geometries_container.size());
					auto& c_soaked_geometry = soaked_geometries[geometry_idx];
					Boost_Box_2 c_geometry_bounds;
					boost::geometry::envelope(c_soaked_geometry.outer(), c_geometry_bounds);
					
					int geom_width_pixel = c_geometry_bounds.max_corner().get<0>() - c_geometry_bounds.min_corner().get<0>();
					int geom_height_pixel = c_geometry_bounds.max_corner().get<1>() - c_geometry_bounds.min_corner().get<1>();
					GeoImage<cv_mat_type> search_geoimage_crop = search_geoimage.get_view_pixel<cv_mat_type>(
						-search_radius_pixels + c_geometry_bounds.min_corner().get<0>(),
						-3 + couple_v_displacement+ c_geometry_bounds.min_corner().get<1>(),
						2*search_radius_pixels+ geom_width_pixel,
						3 + couple_v_displacement + geom_height_pixel
						);
					auto& template_image_crop = c_soaked_geometry.outer_ring_soak.img;
					auto search_result = tms.compute(search_geoimage_crop.image, template_image_crop);
					auto best_val_pos = tms.get_best(search_result);
					Polygon_with_attributes out_polygon(target_shape.geometries_container[geometry_idx]);
					/****NOT SURE CHECK IT AGAIN***/
					double epipolar_disp = (best_val_pos.second.x - search_radius_pixels + couple_v_displacement) * search_raster.geotransform[1];
					double x_weight = std::cos(RADS(-rotation_angle));
					double y_weight = std::sin(RADS(-rotation_angle));

					out_polygon.set_double_attribute("dx", epipolar_disp * x_weight);
					out_polygon.set_double_attribute("dy", epipolar_disp * y_weight);
					out_polygon.set_double_attribute("epi_x", best_val_pos.second.x);
					out_polygon.set_double_attribute("epi_y", best_val_pos.second.y);
					out_polygon.set_double_attribute("tm_val", best_val_pos.first);
					shift_estimated_polygons.push_back(out_polygon);
					/***
					cv::destroyAllWindows();
					cv::namedWindow("img_display", cv::WINDOW_AUTOSIZE); cv::Mat img_display;
					cv::namedWindow("result_display", cv::WINDOW_AUTOSIZE); cv::Mat result_display;
					cv::namedWindow("template_display", cv::WINDOW_AUTOSIZE); cv::Mat template_display;
					if constexpr (std::is_same_v<cv_mat_type, cv::Mat>) {
						img_display = search_geoimage_crop.image;
						result_display = search_result;
						template_display = template_image_crop;
					}
					else {
						search_geoimage_crop.image.download(img_display);
						search_result.download(result_display);
						template_image_crop.download(template_display);
					}
					cv::rectangle(img_display, best_val_pos.second, cv::Point(best_val_pos.second.x + template_image_crop.cols, best_val_pos.second.y + template_image_crop.rows), cv::Scalar(10, 80, 180), 2, 2, 0);
					cv::rectangle(result_display, best_val_pos.second, cv::Point(best_val_pos.second.x + template_image_crop.cols, best_val_pos.second.y + template_image_crop.rows), cv::Scalar(10, 80, 180), 2, 2, 0);
					
					cv::destroyAllWindows();
					cv::imshow("img_display", img_display);
					cv::imshow("template_display", template_display);
					cv::imshow("result_display", result_display);
					cv::waitKey(0);
					***/
				}
				progress_bar.finish();
				return shift_estimated_polygons;
			}


			/**
			*  A method to check requirements before running transformation steps.
			* Example: -Checking vector exsitance, checking output_path overwrite, check algorithm parameters ...
			* @return an bool indicating if can run algorithm
			*/
			bool pre_check();

			private:
				double couple_rotation_angle, couple_v_displacement = 0;
				int search_radius_pixels;

		};
	}
}