#pragma once
#include "defs.h"
#include "custom_image_distances.h"
#include "tm_ssim.h"
#include "immintrin.h"

namespace LxGeo
{
	namespace templateMatchingAlignment
	{

		enum TemplateMatchingMethod
		{
			sq_diff = 0, // Squared distance diffrence
			sq_diff_N = 1 << 0, // Squared distance diffrence normalized
			ccorr = 1 << 1, // Cross correaltion distance
			ccorr_N = 1 << 2, // Cross correaltion distance normalized
			ccoef = 1 << 3, // Pearson coeffecient
			ccoef_N = 1 << 4, // Pearson coeffecient normalized
			ssim = 1 << 5, // Structrual similarity index
			mu_inf_N = 1 << 6 // Mutual inforamtion normalized
		};

		static std::map<TemplateMatchingMethod,int> cv_tm_keys = {
			{TemplateMatchingMethod::sq_diff, cv::TM_SQDIFF},
			{TemplateMatchingMethod::sq_diff_N, cv::TM_SQDIFF_NORMED},
			{TemplateMatchingMethod::ccorr, cv::TM_CCORR},
			{TemplateMatchingMethod::ccorr_N, cv::TM_CCORR_NORMED},
			{TemplateMatchingMethod::ccoef, cv::TM_CCOEFF},
			{TemplateMatchingMethod::ccoef_N, cv::TM_CCOEFF_NORMED},
		};

		// Base class for template matching strategy

		template <typename cv_mat_type>
		class TMStrategy {

		public:

			TemplateMatchingMethod matching_strategy;
			
			TMStrategy(const TemplateMatchingMethod& _matching_strategy): matching_strategy(_matching_strategy){}

			cv_mat_type compute(const cv_mat_type& search_image, const cv_mat_type& template_image) {
				cv_mat_type result_matrix;
				int result_cols = search_image.cols - template_image.cols + 1;
				int result_rows = search_image.rows - template_image.rows + 1;
				result_matrix.create(result_rows, result_cols, CV_32FC1);
				// Switch method strategy
				auto cv_key = cv_tm_keys.find(matching_strategy);
				if (cv_key!=cv_tm_keys.end()){
					compute_cv_available(search_image, template_image, result_matrix, cv_key->second);
					return result_matrix;
				}
				else if (matching_strategy == TemplateMatchingMethod::ssim) {
					compute_ssim(search_image, template_image, result_matrix);
					return result_matrix;
				}
			}

			template <typename cv_mat_type>
			std::pair<double, cv::Point> get_best(const cv_mat_type& match_array) {
				double min_val, max_val;
				cv::Point min_loc, max_loc;
				if constexpr (std::is_same_v<cv_mat_type, cv::Mat>)
					cv::minMaxLoc(match_array, &min_val, &max_val, &min_loc, &max_loc);
				else
					cv::cuda::minMaxLoc(match_array, &min_val, &max_val, &min_loc, &max_loc);
				// if best value is maximum value
				if (
					matching_strategy == TemplateMatchingMethod::sq_diff ||
					matching_strategy == TemplateMatchingMethod::sq_diff_N
					) {
					return std::make_pair(min_val, min_loc);
				}
				else
					return std::make_pair(max_val, max_loc);

			}


		private:
			void compute_cv_available(const cv_mat_type& search_image, const cv_mat_type& template_image, cv_mat_type& result, int method_name) {
				if constexpr (std::is_same_v<cv_mat_type, cv::Mat>) {
					cv::matchTemplate(search_image, template_image, result, method_name);
				}
				else {
					auto matcher = cv::cuda::createTemplateMatching(search_image.type(), method_name);
					matcher->match(search_image, template_image, result);
				}
			}

			void compute_ssim(const cv_mat_type& search_image, const cv_mat_type& template_image, cv_mat_type& result) {
				if constexpr (std::is_same_v<cv_mat_type, cv::Mat>) {
					for (int c_col = 0; c_col < result.cols; c_col++) {
						for (int c_row = 0; c_row < result.rows; c_row++) {
							result.at<float>(c_row, c_col) = getMSSIM(
								search_image(cv::Rect(c_col, c_row, template_image.cols, template_image.rows)),
								template_image
							)[0];
						}
					}
				}
				else {
					cv::Mat aux(result.size(), result.type());
					BufferMSSIM bufferMSSIM;
					for (int c_col = 0; c_col < result.cols; c_col++) {
						for (int c_row = 0; c_row < result.rows; c_row++) {
							aux.at<float>(c_row, c_col) = getMSSIM_optimized(
								search_image(cv::Rect(c_col, c_row, template_image.cols, template_image.rows)),
								template_image,
								bufferMSSIM
							)[0];
						}
					}
					result.upload(aux);
					//template_match_ssim(search_image, template_image, result);
				}
			}

		};

	}
}