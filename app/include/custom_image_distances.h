#pragma once
#include "defs.h"

namespace LxGeo
{
	namespace templateMatchingAlignment
	{

        /******************* Peak signal to noise ratio ******************/
        double getPSNR(const cv::Mat& I1, const cv::Mat& I2);

        double getPSNR(const cv::cuda::GpuMat& gI1, const cv::cuda::GpuMat& gI2);


        struct BufferPSNR                                     // Optimized CUDA versions
        {   // Data allocations are very expensive on CUDA. Use a buffer to solve: allocate once reuse later.
            cv::cuda::GpuMat gs, t1, t2;
            cv::cuda::GpuMat buf;
        };
        double getPSNR_optimized(const cv::cuda::GpuMat& gI1, const cv::cuda::GpuMat& gI2, BufferPSNR& b);



        /****************** Structural similarity index *******************/
        cv::Scalar getMSSIM(const cv::Mat& i1, const cv::Mat& i2);

        cv::Scalar getMSSIM(const cv::cuda::GpuMat& gI1, const cv::cuda::GpuMat& gI2);

        struct BufferMSSIM                                     // Optimized CUDA versions
        {   // Data allocations are very expensive on CUDA. Use a buffer to solve: allocate once reuse later.
            cv::cuda::GpuMat gs, t1, t2;
            cv::cuda::GpuMat I1_2, I2_2, I1_I2;
            std::vector<cv::cuda::GpuMat> vI1, vI2;
            cv::cuda::GpuMat mu1, mu2;
            cv::cuda::GpuMat mu1_2, mu2_2, mu1_mu2;
            cv::cuda::GpuMat sigma1_2, sigma2_2, sigma12;
            cv::cuda::GpuMat t3;
            cv::cuda::GpuMat ssim_map;
            cv::cuda::GpuMat buf;
        };

        cv::Scalar getMSSIM_optimized(const cv::cuda::GpuMat& gI1, const cv::cuda::GpuMat& gI2, BufferMSSIM& b);

        /*
        template <typename _Tp>
        void OLBP(const cv::Mat& src, cv::Mat& dst) {
            dst = cv::Mat::zeros(src.rows - 2, src.cols - 2, CV_8UC1);
            for (int i = 1; i < src.rows - 1; i++) {
                for (int j = 1; j < src.cols - 1; j++) {
                    _Tp center = src.at<_Tp>(i, j);
                    unsigned char code = 0;
                    code |= (src.at<_Tp>(i - 1, j - 1) > center) << 7;
                    code |= (src.at<_Tp>(i - 1, j) > center) << 6;
                    code |= (src.at<_Tp>(i - 1, j + 1) > center) << 5;
                    code |= (src.at<_Tp>(i, j + 1) > center) << 4;
                    code |= (src.at<_Tp>(i + 1, j + 1) > center) << 3;
                    code |= (src.at<_Tp>(i + 1, j) > center) << 2;
                    code |= (src.at<_Tp>(i + 1, j - 1) > center) << 1;
                    code |= (src.at<_Tp>(i, j - 1) > center) << 0;
                    dst.at<unsigned char>(i - 1, j - 1) = code;
                }
            }
        }
        */

	}
}