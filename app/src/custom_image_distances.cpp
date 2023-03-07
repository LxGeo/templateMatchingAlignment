#include "custom_image_distances.h"

namespace LxGeo
{
    namespace templateMatchingAlignment
    {

        /******************* Peak signal to noise ratio ******************/
        double getPSNR(const cv::Mat& I1, const cv::Mat& I2)
        {
            cv::Mat s1;
            absdiff(I1, I2, s1);       // |I1 - I2|
            s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
            s1 = s1.mul(s1);           // |I1 - I2|^2
            cv::Scalar s = cv::sum(s1);         // sum elements per channel
            double sse = s.val[0] + s.val[1] + s.val[2]; // sum channels
            if (sse <= 1e-10) // for small values return zero
                return 0;
            else
            {
                double  mse = sse / (double)(I1.channels() * I1.total());
                double psnr = 10.0 * log10((255 * 255) / mse);
                return psnr;
            }
        }

        double getPSNR(const cv::cuda::GpuMat& gI1, const cv::cuda::GpuMat& gI2)
        {
            cv::cuda::GpuMat gs, t1, t2;
            gI1.convertTo(t1, CV_32F);
            gI2.convertTo(t2, CV_32F);
            cv::cuda::absdiff(t1.reshape(1), t2.reshape(1), gs);
            cv::cuda::multiply(gs, gs, gs);
            cv::Scalar s = cv::cuda::sum(gs);
            double sse = s.val[0] + s.val[1] + s.val[2];
            if (sse <= 1e-10) // for small values return zero
                return 0;
            else
            {
                double  mse = sse / (double)(gI1.channels() * gI1.rows * gI1.cols);
                double psnr = 10.0 * log10((255 * 255) / mse);
                return psnr;
            }
        }

        double getPSNR_optimized(const cv::cuda::GpuMat& gI1, const cv::cuda::GpuMat& gI2, BufferPSNR& b)
        {
            gI1.convertTo(b.t1, CV_32F);
            gI2.convertTo(b.t2, CV_32F);
            cv::cuda::absdiff(b.t1.reshape(1), b.t2.reshape(1), b.gs);
            cv::cuda::multiply(b.gs, b.gs, b.gs);
            double sse = cv::cuda::sum(b.gs, b.buf)[0];
            if (sse <= 1e-10) // for small values return zero
                return 0;
            else
            {
                double mse = sse / (double)(gI1.channels() * gI1.rows * gI1.cols);
                double psnr = 10.0 * log10((255 * 255) / mse);
                return psnr;
            }
        }

        /****************** Structural similarity index *******************/
        cv::Scalar getMSSIM(const cv::Mat& i1, const cv::Mat& i2)
        {
            const double C1 = 6.5025, C2 = 58.5225;
            /***************************** INITS **********************************/
            int d = CV_32F;
            cv::Mat I1, I2;
            i1.convertTo(I1, d);           // cannot calculate on one byte large values
            i2.convertTo(I2, d);
            cv::Mat I2_2 = I2.mul(I2);        // I2^2
            cv::Mat I1_2 = I1.mul(I1);        // I1^2
            cv::Mat I1_I2 = I1.mul(I2);        // I1 * I2
            /*************************** END INITS **********************************/
            cv::Mat mu1, mu2;   // PRELIMINARY COMPUTING
            GaussianBlur(I1, mu1, cv::Size(11, 11), 1.5);
            GaussianBlur(I2, mu2, cv::Size(11, 11), 1.5);
            cv::Mat mu1_2 = mu1.mul(mu1);
            cv::Mat mu2_2 = mu2.mul(mu2);
            cv::Mat mu1_mu2 = mu1.mul(mu2);
            cv::Mat sigma1_2, sigma2_2, sigma12;
            GaussianBlur(I1_2, sigma1_2, cv::Size(11, 11), 1.5);
            sigma1_2 -= mu1_2;
            GaussianBlur(I2_2, sigma2_2, cv::Size(11, 11), 1.5);
            sigma2_2 -= mu2_2;
            GaussianBlur(I1_I2, sigma12, cv::Size(11, 11), 1.5);
            sigma12 -= mu1_mu2;
            cv::Mat t1, t2, t3;
            t1 = 2 * mu1_mu2 + C1;
            t2 = 2 * sigma12 + C2;
            t3 = t1.mul(t2);              // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))
            t1 = mu1_2 + mu2_2 + C1;
            t2 = sigma1_2 + sigma2_2 + C2;
            t1 = t1.mul(t2);               // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))
            cv::Mat ssim_map;
            divide(t3, t1, ssim_map);      // ssim_map =  t3./t1;
            cv::Scalar mssim = mean(ssim_map); // mssim = average of ssim map
            return mssim;
        }

        cv::Scalar getMSSIM(const cv::cuda::GpuMat& gI1, const cv::cuda::GpuMat& gI2)
        {
            const float C1 = 6.5025f, C2 = 58.5225f;
            /***************************** INITS **********************************/
            cv::cuda::GpuMat gs1, tmp1, tmp2;
            gI1.convertTo(tmp1, CV_MAKE_TYPE(CV_32F, gI1.channels()));
            gI2.convertTo(tmp2, CV_MAKE_TYPE(CV_32F, gI2.channels()));
            std::vector<cv::cuda::GpuMat> vI1, vI2;
            cv::cuda::split(tmp1, vI1);
            cv::cuda::split(tmp2, vI2);
            cv::Scalar mssim;
            cv::Ptr<cv::cuda::Filter> gauss = cv::cuda::createGaussianFilter(vI2[0].type(), -1, cv::Size(11, 11), 1.5);
            for (int i = 0; i < gI1.channels(); ++i)
            {
                cv::cuda::GpuMat I2_2, I1_2, I1_I2;
                cv::cuda::multiply(vI2[i], vI2[i], I2_2);        // I2^2
                cv::cuda::multiply(vI1[i], vI1[i], I1_2);        // I1^2
                cv::cuda::multiply(vI1[i], vI2[i], I1_I2);       // I1 * I2
                /*************************** END INITS **********************************/
                cv::cuda::GpuMat mu1, mu2;   // PRELIMINARY COMPUTING
                gauss->apply(vI1[i], mu1);
                gauss->apply(vI2[i], mu2);
                cv::cuda::GpuMat mu1_2, mu2_2, mu1_mu2;
                cv::cuda::multiply(mu1, mu1, mu1_2);
                cv::cuda::multiply(mu2, mu2, mu2_2);
                cv::cuda::multiply(mu1, mu2, mu1_mu2);
                cv::cuda::GpuMat sigma1_2, sigma2_2, sigma12;
                gauss->apply(I1_2, sigma1_2);
                cv::cuda::subtract(sigma1_2, mu1_2, sigma1_2); // sigma1_2 -= mu1_2;
                gauss->apply(I2_2, sigma2_2);
                cv::cuda::subtract(sigma2_2, mu2_2, sigma2_2); // sigma2_2 -= mu2_2;
                gauss->apply(I1_I2, sigma12);
                cv::cuda::subtract(sigma12, mu1_mu2, sigma12); // sigma12 -= mu1_mu2;
                cv::cuda::GpuMat t1, t2, t3;
                mu1_mu2.convertTo(t1, -1, 2, C1); // t1 = 2 * mu1_mu2 + C1;
                sigma12.convertTo(t2, -1, 2, C2); // t2 = 2 * sigma12 + C2;
                cv::cuda::multiply(t1, t2, t3);        // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))
                cv::cuda::addWeighted(mu1_2, 1.0, mu2_2, 1.0, C1, t1);       // t1 = mu1_2 + mu2_2 + C1;
                cv::cuda::addWeighted(sigma1_2, 1.0, sigma2_2, 1.0, C2, t2); // t2 = sigma1_2 + sigma2_2 + C2;
                cv::cuda::multiply(t1, t2, t1);                              // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))
                cv::cuda::GpuMat ssim_map;
                cv::cuda::divide(t3, t1, ssim_map);      // ssim_map =  t3./t1;
                cv::Scalar s = cv::cuda::sum(ssim_map);
                mssim.val[i] = s.val[0] / (ssim_map.rows * ssim_map.cols);
            }
            return mssim;
        }

        cv::Scalar getMSSIM_optimized(const cv::cuda::GpuMat& gI1, const cv::cuda::GpuMat& gI2, BufferMSSIM& b)
        {
            const float C1 = 6.5025f, C2 = 58.5225f;
            /***************************** INITS **********************************/
            cv::cuda::Stream stream;
            gI1.convertTo(b.t1, CV_32F, stream);
            gI2.convertTo(b.t2, CV_32F, stream);
            cv::cuda::split(b.t1, b.vI1, stream);
            cv::cuda::split(b.t2, b.vI2, stream);
            cv::Scalar mssim;
            cv::Ptr<cv::cuda::Filter> gauss = cv::cuda::createGaussianFilter(b.vI1[0].type(), -1, cv::Size(11, 11), 1.5);
            for (int i = 0; i < gI1.channels(); ++i)
            {
                cv::cuda::multiply(b.vI2[i], b.vI2[i], b.I2_2, 1, -1, stream);        // I2^2
                cv::cuda::multiply(b.vI1[i], b.vI1[i], b.I1_2, 1, -1, stream);        // I1^2
                cv::cuda::multiply(b.vI1[i], b.vI2[i], b.I1_I2, 1, -1, stream);       // I1 * I2
                gauss->apply(b.vI1[i], b.mu1, stream);
                gauss->apply(b.vI2[i], b.mu2, stream);
                cv::cuda::multiply(b.mu1, b.mu1, b.mu1_2, 1, -1, stream);
                cv::cuda::multiply(b.mu2, b.mu2, b.mu2_2, 1, -1, stream);
                cv::cuda::multiply(b.mu1, b.mu2, b.mu1_mu2, 1, -1, stream);
                gauss->apply(b.I1_2, b.sigma1_2, stream);
                cv::cuda::subtract(b.sigma1_2, b.mu1_2, b.sigma1_2, cv::cuda::GpuMat(), -1, stream);
                //b.sigma1_2 -= b.mu1_2;  - This would result in an extra data transfer operation
                gauss->apply(b.I2_2, b.sigma2_2, stream);
                cv::cuda::subtract(b.sigma2_2, b.mu2_2, b.sigma2_2, cv::cuda::GpuMat(), -1, stream);
                //b.sigma2_2 -= b.mu2_2;
                gauss->apply(b.I1_I2, b.sigma12, stream);
                cv::cuda::subtract(b.sigma12, b.mu1_mu2, b.sigma12, cv::cuda::GpuMat(), -1, stream);
                //b.sigma12 -= b.mu1_mu2;
                //here too it would be an extra data transfer due to call of operator*(Scalar, Mat)
                cv::cuda::multiply(b.mu1_mu2, 2, b.t1, 1, -1, stream); //b.t1 = 2 * b.mu1_mu2 + C1;
                cv::cuda::add(b.t1, C1, b.t1, cv::cuda::GpuMat(), -1, stream);
                cv::cuda::multiply(b.sigma12, 2, b.t2, 1, -1, stream); //b.t2 = 2 * b.sigma12 + C2;
                cv::cuda::add(b.t2, C2, b.t2, cv::cuda::GpuMat(), -12, stream);
                cv::cuda::multiply(b.t1, b.t2, b.t3, 1, -1, stream);     // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))
                cv::cuda::add(b.mu1_2, b.mu2_2, b.t1, cv::cuda::GpuMat(), -1, stream);
                cv::cuda::add(b.t1, C1, b.t1, cv::cuda::GpuMat(), -1, stream);
                cv::cuda::add(b.sigma1_2, b.sigma2_2, b.t2, cv::cuda::GpuMat(), -1, stream);
                cv::cuda::add(b.t2, C2, b.t2, cv::cuda::GpuMat(), -1, stream);
                cv::cuda::multiply(b.t1, b.t2, b.t1, 1, -1, stream);     // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))
                cv::cuda::divide(b.t3, b.t1, b.ssim_map, 1, -1, stream);      // ssim_map =  t3./t1;
                stream.waitForCompletion();
                cv::Scalar s = cv::cuda::sum(b.ssim_map, b.buf);
                mssim.val[i] = s.val[0] / (b.ssim_map.rows * b.ssim_map.cols);
            }
            return mssim;
        }

        

    }
}