#include <opencv2/opencv.hpp>
extern "C"
{
#include "Lib/Kernels/ref.h"
#include "Lib/Common/types.h"
}
#include "..//DemoEngine.h"
#include "time.h"

///////////////////////////////////////////////////////////////////////////////
class demo_Watershed : public IDemoCase
{
public:
   demo_Watershed(){}
   ///@see IDemoCase::ReplyName
   virtual std::string ReplyName() const override
   {
      return "Watershed";
   }
private:
   ///@see IDemoCase::execute
   virtual void execute() override;
   ///@brief provide interactive demo
   static void applyParameters(void* data);
private:
   ///@source image
   cv::Mat m_srcImage;
};
///////////////////////////////////////////////////////////////////////////////
namespace
{
   const std::string m_openVXWindow = "openVX";
   const std::string m_openCVWindow = "openCV";
   const std::string m_originalWindow = " original ";
   const std::string m_diffWindow = m_openVXWindow + "-" + m_openCVWindow;
}
///////////////////////////////////////////////////////////////////////////////
void demo_Watershed::execute()
{
   cv::destroyAllWindows();
   cv::namedWindow(m_originalWindow, CV_WINDOW_NORMAL);
   cv::namedWindow(m_openVXWindow, CV_WINDOW_NORMAL);
   cv::namedWindow(m_openCVWindow, CV_WINDOW_NORMAL);
   cv::namedWindow(m_diffWindow, CV_WINDOW_NORMAL);
   const std::string imgPath = "..\\Image\\apple.png";
   m_srcImage = cv::imread(imgPath, CV_LOAD_IMAGE_COLOR);
   cv::imshow(m_originalWindow, m_srcImage);
   applyParameters(this);
   cv::waitKey(0);
}
///////////////////////////////////////////////////////////////////////////////
void demo_Watershed::applyParameters(void* data)
{
   auto demo = static_cast<demo_Watershed*>(data);
   const cv::Size imgSize(demo->m_srcImage.cols, demo->m_srcImage.rows);  //width , height
   ///@{ OPENCV
   cv::Mat bw1;
   cv::cvtColor(demo->m_srcImage, bw1, CV_BGR2GRAY);
   threshold(bw1, bw1, 40, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
   cv::Mat dist1;
   distanceTransform(bw1, dist1, CV_DIST_L2, 3);
   // Normalize the distance image for range = {0.0, 1.0}
   // so we can visualize and threshold it
   normalize(dist1, dist1, 0, 1., cv::NORM_MINMAX);
   // Threshold to obtain the peaks
   // This will be the markers for the foreground objects
   threshold(dist1, dist1, .4, 1., CV_THRESH_BINARY);
   // Dilate a bit the dist image
   cv::Mat kernel1 = cv::Mat::ones(3, 3, CV_8UC1);
   dilate(dist1, dist1, kernel1);
   // Create the CV_8U version of the distance image
   // It is needed for findContours()
   cv::Mat dist1_8u;
   dist1.convertTo(dist1_8u, CV_8U);
   // Find total markers
   std::vector<std::vector<cv::Point> > contours1;
   findContours(dist1_8u, contours1, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
   // Create the marker image for the watershed algorithm
   cv::Mat markers1 = cv::Mat::zeros(dist1.size(), CV_32S);
   // Draw the foreground markers
   for (size_t i = 0; i < contours1.size(); i++)
      drawContours(markers1, contours1, static_cast<int>(i), cv::Scalar::all(static_cast<int>(i)+1), -1);
   // Draw the background marker
   circle(markers1, cv::Point(5, 5), 3, CV_RGB(255, 255, 255), -1);
   // Perform the watershed algorithm
   cv::Mat markersfromCVwatershed = cv::Mat(markers1.size(), markers1.type(), markers1.data);
   //var for count time work cv::watershed
   clock_t t1;
   // time before working cv::watershed
   t1 = clock();
   cv::watershed(demo->m_srcImage, markersfromCVwatershed);
   //time after working cv::watershed
   t1 = clock() - t1;
   std::cout << ((double)t1) / CLOCKS_PER_SEC << ": sec for cv::watershed" << std::endl;
   //show result image after cv::watershed
   cv::imshow(m_openCVWindow, markersfromCVwatershed * 100);
   ///@}
   ///@{ OPENVX
   cv::Mat bw2;
   cv::cvtColor(demo->m_srcImage, bw2, CV_BGR2GRAY);
   threshold(bw2, bw2, 40, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
   cv::Mat dist2;
   distanceTransform(bw2, dist2, CV_DIST_L2, 3);
   // Normalize the distance image for range = {0.0, 1.0}
   // so we can visualize and threshold it
   normalize(dist2, dist2, 0, 1., cv::NORM_MINMAX);
   // Threshold to obtain the peaks
   // This will be the markers for the foreground objects
   threshold(dist2, dist2, .4, 1., CV_THRESH_BINARY);
   // Dilate a bit the dist image
   cv::Mat kernel2 = cv::Mat::ones(3, 3, CV_8UC1);
   dilate(dist2, dist2, kernel2);
   // Create the CV_8U version of the distance image
   // It is needed for findContours()
   cv::Mat dist2_8u;
   dist2.convertTo(dist2_8u, CV_8U);
   // Find total markers
   std::vector<std::vector<cv::Point> > contours2;
   findContours(dist2_8u, contours2, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
   // Create the marker image for the watershed algorithm
   cv::Mat markers2 = cv::Mat::zeros(dist2.size(), CV_32S);
   // Draw the foreground markers
   for (size_t i = 0; i < contours2.size(); i++)
      drawContours(markers2, contours2, static_cast<int>(i), cv::Scalar::all(static_cast<int>(i)+1), -1);
   // Draw the background marker
   circle(markers2, cv::Point(5, 5), 3, CV_RGB(255, 255, 255), -1);
   cv::Mat srcImageForVX, maskImageForVX;
   demo->m_srcImage.convertTo(srcImageForVX, CV_8UC3);
   cv::cvtColor(srcImageForVX, srcImageForVX, CV_BGR2RGB);
   markers2.convertTo(maskImageForVX, CV_32S);
   _vx_image srcVXImage = {
      srcImageForVX.data,
      imgSize.width,
      imgSize.height,
      VX_DF_IMAGE_RGB,
      VX_COLOR_SPACE_DEFAULT
   };
   _vx_image markersVXImage = {
      maskImageForVX.data,
      imgSize.width,
      imgSize.height,
      VX_DF_IMAGE_S32,
      VX_COLOR_SPACE_DEFAULT
   };
   // Perform the watershed algorithm
   //var for count time work ref_WatershedSegmentation
   clock_t t;
   //time before working ref_WatershedSegmentation
   t = clock();
   // works ref_WatershedSegmentation
   ref_WatershedSegmentation(&srcVXImage, &markersVXImage);
   //time after working ref_WatershedSegmentation
   t = clock() - t;
   std::cout << (double(t)) / CLOCKS_PER_SEC << "sec for ref_watershed" << std::endl;
   // show image after ref_WatershedSegmentation
   cv::Mat markersfromVXwatershed = cv::Mat(markers2.size(), markers2.type(), markersVXImage.data);
   cv::imshow(m_openVXWindow, markersfromVXwatershed * 100);
   ///@}
   //// Show difference of OpenVX and OpenCV
   const cv::Mat diffImage(imgSize, CV_32S);
   cv::absdiff(markersfromVXwatershed, markersfromCVwatershed, diffImage);
   cv::imshow(m_diffWindow, diffImage);
}
///////////////////////////////////////////////////////////////////////////////
IDemoCasePtr CreateWatershedDemo()
{
   return std::make_unique<demo_Watershed>();
}
