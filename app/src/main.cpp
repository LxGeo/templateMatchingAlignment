#include "gdal.h"
#include "parameters.h"
#include "GDAL_OPENCV_IO.h"
#include "template_matcher_alignment.h"

using namespace LxGeo::templateMatchingAlignment;

int main(int argc, char* argv[])
{
	clock_t t_begin = clock();
	//GDALAllRegister();
	KGDAL2CV* kgdal2cv = new KGDAL2CV();
	
	params = new Parameters(argc, argv);
	if (!params->initialized()) {
		delete params;
		return 1;
	}
	
	// Runs process
	polygonsMatcher v_t_app = polygonsMatcher();
	if (v_t_app.pre_check())
		v_t_app.run();

	// Quits
	delete params;	
	delete kgdal2cv;
	

	clock_t t_end = clock();
	std::cout << "** Elapsed time : " << double(t_end - t_begin) / CLOCKS_PER_SEC << " s." << std::endl;

	return 0;
}