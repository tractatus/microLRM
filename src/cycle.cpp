#include<opencv2/opencv.hpp>

#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <cmath>        // std::abs

#include <vector>

#include<iostream>
using namespace std;
using namespace cv;

#include "cycle.hpp"

ImageCycle::ImageCycle()
{

}

ImageCycle::~ImageCycle()
{

}

void ImageCycle::load_image (const String & filename, string folder, int start, int stop, int leadingzeros) {
    //output showing that it is loading
    //std::cout << "LOADING RAW IMAGES" << std::endl;
    if (stop == 0){
        cv::imreadmulti(filename, multipageImage, cv::IMREAD_ANYDEPTH);
    }else{
        for (int i = start; i <= stop; i++) {
            std::vector<Mat> B;
            std::string zplane = std::to_string(i);
            std::string newfilename = folder + "/" + filename + std::string(leadingzeros - zplane.length(), '0') + zplane + ".tif";
            cv::imreadmulti(newfilename, B, cv::IMREAD_ANYDEPTH);
            
            multipageImage.reserve(multipageImage.size() + B.size());                // preallocate memory without erase original data
            multipageImage.insert(multipageImage.end(), B.begin(), B.end());         // add B;
        }
        
    }
}

void ImageCycle::load_mask (const String & filename, string folder, int start, int stop, bool seperate) {
    //std::cout << "LOADING MASK IMAGES" << std::endl;
    if (seperate){
        cv::imreadmulti(filename, multipageMask, cv::IMREAD_ANYDEPTH);
    }else{
        std::string channelImage;
        std::vector<Mat> A;
        channelImage = folder + "/labels_" + std::to_string(0) + "_" + filename;
        cv::imreadmulti(channelImage, A, cv::IMREAD_ANYDEPTH);
        std::vector<Mat> B;
        channelImage = folder + "/labels_" + std::to_string(1) + "_" + filename;
        cv::imreadmulti(channelImage, B, cv::IMREAD_ANYDEPTH);
        std::vector<Mat> C;
        channelImage = folder + "/labels_" + std::to_string(2) + "_" + filename;
        cv::imreadmulti(channelImage, C, cv::IMREAD_ANYDEPTH);
        std::vector<Mat> D;
        channelImage = folder + "/labels_" + std::to_string(3) + "_" + filename;
        cv::imreadmulti(channelImage, D, cv::IMREAD_ANYDEPTH);
        stop = stop - start;
        start = start - start;
        for (int i = start; i <= stop; i++) {
            multipageMask.push_back(A[i]); //removed subtraction of 32768
            multipageMask.push_back(B[i]); //removed subtraction of 32768
            multipageMask.push_back(C[i]); //removed subtraction of 32768
            multipageMask.push_back(D[i]); //removed subtraction of 32768
        }
    }
    
}

void ImageCycle::save_image (const String & filename) {
  cv::imwrite(filename, multipageImage);
}

void ImageCycle::hist_norm (bool display) {
  entropy.clear();
  cdfs.clear();
  cv::Mat src_hist = cv::Mat::zeros(1, 65536, CV_64F); 
  cv::Mat src_cdf = cv::Mat::zeros(1, 65536, CV_64F);
  counter++;

  //std::cout << "dat" << counter << "<- list(";
  for(int i = 0; i < 4; i++){
   	src_hist.setTo(cv::Scalar(0));
  	src_cdf.setTo(cv::Scalar(0)); 	
	double* _src_hist = src_hist.ptr<double>();
	double* _src_cdf = src_cdf.ptr<double>();

	cv::Mat channel;
	cv::Mat mask;
	for (int l = i; l < multipageImage.size(); l = l + 4) {
    	channel.push_back(multipageImage[l]);
    	mask.push_back(multipageMask[l]);
  	}

	this->histogram(channel, mask, _src_hist, _src_cdf, entropy);
	cdfs.push_back(src_cdf.clone());
	//std::cout << "ch0" << i << " = c(";
  	for(int j = 0; j  < src_cdf.cols; j++){
  		

  		if(j == (src_cdf.cols - 1) ){
  			//std::cout << cdfs[i].at<double>(0,j);//_src_cdf[j];
  		}else{
  			//std::cout << cdfs[i].at<double>(0,j) << ", ";
  		}
  	}
  	
  	if(i == 3){
  		//std::cout << ")" << std::endl;
  	}else{
  		//std::cout << "), " << std::endl;
  	}
   
   }

   //std::cout << ", entropy" << " = c(";
   
   for(int k = 0; k < entropy.size(); k++){
   	if(k == 3){
  		//std::cout << entropy[k];
  	}else{
  		//std::cout << entropy[k] << ",";
  	}
   }


   //std::cout << "), cdfs = " <<  cdfs[0].cols << ")" << std::endl;

   refChannel = std::max_element(entropy.begin(),entropy.end()) - entropy.begin(); 



 }

 void ImageCycle::getIntensities(const std::vector<cv::Mat>& multipageImage, const std::vector<cv::Mat>& multipageMask){

 	int pixels = multipageImage[0].total() * multipageImage.size()/4;

	int ampliconID[4][65536] = {0};
	int noiseID[4][65536] = {0};
	int perimeterVolume[4][65536] = {0};
	//double Z[4][65536] = {0.0};
  	//double X[4][65536] = {0.0};
  	//double Y[4][65536] = {0.0};

	int strides[4][4] ={ {0, 1, 2, 3}, {-1, 0, 1, 2}, {-2, -1, 0, 1}, {-3, -2, -1, 0} };

	//Initialize m
	double minVal; 
	double maxVal; 
	cv::Point minLoc; 
	cv::Point maxLoc;

	int numberOfobjects = 0;
	for(int ch = 0; ch < 4; ch++){
		for (int l = ch; l < multipageImage.size(); l = l + 4) {
			cv::minMaxLoc(multipageMask[l], &minVal, &maxVal, &minLoc, &maxLoc );
			if(numberOfobjects < maxVal){
				numberOfobjects = (int)maxVal;
			}
		}
	}

	double X[4][numberOfobjects];
  	double Y[4][numberOfobjects];
  	double Z[4][numberOfobjects];
  	for(int ch = 0; ch < 4; ch++){
  	for(int i = 0; i < numberOfobjects; i++){
    	X[ch][i] = 0;
    	Y[ch][i] = 0;
    	Z[ch][i] = 0;
  	}
  	}

	for(int ch = 0; ch < 4; ch++){
	for (int l = ch; l < multipageImage.size(); l = l + 4) {
		for(int i = 0; i < multipageMask[0].rows; i++){
			for(int j = 0; j < multipageMask[0].cols; j++){
				ushort m = multipageMask[l].at<ushort>(i, j);
				ushort bs = backgroundSignal[l].at<ushort>(i, j);
				if(m > 0){
					pixelVolume[ch][m]++;
					double zvalue = (double)(l/4.0)/(multipageImage.size()/4);
					Z[ch][m] += zvalue;
					double xvalue = (double)j/multipageMask[0].cols;
					X[ch][m] += xvalue;
					double yvalue = (double)i/multipageMask[0].rows;
					Y[ch][m] += yvalue;
    
					if(ampliconID[ch][m] == 0){
						detectedAmplicons++;
					 	

					 	std::vector<int> values;

					 	values.push_back(multipageImage[l+strides[ch][0]].at<ushort>(i,j));
					 	values.push_back(multipageImage[l+strides[ch][1]].at<ushort>(i,j));
					 	values.push_back(multipageImage[l+strides[ch][2]].at<ushort>(i,j));
					 	values.push_back(multipageImage[l+strides[ch][3]].at<ushort>(i,j));

					 	ampliconMax.push_back( values );

					 	ampliconID[ch][m] = (int)ampliconMax.size() - 1;
					 	

					}else{

					 	ampliconMax[ampliconID[ch][m]].at(0) = ampliconMax[ampliconID[ch][m]].at(0) + multipageImage[l+strides[ch][0]].at<ushort>(i,j);
					 	ampliconMax[ampliconID[ch][m]].at(1) = ampliconMax[ampliconID[ch][m]].at(1) + multipageImage[l+strides[ch][1]].at<ushort>(i,j);
					 	ampliconMax[ampliconID[ch][m]].at(2) = ampliconMax[ampliconID[ch][m]].at(2) + multipageImage[l+strides[ch][2]].at<ushort>(i,j);
					 	ampliconMax[ampliconID[ch][m]].at(3) = ampliconMax[ampliconID[ch][m]].at(3) + multipageImage[l+strides[ch][3]].at<ushort>(i,j);
					 	
					}
					
				}//end checking mask pixel

				if(bs > 0){
					perimeterVolume[ch][bs]++;
					if(noiseID[ch][bs] == 0){
					 	

					 	std::vector<int> values;

					 	values.push_back(multipageImage[l+strides[ch][0]].at<ushort>(i,j));
					 	values.push_back(multipageImage[l+strides[ch][1]].at<ushort>(i,j));
					 	values.push_back(multipageImage[l+strides[ch][2]].at<ushort>(i,j));
					 	values.push_back(multipageImage[l+strides[ch][3]].at<ushort>(i,j));

					 	noise.push_back( values );

					 	noiseID[ch][m] = (int)noise.size() - 1;
					 	

					}else{


					 	noise[noiseID[ch][bs]].at(0) = noise[noiseID[ch][bs]].at(0) + multipageImage[l+strides[ch][0]].at<ushort>(i,j);
					 	noise[noiseID[ch][bs]].at(1) = noise[noiseID[ch][bs]].at(1) + multipageImage[l+strides[ch][1]].at<ushort>(i,j);
					 	noise[noiseID[ch][bs]].at(2) = noise[noiseID[ch][bs]].at(2) + multipageImage[l+strides[ch][2]].at<ushort>(i,j);
					 	noise[noiseID[ch][bs]].at(3) = noise[noiseID[ch][bs]].at(3) + multipageImage[l+strides[ch][3]].at<ushort>(i,j);

					 	
					}
					
				}//end checking mask pixel



			}
		}
		

  	}
  	}
     
  	std::vector<arma::mat> channel;
	std::cout << "intens = c(" << std::endl;
	for(int ch = 0; ch < 4; ch++){
		arma::mat singleChannel;
		arma::mat rowMat;
	for(int c=0; c < 65536; c++){
			if(pixelVolume[ch][c] > 0){
				rowMat << (double)Z[ch][c]/pixelVolume[ch][c] << (double)Y[ch][c]/pixelVolume[ch][c] << (double)X[ch][c]/pixelVolume[ch][c] << arma::endr;

				std::cout << ch << ", ";
				std::cout << pixelVolume[ch][c] << ", ";
				std::cout << ampliconMax[ampliconID[ch][c]].at(0)/pixelVolume[ch][c] << ", ";
				std::cout << ampliconMax[ampliconID[ch][c]].at(1)/pixelVolume[ch][c] << ", ";
				std::cout << ampliconMax[ampliconID[ch][c]].at(2)/pixelVolume[ch][c] << ", ";
				std::cout << ampliconMax[ampliconID[ch][c]].at(3)/pixelVolume[ch][c] << ", ";
                
                if(perimeterVolume[ch][c] == 0){
                    std::cout << 0 << ", ";
                    std::cout << 0 << ", ";
                    std::cout << 0 << ", ";
                    std::cout << 0 << ", ";
                }else{
                    std::cout << noise[noiseID[ch][c]].at(0)/perimeterVolume[ch][c] << ", ";
                    std::cout << noise[noiseID[ch][c]].at(1)/perimeterVolume[ch][c] << ", ";
                    std::cout << noise[noiseID[ch][c]].at(2)/perimeterVolume[ch][c] << ", ";
                    std::cout << noise[noiseID[ch][c]].at(3)/perimeterVolume[ch][c] << ", ";
                }

				std::cout << (double)X[ch][c]/pixelVolume[ch][c] << ", ";
				std::cout << (double)Y[ch][c]/pixelVolume[ch][c] << ", ";
				std::cout << (double)Z[ch][c]/pixelVolume[ch][c] << ", ";
				
				singleChannel = arma::join_cols(singleChannel, rowMat);
			}
		}
		channel.push_back(singleChannel);
	}

	std::cout << "NA" << ")" << std::endl;
	cycle_channel_Vec.push_back(channel);
 }


int ImageCycle::chastity(const std::vector<int>& values) 
{
	int largest, second;
	for (int i = 0; i< values.size() ; i ++) {
      /* If the current array element is greater than largest*/
      if (values[i] > largest) {
         second = largest;
         largest = values[i];
      }
      /* If current array element is less than largest but greater
       * then second largest ("second" variable) then copy the
       * element to "second"
       */
      else if (values[i] > second && values[i] != largest) {
         second = values[i];
      }
   }
   double chastity = (double)largest/(double)(largest+second);
   return(chastity);
} 




void ImageCycle::histogram(const cv::Mat& _i, const cv::Mat& mask, double* h, double* cdf, std::vector<double> &entropy) {
	cv::Mat _t = _i.reshape(1,1);

	if (!mask.empty()) {
		cv::Mat _tm;
		mask.copyTo(_tm);
		_tm = _tm.reshape(1, 1);
        //to handle
		for(int p=0; p<_t.cols; p++){
			ushort m = _tm.at<ushort>(0, p);
			if(m > 0){
				ushort c = _t.at<ushort>(0, p);
				h[c] += 1.0;
			}
		}
	}else{
		for(int p=0; p<_t.cols; p++){
			ushort c = _t.at<ushort>(0, p);
			h[c] += 1.0;
		}
	}

	//normalize histogram to have maximum of 1.0
	cv::Mat _tmp(1, 65536, CV_64FC1, h);
	double minVal, maxVal;
	cv::minMaxLoc(_tmp, &minVal, &maxVal);
	_tmp = _tmp / maxVal;

	cv::Mat logP;
	cv::Mat eTmp;
	_tmp.copyTo(eTmp);
	eTmp += 1e-4;
	cv::log(eTmp, logP);
	entropy.push_back( (-1*sum(eTmp.mul(logP)).val[0]) );

	// Calculate the CDF
	cdf[0] = h[0];
	for(int j=1; j<65536; j++)
		cdf[j] = cdf[j-1]+h[j];

	//normalize CDF to have maximum of 1.0
	_tmp.data = (uchar*)cdf;
	cv::minMaxLoc(_tmp, &minVal, &maxVal);
	_tmp = _tmp / maxVal;

}

void ImageCycle::hist_match(int refChannel){

	double const HISTMATCH_EPSILON = 0.000001;
	for(int i = 0; i < 4; i++){
	cv::Mat Mv(1, 65536, CV_16U);
	ushort* M = Mv.ptr<ushort>();		

		ushort last = 0;	
	for(int j=0; j < cdfs[i].cols; j++){
		double F1j = cdfs[i].at<double>(0, j);
		for(ushort k = last; k < cdfs[refChannel].cols; k++) {
			double F2k = cdfs[refChannel].at<double>(0, k);
			if (F2k > F1j - HISTMATCH_EPSILON) {
				M[j] = k;
				last = k;
				break;
			}
		}
	}

	

		cv::Mat lut(1, 65536, CV_16U, M);
		//std::cout << "lut" << i << "<- c(";
		//for(int u = 0; u < 65536; u++)
			//std::cout << lut.at<ushort>(0, u) << ", ";
			//std::cout << "NA" << ")" << std::endl;
		for (int l = i; l < multipageImage.size(); l = l + 4) {
			for(int p = 0; p < multipageImage[l].rows; p++){
				for(int q = 0; q < multipageImage[l].cols; q++){
					multipageImage[l].at<ushort>(p, q) = lut.at<ushort>(0, multipageImage[l].at<ushort>(p, q) );
				}
			}
	  	}
	 }

}

std::vector<cv::Mat> reorder(const std::vector<cv::Mat>& vA, const std::vector<size_t>& vOrder)  
{   
    std::vector<cv::Mat> vCopy(vA.size()); 
    for(int i = 0; i < vOrder.size(); ++i)  
        vCopy[i] = vA[ vOrder[i] ];  
    return vCopy;
} 

void ImageCycle::morph3d(bool isDilation)
{

	std::vector<size_t> order;
	for(int i = 0; i <  multipageMask.size()/4; i++ ){
        for(int j = 0; j < 4; j++){
                order.push_back(i+ (multipageMask.size()/4)*j);
        }
    }

	std::vector<cv::Mat> inData;
	for (int i = 0; i < 4; i++) {
		for (int l = i; l < multipageMask.size(); l = l + 4) {
    		backgroundSignal.push_back(multipageMask[l].clone());
    		inData.push_back(multipageMask[l].clone());
  		}
  	}

  	const int dims[] = {multipageMask[0].rows, multipageMask[0].cols, (int)multipageMask.size()/4, 4};
  	const int window[] = {7, 7, 7};

    //const int dims       = in.dims();
    //const int window     = mask.dims();
    const int R0     = window[0]/2;
    const int R1     = window[1]/2;
    const int R2     = window[2]/2;
    //const int istrides   = in.strides();
    //const int fstrides   = mask.strides();
    const int bCount = dims[3];

    //create copy of array to serve as output
    /*Array<T> out         = createEmptyArray<T>(dims);
    const dim4 ostrides   = out.strides();

    T* outData            = out.get();
    const T*   inData     = in.get();
    const T*   filter     = mask.get();*/
    for(int channelId=0; channelId<bCount; ++channelId) {
        // channels are handled by outer most loop
        for(int k=0; k<dims[2]; ++k) {
            // k steps along 3rd dimension
            for(int j=0; j<dims[1]; ++j) {
                // j steps along 2nd dimension
                for(int i=0; i<dims[0]; ++i) {
                    // i steps along 1st dimension
                    ushort filterResult = inData[k+dims[2]*channelId].at<ushort>(i, j);

                    // wk, wj,wi steps along 2nd & 1st dimensions of filter window respectively
                    for(int wk=0; wk<window[2]; wk++) {
                        for(int wj=0; wj<window[1]; wj++) {
                            for(int wi=0; wi<window[0]; wi++) {

                                int offk = k+wk-R2;
                                int offj = j+wj-R1;
                                int offi = i+wi-R0;

                                ushort maskValue = 1;//;multipageMask[wk].at<ushort>(wi, wj);// filter[ getIdx(fstrides, wi, wj, wk) ];

                                if ((maskValue > 0) && offi>=0 && offj>=0 && offk>=0 &&
                                        offi<dims[0] && offj<dims[1] && offk<dims[2]) {

                                    ushort inValue   = inData[offk+dims[2]*channelId].at<ushort>(offi, offj); // inData[ getIdx(istrides, offi, offj, offk) ];

                                    if (isDilation)
                                        filterResult = std::max(filterResult, inValue);
                                    else
                                        filterResult = std::min(filterResult, inValue);
                                }

                            } // window 1st dimension loop ends here
                        }  // window 1st dimension loop ends here
                    }// filter window loop ends here
                    backgroundSignal[k+dims[2]*channelId].at<ushort>(i, j) = filterResult;
                    //outData[ getIdx(ostrides, i, j, k) ] = filterResult;
                } //1st dimension loop ends here
            } // 2nd dimension loop ends here
        } // 3rd dimension loop ends here
        // next iteration will be next channel if any
        //outData += ostrides[3];
        //inData  += istrides[3];
    }
    


for(int channelId=0; channelId<bCount; ++channelId) {
    for(int k=0; k<dims[2]; ++k) {
        // k steps along 3rd dimension
        for(int j=0; j<dims[1]; ++j) {
            // j steps along 2nd dimension
            for(int i=0; i<dims[0]; ++i) {
            	if(inData[k+dims[2]*channelId].at<ushort>(i, j) ==  backgroundSignal[k+dims[2]*channelId].at<ushort>(i, j) )
            		backgroundSignal[k+dims[2]*channelId].at<ushort>(i, j) = backgroundSignal[k+dims[2]*channelId].at<ushort>(i, j) - inData[k+dims[2]*channelId].at<ushort>(i, j);
            }
        }
    }
   }

    backgroundSignal = reorder(backgroundSignal, order);

}
