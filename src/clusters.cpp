#include "clusters.hpp"

Clusters::Clusters(std::vector<std::vector<int>> ampliconSignal,
            std::vector<std::vector<int>> ampliconNoise,
             int correction_cycle,                            ///< Build crosstalk correction based on this cycle
             int iteration_threshold,                   ///< Keep trying at most this number of times.
             double slope_threshold ,                   ///< Iteration stops when slopes are less than this
             double crosstalk_lowerpercentile,
             double crosstalk_upperpercentile) : ampliconSignal(ampliconSignal),
             ampliconNoise(ampliconNoise),
             correction_cycle(correction_cycle),
             iteration_threshold(iteration_threshold),
             slope_threshold(slope_threshold),
             crosstalk_lowerpercentile(crosstalk_lowerpercentile),
             crosstalk_upperpercentile(crosstalk_upperpercentile) {


};

Clusters::~Clusters()
{

}


double Clusters::get_percentile(std::vector<double> c, double prob) {


  double index = 1 + ((double)c.size() - 1) * prob;
  int lo = std::floor(index) - 1;
  int hi = std::ceil(index) - 1;

  sort(c.begin(),c.end()); //change this to partial sort using hi and lo when you have time

  double qs = c[lo];

  int h = index - lo;
  qs = (1 - h) * qs + h * c[hi];

  return(qs);
  //return c[static_cast<int>((static_cast<double>(prob))*static_cast<double>(c.size()))];
}


void Clusters::make_bins(std::vector<double> &x,             /// X values
                         std::vector<double> &y,             /// Y values
                         std::vector<double> &xo,
                         std::vector<double> &yo) { /// Bins must contain less than this many items
  // 1. Bin 65th -> 95th Percentile along each axis

  // 1.1 Find 65th and 95th percentile locations

  // 1.1.3.1 Find position in count vector where we go over lower percentile limit
  // 1.1.3.2 Find position in count vector where we go over upper percentile limit
  
  double percentile_lower_limit_position = get_percentile(x,crosstalk_lowerpercentile); // 65
  double percentile_upper_limit_position = get_percentile(x,crosstalk_upperpercentile); // 95

  //std::cout << "Lower percentile Limit: " << percentile_lower_limit_position << std::endl;
  //std::cout << "Upper percentile Limit: " << percentile_upper_limit_position << std::endl;

  // 2. Take lowest value in every bin (num_bins bins)

  std::vector<double> bins;
  std::vector<double> bins_x;

  //
  int current_num_bins = 50;
  std::vector<bool>  bins_first(current_num_bins,true);
  std::vector<int>  num_in_bins(current_num_bins,0);

  bins.clear();
  bins_x.clear();
  for(int n=0;n<=current_num_bins;n++) bins.push_back(0);
  for(int n=0;n<=current_num_bins;n++) bins_x.push_back(0);
  
  double del = percentile_upper_limit_position - percentile_lower_limit_position;  
  double binsize = (del/(current_num_bins - 1));

  //double binsize = std::ceil( (percentile_upper_limit_position-percentile_lower_limit_position)/current_num_bins );
  for(unsigned int n = 0;n < x.size();n++) {
      // this is very inefficient...
      for(int current_bin_number = 0; current_bin_number < current_num_bins; current_bin_number++) {
        
          double current_bin_start = std::round( percentile_lower_limit_position+(binsize*current_bin_number) ) ;
          double current_bin_end   = std::round( percentile_lower_limit_position+(binsize*current_bin_number)+binsize);

          //std::cout << "bin start: " << current_bin_start << std::endl;
          //std::cout << "bin   end: " << current_bin_end << std::endl;

          

          if((x[n] > current_bin_start) && (x[n] <= current_bin_end)) {
            if(bins_first[current_bin_number] || (bins[current_bin_number] > y[n])) { bins  [current_bin_number] = static_cast<double>(y[n]);
                                                                                      bins_x[current_bin_number] = static_cast<double>(x[n]); }

                    
            bins_first[current_bin_number] = false;

            num_in_bins[current_bin_number]++;
          }
        
      }
    }/*

    std::vector<double> bins_real;
    std::vector<double> bins_x_real;

    bins_real.clear();
    bins_x_real.clear();
    int real_bins=0;
    for(size_t n=0;n<bins_first.size();n++) {
      if((bins_first[n] == false) && (((bins[n] > 1) || (bins[n] < -1)) && ((bins_x[n] > 1) || (bins_x[n] < -1)))) {

        bins_real.push_back(bins[n]);
        bins_x_real.push_back(bins_x[n]);
        real_bins++;
      } else {
        //std::cout << "Removing empty or zero bin: " << n << std::endl;
      }
    }
    bins   = bins_real;
    bins_x = bins_x_real;  
    */

  xo = bins_x;
  yo = bins;
  // plotxy(yo,xo,"Crosstalk");
}

void Clusters::slope(double &slope, const std::vector<double>& x, const std::vector<double>& y, int origin = 0) {
    const auto n    = x.size();
    const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    double b;
    switch(origin) {
      case 0:
        b = s_xy/s_xx;
        break;
      case 1:
        b    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
        break;
      default:
        b = s_xy/s_xx;
    }
    slope = b;
    //return b;
}


void Clusters::apply_correction_values() {
  for(size_t n=0;n<ch01_noise.size();n++) {
    double ch1_val = ch01_noise[n];
    double ch2_val = ch02_noise[n];
    double ch3_val = ch03_noise[n];
    double ch4_val = ch04_noise[n];

    apply_correction(ch1_val,ch2_val,ch3_val,ch4_val);

    ch01_noise[n] = ch1_val;
    ch02_noise[n] = ch2_val;
    ch03_noise[n] = ch3_val;
    ch04_noise[n] = ch4_val;
  }


  for(size_t n=0;n<ch01_values.size();n++) {
    double ch1_val = ch01_values[n];
    double ch2_val = ch02_values[n];
    double ch3_val = ch03_values[n];
    double ch4_val = ch04_values[n];

    apply_correction(ch1_val,ch2_val,ch3_val,ch4_val);

    ch01_values[n] = ch1_val;
    ch02_values[n] = ch2_val;
    ch03_values[n] = ch3_val;
    ch04_values[n] = ch4_val;
  }

}    


void Clusters::apply_correction(double &ch1_val, double &ch2_val, double &ch3_val, double &ch4_val) {

  //Correct channel01
  double old_ch1 = ch1_val;
  double old_ch3 = ch3_val;
  double old_ch4 = ch4_val;
  double old_ch2 = ch2_val;


  double new_ch1 = old_ch1 -ch02ch01_m*old_ch2 -ch03ch01_m*old_ch3 -ch04ch01_m*old_ch4;
  double new_ch3 = old_ch3 -ch01ch03_m*old_ch1 -ch02ch03_m*old_ch2 -ch04ch03_m*old_ch4;
  double new_ch4 = old_ch4 -ch01ch04_m*old_ch1 -ch02ch04_m*old_ch2 -ch03ch04_m*old_ch3; 
  double new_ch2 = old_ch2 -ch01ch02_m*old_ch1 -ch03ch02_m*old_ch3 -ch04ch02_m*old_ch4;


  ch1_val = new_ch1;
  ch2_val = new_ch2;
  ch3_val = new_ch3;
  ch4_val = new_ch4;


}

void Clusters::replace_negative_values() {
  for(size_t n=0;n<ch01_noise.size();n++) {
    double ch1_val = ch01_noise[n];
    double ch2_val = ch02_noise[n];
    double ch3_val = ch03_noise[n];
    double ch4_val = ch04_noise[n];

    if(ch1_val < 0){
        ch1_val = 0;
    }
    if(ch2_val < 0){
        ch2_val = 0;
    }
    if(ch3_val < 0){
        ch3_val = 0;
    }
    if(ch4_val < 0){
        ch4_val = 0;
    }

    ch01_noise[n] = ch1_val;
    ch02_noise[n] = ch2_val;
    ch03_noise[n] = ch3_val;
    ch04_noise[n] = ch4_val;
  }


  for(size_t n=0;n<ch01_values.size();n++) {
    double ch1_val = ch01_values[n];
    double ch2_val = ch02_values[n];
    double ch3_val = ch03_values[n];
    double ch4_val = ch04_values[n];

    if(ch1_val < 0){
        ch1_val = 0;
    }
    if(ch2_val < 0){
        ch2_val = 0;
    }
    if(ch3_val < 0){
        ch3_val = 0;
    }
    if(ch4_val < 0){
        ch4_val = 0;
    }

    ch01_values[n] = ch1_val;
    ch02_values[n] = ch2_val;
    ch03_values[n] = ch3_val;
    ch04_values[n] = ch4_val;
  }

} 

void Clusters::subtract_median_for_each_channel() {
    //get median and subtract
    double ch01_median = get_percentile(ch01_values, 0.5);
    double ch02_median = get_percentile(ch02_values, 0.5);
    double ch03_median = get_percentile(ch03_values, 0.5);
    double ch04_median = get_percentile(ch04_values, 0.5);

    for(size_t n=0;n<ch01_values.size();n++) {

        double ch1_val = ch01_values[n];
        double ch2_val = ch02_values[n];
        double ch3_val = ch03_values[n];
        double ch4_val = ch04_values[n];

        double new_ch1 = ch1_val - ch01_median;
        double new_ch2 = ch2_val - ch02_median;
        double new_ch3 = ch3_val - ch03_median;
        double new_ch4 = ch4_val - ch04_median;

        ch01_values[n] = new_ch1;
        ch02_values[n] = new_ch2;
        ch03_values[n] = new_ch3;
        ch04_values[n] = new_ch4;

    }

}


int Clusters::max2nd(std::vector<int> arr) 
{ 
  /*
    int i, first, second; 

  
    first = second = INT_MIN; 
    for (i = 0; i < arr.size() ; i ++) 
    { 

        if (arr[i] > first) 
        { 
            second = first; 
            first = arr[i]; 
        } 
  
        else if (arr[i] > second && arr[i] != first) 
            second = arr[i]; 
    } */
  sort(arr.begin(),arr.end());

  return(arr[2]);
}

void Clusters::recompute_chastity() {
  chastityValue.clear();  

  for(size_t n=0;n<ch01_values.size();n++) {

    std::vector<int> singleAmplicon;
    singleAmplicon.push_back((int)ch01_values[n]);
    singleAmplicon.push_back((int)ch02_values[n]);
    singleAmplicon.push_back((int)ch03_values[n]);
    singleAmplicon.push_back((int)ch04_values[n]);

    chastityValue.push_back( chastity(singleAmplicon) );

  }

}   

double Clusters::chastity(std::vector<int> arr)
{
    int max = *max_element(std::begin(arr), std::end(arr));
    int secondMax = max2nd(arr);
    double chastityValue = (double)max/((double)max + (double)secondMax);
    return(chastityValue);
}

double Clusters::l2_norm(std::vector<double> const& u) {
    double accum = 0.;
    for (double x : u) {
        accum += x * x;
    }
    return std::sqrt(accum);
}



void Clusters::basecall(double chastityColoc, double chastityThresh)
{
  std::vector<int> base;
  //std::vector<double> quality;
  std::vector<int> ampliconID;
  std::vector<double> angle;
  std::vector<double> intensity;

  std::vector<std::vector<double>> polar;

  double pi = 2 * std::acos(0.0);

  int direction[] = {-1,1};
  double angleAdjust[] = {1.0/4.0*pi, 5.0/4.0*pi, 3.0/4.0*pi, 7.0/4.0*pi};

  std::vector<double> coord;
  std::vector<double> singleAmplicon;
  std::vector<double> reverseRead;
  std::vector<double> forwardRead;


  for(size_t n=0;n<ch01_values.size();n++) {
    if(chastityValue[n] >= chastityThresh) {
      singleAmplicon.clear();
      singleAmplicon.push_back(ch01_values[n]);
      singleAmplicon.push_back(ch02_values[n]);
      singleAmplicon.push_back(ch03_values[n]);
      singleAmplicon.push_back(ch04_values[n]);

      int maxChannel = std::distance(singleAmplicon.begin(), std::max_element(singleAmplicon.begin(), singleAmplicon.end()));
      base.push_back(maxChannel);
      ampliconID.push_back(n);

      double Cx = singleAmplicon.at(maxChannel) / l2_norm(singleAmplicon);

      //quality.push_back(q);

      int randomNumber = rand() % 2;

      double theta = (std::acos(1.0 - Cx) - pi/2) * direction[randomNumber] + angleAdjust[maxChannel];

      angle.push_back( theta );
      double value = singleAmplicon.at(maxChannel);
      intensity.push_back( value );

      coord.clear();
      coord.push_back(value*std::cos(theta));
      coord.push_back(value*std::sin(theta));
      polar.push_back(coord);

      std::cout << ch01_values[n] << ", " << ch02_values[n] << ", " << ch03_values[n] << ", " << ch04_values[n] << ", " << value <<  ", " << theta  << ", " << chastityValue[n] << ", " << maxChannel << std::endl;
      //std::cout << singleAmplicon.at(0) << ", " << singleAmplicon.at(1) << ", " <<singleAmplicon.at(2) << ", " << singleAmplicon.at(3) << " | " << l2_norm(singleAmplicon) << " | " << maxChannel << std::endl;

    }
    if(chastityValue[n] < chastityColoc) {
      reverseRead.clear();

      if(ch01_values[n] > ch04_values[n]){
        double x = ch01_values[n];
        double y = ch03_values[n];
        double x_new  = x * (std::cos(pi/4.0))  - y * (std::sin(pi/4.0));
        double y_new = x * (std::sin(pi/4.0))  + y * (std::cos(pi/4.0)); 
    

        double theta = atan2(y_new, x_new) - atan2(65536, 0) - pi/2;
        angle.push_back( theta );
        intensity.push_back( y_new );

        coord.clear();
        coord.push_back(y_new*std::cos(theta));
        coord.push_back(y_new*std::sin(theta));
        polar.push_back(coord);

        std::cout << ch01_values[n] << ", " << ch02_values[n] << ", " << ch03_values[n] << ", " << ch04_values[n] << ", " << y_new <<  ", " << theta  << ", " << chastityValue[n] << ", " << -2 << std::endl;

      }else{
        double x = ch04_values[n];
        double y = ch02_values[n];
        double x_new  = x * (std::cos(pi/4.0))  - y * (std::sin(pi/4.0));
        double y_new = x * (std::sin(pi/4.0))  + y * (std::cos(pi/4.0)); 
    

        double theta = atan2(y_new, x_new) - atan2(65536, 0) + pi/2;
        angle.push_back( theta );
        intensity.push_back( y_new );

        coord.clear();
        coord.push_back(y_new*std::cos(theta));
        coord.push_back(y_new*std::sin(theta));
        polar.push_back(coord);

        std::cout << ch01_values[n] << ", " << ch02_values[n] << ", " << ch03_values[n] << ", " << ch04_values[n] << ", " << y_new <<  ", " << theta  << ", " << chastityValue[n] << ", " << -1 << std::endl;

      }




    }

  }

}


bool Clusters::initialise(const std::vector<std::vector<int>>& ampliconSignal, const std::vector<std::vector<int>>& ampliconNoise, double chastityThresh) {

  ch01_values.clear();
  ch02_values.clear();
  ch03_values.clear();
  ch04_values.clear();

  ch01_noise.clear();
  ch02_noise.clear();
  ch03_noise.clear();
  ch04_noise.clear();

  chastityValue.clear();

  for(int i = 0; i < ampliconSignal.size(); i++){
    chastityValue.push_back( chastity(ampliconSignal[i]) );
    if( chastityValue[i] > chastityThresh){
        ch01_values.push_back((double)ampliconSignal[i].at(0));
        ch02_values.push_back((double)ampliconSignal[i].at(1));
        ch03_values.push_back((double)ampliconSignal[i].at(2));
        ch04_values.push_back((double)ampliconSignal[i].at(3));

        ch01_noise.push_back((double)ampliconNoise[i].at(0));
        ch02_noise.push_back((double)ampliconNoise[i].at(1));
        ch03_noise.push_back((double)ampliconNoise[i].at(2));
        ch04_noise.push_back((double)ampliconNoise[i].at(3));
      }

    
  }

  return true;
}
