#include "registration.hpp"
#include<iostream>

Registration::Registration(){
	
}

Registration::~Registration()
{

}

void Registration::loadCSV(std::string segmentedCSV){
	segmented.load(segmentedCSV);
	arma::umat columns;
	columns << 2 << 3 << 4 << arma::endr;
	int k = 1;
	cycle_channel_Vec.clear();
	for(int i = 1; i < segmented.n_rows; i++){
		if(segmented(i, 1) != segmented(i-1, 1)) {  
			std::vector<arma::mat> channel;
			for(int ch = 0; ch < 4; ch++)
				channel.push_back( segmented.submat( find( segmented.col(0) == k && segmented.col(1) == ch ), columns ) ); 
			cycle_channel_Vec.push_back(channel);
			k++;
		};
	}

}

void Registration::loadCYCLE(const std::vector<std::vector<arma::mat>>& cycle){
	cycle_channel_Vec.clear();
	cycle_channel_Vec = cycle;
}

void Registration::baseCallG(int cycle)
{
	PointReg regi;
	std::cout << "ch01_matched_ch03" << std::endl;
	if( cycle_channel_Vec.at(cycle).at(2).n_rows   > cycle_channel_Vec.at(cycle).at(0).n_rows  ){
		regi.run(cycle_channel_Vec.at(cycle).at(2), cycle_channel_Vec.at(cycle).at(0) );
	}else{
		regi.run(cycle_channel_Vec.at(cycle).at(0), cycle_channel_Vec.at(cycle).at(2) );
	}

    PointReg regi2;
    std::cout << "ch02_matched_ch04" << std::endl;
    if( cycle_channel_Vec.at(cycle).at(3).n_rows   > cycle_channel_Vec.at(cycle).at(1).n_rows  ){
		regi2.run(cycle_channel_Vec.at(cycle).at(3), cycle_channel_Vec.at(cycle).at(1) );
	}else{
		regi2.run(cycle_channel_Vec.at(cycle).at(1), cycle_channel_Vec.at(cycle).at(3) );
	}
	

}

void Registration::matchFwdRev(int cycle)
{
	PointReg regi;
	    std::cout << "ch01ch02" << std::endl;

    regi.run(cycle_channel_Vec.at(cycle).at(0), cycle_channel_Vec.at(cycle).at(1) );
	
	PointReg regi1;
	    std::cout << "ch01ch04" << std::endl;

    regi1.run(cycle_channel_Vec.at(cycle).at(0), cycle_channel_Vec.at(cycle).at(3) );
	
	PointReg regi2;	   
	    std::cout << "ch03ch02" << std::endl;
 
    regi2.run(cycle_channel_Vec.at(cycle).at(2), cycle_channel_Vec.at(cycle).at(1) );
	
	PointReg regi3;
	    std::cout << "ch03ch04" << std::endl;

    regi3.run(cycle_channel_Vec.at(cycle).at(2), cycle_channel_Vec.at(cycle).at(3) );

    
}