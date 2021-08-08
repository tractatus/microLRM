#include "pointreg.hpp"

PointReg::~PointReg()
{

}


void PointReg::run(arma::mat FIXED, arma::mat MOVING){

    
	std::unique_ptr<cpdArma::Registration> reg;
	//std::cout << "started" << std::endl;

    reg = std::unique_ptr<cpdArma::Registration>(new cpdArma::NonrigidLowrank(
                    m_tol,
                    m_max_it,
                    m_outliers,
                    m_usefgt,
                    m_epsilon,
                    m_beta,
                   	m_lambda,
                    m_numeig
                ));

    cpdArma::Registration::ResultPtr result = reg->run(FIXED, MOVING);
        	//std::cout << "finished" << std::endl;

    cpd::Matrix moving = cast_arma_to_eigen(result->Y);
    cpd::Matrix fixed = cast_arma_to_eigen(FIXED);
            	//std::cout << "casted" << std::endl;
    corr(fixed, moving);
    //cpd::Vector correspondence = resolveAmbigious(fixed, moving);
            	//std::cout << "resolved" << std::endl;

    //onlyLeaveShortestDist(fixed, moving, correspondence);
            	//std::cout << "shortest" << std::endl;
    //std::cout << correspondence.size() << std::endl;


    //return(correspondence);
}

std::tuple<std::vector<std::size_t>, std::vector<std::size_t>> PointReg::findDuplicateIndices(std::vector<int> const & v)
{
    std::vector<std::size_t> firstoccurance;
    std::vector<std::size_t> duplicates;
    std::map<int, std::pair<int, std::size_t>> counts; // pair<amount, first position seened>

    for (std::size_t i = 0 ; i < v.size() ; ++i)
    {
        std::size_t const amount = ++counts[v[i]].first;
        if (amount == 1)
        {
            counts[v[i]].second = i;
            continue;
        }
        else if (amount == 2) {
            firstoccurance.push_back(counts[v[i]].second);//comment this out if only second element is desired
            duplicates.push_back(i);
        }
        if (amount > 2) {
            duplicates.push_back(i);
        }

    }
    return {firstoccurance, duplicates};
}

bool sortbysec(const std::tuple<int, int, double>& a,  
               const std::tuple<int, int, double>& b) 
{ 
    return (std::get<1>(a) < std::get<1>(b)); 
} 

void PointReg::corr(const cpd::Matrix & fixed, const cpd::Matrix & moving)
{

    double ksig = -2.0 * 0.01;
    cpd::Vector p = cpd::Vector::Zero(moving.rows());
    cpd::Vector p1 = cpd::Vector::Zero(moving.rows());
    cpd::Vector p1_max = cpd::Vector::Zero(moving.rows());
    cpd::Vector pt1 = cpd::Vector::Zero(fixed.rows());

    cpd::Vector correspondence = - cpd::Vector::Ones(moving.rows());
    for (cpd::Matrix::Index i = 0; i < fixed.rows(); ++i) {
        double sp = 0;
        for (cpd::Matrix::Index j = 0; j < moving.rows(); ++j) {
            double razn = (fixed.row(i) - moving.row(j)).array().pow(2).sum();
            razn = std::sqrt(razn);
            p(j) = std::exp(razn / ksig);
            sp += p(j);
        }
        for (cpd::Matrix::Index j = 0; j < moving.rows(); ++j) {
            p1(j) += p(j) / sp;
            if (p(j) / sp > p1_max(j)) {
                correspondence(j) = i;
                p1_max(j) = p(j) / sp;
            }
        }
    }

    std::vector<int> vec;
    for(int i = 0; i < correspondence.size(); i++){
        vec.push_back( (int)correspondence(i) );
    }
 

    auto [matched_to_already_matched, duplicated] = findDuplicateIndices(vec);

    std::vector<int> diff;

    std::vector<int> orig(correspondence.size()) ;
    std::iota (std::begin(orig), std::end(orig), 0);
    sort(vec.begin(),vec.end());

    std::set_difference(orig.begin(), orig.end(), vec.begin(), vec.end(), 
                        std::inserter(diff, diff.begin()));

    for(int i = 0; i < duplicated.size(); i++){
    	int k = 0;
    	for(int j = 0; j < diff.size(); j++){
    		if(diff.at(j) == correspondence(duplicated[i]) ){
    			k++;
    		}
    	}
    	if(k == 0){
    		diff.push_back( correspondence(duplicated[i]) );
    	}
    }


    sort(diff.begin(),diff.end());
     //for(int i = 0; i < diff.size(); i++){

    	//std::cout << diff.at(i) << std::endl;
    //}   


    cpd::Matrix fixedSubset = fixed(diff, Eigen::all);
    cpd::Matrix movingSubset = moving(matched_to_already_matched, Eigen::all);

    cpd::Vector p2 = cpd::Vector::Zero(movingSubset.rows());
    cpd::Vector p12 = cpd::Vector::Zero(movingSubset.rows());
    cpd::Vector p1_max2 = cpd::Vector::Zero(movingSubset.rows());
    cpd::Vector pt12 = cpd::Vector::Zero(fixedSubset.rows());

    cpd::Vector correspondence2 = cpd::Vector::Zero(movingSubset.rows());

    for (cpd::Matrix::Index i = 0; i < fixedSubset.rows(); ++i) {
        double sp = 0;
        for (cpd::Matrix::Index j = 0; j < movingSubset.rows(); ++j) {
            double razn = (fixedSubset.row(i) - movingSubset.row(j)).array().pow(2).sum();
            razn = std::sqrt(razn);
            p2(j) = std::exp(razn / ksig);
            sp += p2(j);
        }
        for (cpd::Matrix::Index j = 0; j < movingSubset.rows(); ++j) {
            p12(j) += p2(j) / sp;
            if (p2(j) / sp > p1_max2(j)) {
                correspondence2(j) = i;
                p1_max2(j) = p2(j) / sp;
            }
        }
    }

    for(int i = 0; i < matched_to_already_matched.size(); i++){

        correspondence(matched_to_already_matched[i]) = diff[correspondence2[i]]; 

    }

    

    //shortest dist
        //tuple with (i) index in moving, (ii) index in fixed and (iii) Euclidean distance
    std::vector<std::tuple<int, int, double> > vec2; 
    double dist;
    for(int i = 0; i < correspondence.size(); i++){ //for(int i = 0; i < correspondance.size(); i++){
    	if(correspondence(i) > fixed.rows()){
    		//this should not happen so probably some memory issues.
    		//std::cout << "TOO LONG: "<< i << " corr: " << correspondance[i] << std::endl;
    		correspondence(i) = -1;
    	}
        		     
        if(correspondence(i) != -1){
            dist = std::sqrt( (fixed.row(correspondence(i)) - moving.row(i)).array().pow(2).sum() );
        } else {
            dist = 0;
        }
        
        vec2.emplace_back(i, (int)correspondence(i), dist);
    }

    std::sort(vec2.begin(), vec2.end(), sortbysec);


    std::map<int, std::pair<int, std::size_t>> counts; // pair<amount, first position seened>


    for (std::size_t i = 0 ; i < vec2.size() ; ++i)
    {
        std::size_t const amount = ++counts[std::get<1>(vec2.at(i))].first;
        if (amount == 1)
        {
            counts[std::get<1>(vec2.at(i))].second = i;
            continue;
        }
        else if (amount >= 2) {
            int k = counts[std::get<1>(vec2.at(i))].second;
            if( std::get<2>(vec2.at(k)) > std::get<2>(vec2.at(i)) ) {
                std::get<1>(vec2.at(k)) = -1;
                counts[std::get<1>(vec2.at(i))].second = i;
            }else{
                std::get<1>(vec2.at(i)) = -1;
            }
        }
    }
    
    std::sort(vec2.begin(), vec2.end());
    for(int i = 0; i < correspondence.size(); i++){
        correspondence(i) = std::get<1>(vec2.at(i));
    }

    std::cout << correspondence << std::endl;

    //return(correspondence);
    
    
}

/*

cpd::Vector PointReg::computeCorrIndex(const cpd::Matrix & fixed, const cpd::Matrix & moving)
{
    double ksig = -2.0 * 0.01;
    cpd::Vector p = cpd::Vector::Zero(moving.rows());
    cpd::Vector p1 = cpd::Vector::Zero(moving.rows());
    cpd::Vector p1_max = cpd::Vector::Zero(moving.rows());
    cpd::Vector pt1 = cpd::Vector::Zero(fixed.rows());
        size_t cols = fixed.cols();

    cpd::Matrix px = cpd::Matrix::Zero(moving.rows(), cols);
    cpd::Vector correspondence = - cpd::Vector::Ones(moving.rows());
    double l = 0.0;
    for (cpd::Matrix::Index i = 0; i < fixed.rows(); ++i) {
        double sp = 0;
        for (cpd::Matrix::Index j = 0; j < moving.rows(); ++j) {
            double razn = (fixed.row(i) - moving.row(j)).array().pow(2).sum();
            razn = std::sqrt(razn);
            p(j) = std::exp(razn / ksig);
            sp += p(j);
        }
        for (cpd::Matrix::Index j = 0; j < moving.rows(); ++j) {
            p1(j) += p(j) / sp;
            px.row(j) += fixed.row(i) * p(j) / sp;
            if (p(j) / sp > p1_max(j)) {
                correspondence(j) = i;
                p1_max(j) = p(j) / sp;
            }
        }
    }
    return(correspondence);
}

cpd::Vector PointReg::resolveAmbigious(const cpd::Matrix & fixed, const cpd::Matrix & moving)
{
    cpd::Vector oldCorr = computeCorrIndex(fixed, moving);

    std::vector<int> vec;
    for(int i = 0; i < oldCorr.size(); i++){
        vec.push_back( (int)oldCorr[i] );
    }
 

    std::vector<std::size_t> matched_to_already_matched = findDuplicateIndices(vec);

    std::vector<int> diff;

    std::vector<int> orig(oldCorr.size()) ;
    std::iota (std::begin(orig), std::end(orig), 0);
    sort(vec.begin(),vec.end());

    std::set_difference(orig.begin(), orig.end(), vec.begin(), vec.end(), 
                        std::inserter(diff, diff.begin()));

    for(int i = 0; i < matched_to_already_matched.size(); i+=2)
        diff.push_back( oldCorr[matched_to_already_matched[i]]);

    sort(diff.begin(),diff.end());

    cpd::Matrix fixedSubset = fixed(diff, Eigen::all);
    cpd::Matrix movingSubset = moving(matched_to_already_matched, Eigen::all);

    double ksig = -2.0 * 0.01;
    cpd::Vector p = cpd::Vector::Zero(movingSubset.rows());
    cpd::Vector p1 = cpd::Vector::Zero(movingSubset.rows());
    cpd::Vector p1_max = cpd::Vector::Zero(movingSubset.rows());
    cpd::Vector pt1 = cpd::Vector::Zero(fixedSubset.rows());
    size_t cols = fixedSubset.cols();

    cpd::Matrix px = cpd::Matrix::Zero(movingSubset.rows(), cols);
    cpd::Vector correspondence = cpd::Vector::Zero(movingSubset.rows());
    double l = 0.0;
    for (cpd::Matrix::Index i = 0; i < fixedSubset.rows(); ++i) {
        double sp = 0;
        for (cpd::Matrix::Index j = 0; j < movingSubset.rows(); ++j) {
            double razn = (fixedSubset.row(i) - movingSubset.row(j)).array().pow(2).sum();
            razn = std::sqrt(razn);
            p(j) = std::exp(razn / ksig);
            sp += p(j);
        }
        for (cpd::Matrix::Index j = 0; j < movingSubset.rows(); ++j) {
            p1(j) += p(j) / sp;
            px.row(j) += fixedSubset.row(i) * p(j) / sp;
            if (p(j) / sp > p1_max(j)) {
                correspondence(j) = i;
                p1_max(j) = p(j) / sp;
            }
        }
    }

    for(int i = 0; i < matched_to_already_matched.size(); i++){
        for(int j = 0; j < oldCorr.size(); j++)
            oldCorr[matched_to_already_matched[i]] = diff[correspondence[i]];
            
    }


    return(oldCorr);

}

bool sortbysec(const std::tuple<int, int, double>& a,  
               const std::tuple<int, int, double>& b) 
{ 
    return (std::get<1>(a) < std::get<1>(b)); 
} 

void PointReg::onlyLeaveShortestDist(const cpd::Matrix & fixed, const cpd::Matrix & moving, Eigen::Ref<cpd::Vector> correspondance)
{   
    //tuple with (i) index in moving, (ii) index in fixed and (iii) Euclidean distance
    std::vector<std::tuple<int, int, double> > vec; 
    double dist;
    for(int i = 0; i < correspondance.size(); i++){ //for(int i = 0; i < correspondance.size(); i++){
    	if(correspondance[i]>fixed.rows()){
    		//this should not happen so probably some memory issues.
    		//std::cout << "TOO LONG: "<< i << " corr: " << correspondance[i] << std::endl;
    		correspondance[i] = -1;
    	}
        		     
        if(correspondance[i] != -1){
            dist = std::sqrt( (fixed.row(correspondance[i]) - moving.row(i)).array().pow(2).sum() );
        } else {
            dist = 0;
        }
        
        vec.emplace_back(i, (int)correspondance[i], dist);
    }

    std::sort(vec.begin(), vec.end(), sortbysec);


    std::map<int, std::pair<int, std::size_t>> counts; // pair<amount, first position seened>


    for (std::size_t i = 0 ; i < vec.size() ; ++i)
    {
        std::size_t const amount = ++counts[std::get<1>(vec.at(i))].first;
        if (amount == 1)
        {
            counts[std::get<1>(vec.at(i))].second = i;
            continue;
        }
        else if (amount >= 2) {
            int k = counts[std::get<1>(vec.at(i))].second;
            if( std::get<2>(vec.at(k)) > std::get<2>(vec.at(i)) ) {
                std::get<1>(vec.at(k)) = -1;
                counts[std::get<1>(vec.at(i))].second = i;
            }else{
                std::get<1>(vec.at(i)) = -1;
            }
        }
    }
    
    std::sort(vec.begin(), vec.end());
    for(int i = 0; i < correspondance.size(); i++){
        correspondance[i] = std::get<1>(vec.at(i));
    }

}

*/


Eigen::MatrixXd PointReg::cast_arma_to_eigen(arma::mat inputArma) {

  Eigen::MatrixXd eigenFromArma = Eigen::Map<Eigen::MatrixXd>(inputArma.memptr(),
                                                        inputArma.n_rows,
                                                        inputArma.n_cols);

  return eigenFromArma;
}

