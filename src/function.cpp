#include <Rcpp.h>
#include <limits>
#include <RcppEigen.h>
#include <math.h>
#include <sstream>

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

using namespace Rcpp;

// Calculates square of euclidean distance of all the rows of mat1 with rth row of mat2
NumericVector obsDistance(NumericMatrix& mat1, NumericMatrix& mat2, int r) {
  NumericVector out(mat1.rows());
  for(int i=0; i<mat1.rows(); i++) {
    out[i] = sum( (mat1(i,_)-mat2(r,_))*(mat1(i,_)-mat2(r,_)) );
  }
  return out;
}

// Calculates distance matrix for every row of coords with every other row
NumericMatrix calculateDistanceMatrix(NumericMatrix coords){
  NumericMatrix out(coords.rows(), coords.rows());
  for(int i=0; i<coords.rows(); i++) {
    out(i,_) = sqrt( obsDistance(coords, coords, i) );
  }
  return out;
}

// Returns a numeric vector with the same length as of grid units
// For gaussian radius_type, it returns coefficient vector like (0.95, 0.86, 0.23, 0.01)
// For letremy radius_type, it returns coefficient vector like (1,1,0,0,0)
NumericVector selectNei_f(int the_neuron, int iteration, double radius, NumericMatrix& distMatrix, std::string radius_type){
  NumericVector out( distMatrix.rows() ) ;
  if(radius_type.compare("letremy")==0) {
    for(int i=0; i<distMatrix.rows(); i++) {
      if(distMatrix(the_neuron,i)<=radius) 
        out[i] = 1.0;
      else 
        out[i] = 0.0;
    }
  } else {
    int n = distMatrix.rows();
    int len = (n*n-n)/2;
    int upperTri[len];
    for(int row=0, i=0; row<n; row++) {
      for(int col=row+1; col<n; col++,i++) {
        upperTri[i] = distMatrix(row,col)*distMatrix(row,col);
      }
    }
    std::sort( upperTri,upperTri+len );
    int sigma = upperTri[(len-1)/2];
    for(int i=0; i<n; i++) {
      out[i] = std::exp( -distMatrix(the_neuron,i)*distMatrix(the_neuron,i)/(sigma*(radius+1)*(radius+1)) );
    }
  }
  return out;
}

// For dist_type = letremy
NumericVector selectNei_f(int the_neuron, int iteration, double radius, NumericMatrix& distMatrix1, NumericMatrix& distMatrix2, std::string radius_type){
  NumericVector out( distMatrix1.rows() ) ;
  for(int i=0; i<distMatrix1.rows(); i++) {
    if((radius==0.5 && distMatrix1(the_neuron,i)<=1) || distMatrix2(the_neuron,i)<=radius) 
      out[i] = 1.0;
    else 
      out[i] = 0.0;
  }
  return out;
}

// Returns radius
double calculateRadius_f(IntegerVector dim, std::string radius_type, int ind, int maxit) {
  double r = 0;
  if(radius_type.compare("letremy")==0) {
    int r0 = max(floor(dim/2));
    double k = 4*(r0-1.0)/maxit;
    int a = std::floor(maxit/2);
    int b = std::floor(maxit*3/4);
    r = std::ceil(r0/(1+k*ind));
    if(ind==1) r = r0;
    else if(ind>=a && ind<b) r = 0.5;
    else if (ind>=b) r = 0;
  } else if(radius_type.compare("gaussian")==0) {
    double r0 = 1 + (2/3.0)*std::sqrt(dim[0]*dim[0] + dim[1]*dim[1]);
    r = r0*std::pow(1/r0, (ind-1)/(maxit-1.0))-1;
  }
  return r; 
}

// Returns closest prototype for a input observations CC
// [[Rcpp::export]]
IntegerVector obsAffectation_F(NumericMatrix& BB, NumericMatrix& CC) {
  const MapMat Prototypes(as<MapMat>(BB));
  const MapMat New(as<MapMat>(CC));
  Eigen::MatrixXd dist1((Prototypes * New.transpose()));
  Eigen::VectorXd dist2((Prototypes * Prototypes.transpose()).diagonal());
  IntegerVector winners(New.rows());
  for(int i=0; i<dist1.cols(); i++) {
    int current_min = 0;
    double min_dist = -2*dist1(0,i) + dist2(0);
    for(int j=1; j<dist1.rows(); j++){
      double dist = -2*dist1(j,i) + dist2(j);
      if( dist < min_dist ) {
        current_min = j;
        min_dist = dist;
      }
    }
    winners[i] = current_min + 1;
  }
  return winners;
}

IntegerVector selectObs_F(int ind, IntegerVector ddim, std::string type) {
  if(type.compare("korresp")==0) {
    if( ind%2 == 0 ){
      return sample((IntegerVector)Range(1,ddim[0]), 1);
    } else {
      return sample((IntegerVector)Range((ddim[0]+1),ddim[0]+ddim[1]), 1);
    }
  } else {
    return sample((IntegerVector)Range(1, ddim[0]), 1);
  }
}

IntegerVector min_element(Eigen::VectorXd& vec) {
  double min = vec[0];
  IntegerVector min_element;
  min_element[0] = 0;
  for (int i=1; i<vec.rows(); i++) {
    if(vec[i]<min){
      min = vec[i];
      min_element[0] = i;
    }
  }
  return min_element;
}

IntegerVector oneObsAffectation_F(NumericMatrix& x_new, NumericMatrix& prototypes, 
                          std::string type, std::string affectation, 
                          std::string radius_type, NumericVector radius, 
                          List the_grid) {
  // obs are indexed from 0
  IntegerVector the_neuron;
  const MapMat Prototypes(as<MapMat>(prototypes));
  const MapMat New(as<MapMat>(x_new));
  if(affectation.compare("standard")==0) {
    if(type.compare("relational")==0) {
      Eigen::VectorXd dist( Prototypes*New - 0.5*(Prototypes*New*Prototypes.transpose()).diagonal() );
      the_neuron = min_element(dist);
    } else {
      the_neuron = obsAffectation_F(prototypes, x_new);
    }
  } else {
    // Heske's Soft Affectation
    /*
    if(type.compare("relational")==0)  {
      
      if(radius_type.compare("letremy")!=0) {
        Eigen::VectorXd the_dist( Prototypes*New - 0.5*(Prototypes*New*Prototypes.transpose()).diagonal() );
        double min_d = std::numeric_limits<double>::max();
        int min_i = 0;
        for(int i=0; i<Prototypes.rows(); i++) {
          Eigen::VectorXd the_nei( as<MapVec>(selectNei(IntegerVector::create(i), the_grid, radius, CharacterVector::create(radius_type), as<CharacterVector>(the_grid["dist.type"]))) );
          double current_d = (the_dist*the_nei).sum();
          if ( current_d < min_d) {
            min_i = i;
            min_d = current_d;
          }
        }
        the_neuron = IntegerVector::create(min_i);
      } else {
        Eigen::VectorXd the_dist( Prototypes*New - 0.5*(Prototypes*New*Prototypes.transpose()).diagonal() );
        double min_d = std::numeric_limits<double>::max();
        int min_i = 0;
        for(int i=0; i<Prototypes.rows(); i++) {
          IntegerVector the_nei = as<IntegerVector>( selectNei(IntegerVector::create(i), the_grid, radius, CharacterVector::create(radius_type), as<CharacterVector>(the_grid["dist.type"])) );
          double current_d = 0;
          for(int j=0; j<the_nei.size(); j++) {
            if(the_nei[j] == true) current_d += the_dist[j];
          }
          if ( current_d < min_d) {
            min_i = i;
            min_d = current_d;
          }
        }
        the_neuron = IntegerVector::create(min_i);
      }
       
    } else {
      if(radius_type.compare("letremy")!=0) {
        Eigen::VectorXd the_dist( (Prototypes-New).array().square().rowwise().sum().sqrt() );
        double min_d = std::numeric_limits<double>::max();
        int min_i = 0;
        for(int i=0; i<Prototypes.rows(); i++) {
          Eigen::VectorXd the_nei( as<MapVec>(selectNei(IntegerVector::create(i), the_grid, radius, CharacterVector::create(radius_type), as<CharacterVector>(the_grid["dist.type"]))) );
          double current_d = (the_dist*the_nei).sum();
          if ( current_d < min_d) {
            min_i = i;
            min_d = current_d;
          }
        }
        the_neuron = IntegerVector::create(min_i);
      } else {
        Eigen::VectorXd the_dist( (Prototypes-New).array().square().rowwise().sum().sqrt() );
        double min_d = std::numeric_limits<double>::max();
        int min_i = 0;
        for(int i=0; i<Prototypes.rows(); i++) {
          LogicalVector the_nei = selectNei(IntegerVector::create(i), the_grid, radius, CharacterVector::create(radius_type), as<CharacterVector>(the_grid["dist.type"]));
          double current_d = 0;
          for(int j=0; j<the_nei.size(); j++) {
            if(the_nei[j] == true) current_d += the_dist[j];
          }
          if ( current_d < min_d) {
            min_i = i;
            min_d = current_d;
          }
        }
        the_neuron = IntegerVector::create(min_i);
      }
    }
     */
  } 
  return the_neuron;
}

NumericMatrix subsetNumMat(NumericMatrix x, IntegerVector index, int row_column) {
  // row_column : 1 for rowwise, 2 for columnwise
  int n = index.size();
  if (row_column == 1) {
    // rowwise
    NumericMatrix out(n, x.cols());
    index = index - 1;
    for (int i = 0; i < n; i++) {
      out(i, _) = x(index[i], _);
    }
    return out;
  } else {
    // columnwise
    NumericMatrix out(x.rows(), n);
    index = index - 1;
    for (int i = 0; i < n; i++) {
      out(_, i) = x(_, index[i]);
    }
    return out;
  }
}

Eigen::MatrixXd subsetEigenMat(Eigen::MatrixXd x, IntegerVector index, int row_column) {
  // row_column : 1 for rowwise, 2 for columnwise
  int n = index.size();
  index = index - 1;
  if (row_column == 1) {
    // rowwise
    Eigen::MatrixXd out(n, x.cols());
    index = index - 1;
    for (int i = 0; i < n; i++) {
      out.row(i) = x.row(index[i]);
    }
    return out;
  } else {
    // columnwise
    Eigen::MatrixXd out(x.rows(), n);
    index = index - 1;
    for (int i = 0; i < n; i++) {
      out.col(i) = x.col(index[i]);
    }
    return out;
  }
}

NumericVector subsetNumVec(NumericVector x, IntegerVector index) {
  int n = index.size();
  NumericVector out(n);
  index = index - 1;
  for (int i = 0; i < n; i++) {
    out[i] = x[index[i]];
  }
  return out;
}

double col_mean( NumericMatrix::Column X ) {
  double result = 0;
  for( int i=0; i<X.size(); i++ ) {
    result += X[i];
  }
  return result/X.size();
}

double col_sd( NumericMatrix::Column X ) {
  double sqsum = 0;
  double mean = col_mean(X);
  for( int i=0; i < X.size(); i++ ) {
    sqsum += (X[i]-mean)*(X[i]-mean);
  }
  return sqsum/(X.size()-1);
}

// [[Rcpp::export]]
NumericMatrix scaleMatrix(NumericMatrix& mat, std::string scaling) {
  NumericMatrix out(mat.rows(), mat.cols());
  if(scaling.compare("unitvar")==0) {
    for(int i=0; i<mat.cols(); i++) {
      out(_,i) = (mat(_,i)-col_mean(mat(_,i)))/col_sd(mat(_,i));
    }
  } else if(scaling.compare("center")==0) {
    for(int i=0; i<mat.cols(); i++) {
      out(_,i) = mat(_,i)-col_mean(mat(_,i));
    }
  } else {
    return mat;
  }
  return out;
}

std::set<int> generateEquiSteps(int start, int end, int n_steps) {
  std::set<int> steps;
  double d = (end-start)/(n_steps-1.0);
  for(int i=0; i<n_steps; i++){
    steps.insert((int)round(start+i*d)); 
  }
  return steps;
}

std::string int_to_str(int num){
  std::stringstream ss;
  ss << num;
  return ss.str();
}

IntegerVector concatenate(IntegerVector vec1, IntegerVector vec2) {
  IntegerVector out(vec1.size()+vec2.size());
  for(int i=0; i<vec1.size(); i++) {
    out[i] = vec1[i];
  }
  for(int i=0; i<vec2.size(); i++) {
    out[i+vec1.size()] = vec2[i];
  }
  return out;
}

// [[Rcpp::export]]
void trainNumeric(NumericMatrix prototypes, List parameters, NumericMatrix x_data, NumericMatrix norm_x_data, List backup, List the_dist){
  
  // Converting arguments to Rcpp
  std::string type = as< std::vector<std::string> >(parameters["type"])[0];
  std::string affectation = as< std::vector<std::string> >(parameters["affectation"])[0];
  std::string radius_type = as< std::vector<std::string> >(parameters["radius.type"])[0];
  List the_grid = as<List>(parameters["the.grid"]);
  std::string dist_type = as< std::vector<std::string> >(the_grid["dist.type"])[0];
  NumericMatrix coords = as<NumericMatrix>(the_grid["coord"]);
  std::string scaling = as< std::vector<std::string> >(parameters["scaling"])[0];
  IntegerVector nb_save = as<IntegerVector>(parameters["nb.save"]);
  IntegerVector dim = as<IntegerVector>(the_grid["dim"]);
  NumericVector eps0 = as<NumericVector>(parameters["eps0"]);
  int maxit = as<IntegerVector>(parameters["maxit"])[0];
  std::string verbose = as< std::vector<std::string> >(parameters["verbose"])[0];
  
  // Setting up iterations for backup and progess
  std::set<int> backup_steps = generateEquiSteps(1,maxit,nb_save[0]) ;
  std::set<int> completion_steps = generateEquiSteps(1,maxit,10);
  
  // Backup Initialization
  List backup_proto;
  NumericMatrix backup_clustering(x_data.rows(), nb_save[0]);
  NumericVector backup_energy(nb_save[0],0.0);
  
  //NumericMatrix coordDists = calculateDistanceMatrix(coords);
  NumericMatrix coordDists, coordDists2;
  coordDists = as<NumericMatrix>(the_dist["one"]);
  if(dist_type.compare("letremy")==0) {
    coordDists2 = as<NumericMatrix>(the_dist["two"]);
  }
  
  for(int ind=1; ind<= maxit; ind++) {
    
    // Step 6 : Randomly choose an observation (indexed from 1)
    IntegerVector rand_ind = selectObs_F(ind, IntegerVector::create(x_data.rows(), x_data.cols()), type );
    NumericMatrix sel_obs(1, x_data.cols());
    sel_obs(0, Rcpp::_) = norm_x_data(rand_ind[0]-1, Rcpp::_ );
    
    // Step 7 : Assignment step
    NumericMatrix cur_obs = sel_obs;
    NumericMatrix cur_prototypes = prototypes;
    
    double radius = calculateRadius_f(dim, radius_type, ind, maxit);
    
    // winner indexed at 1
    IntegerVector winner = as<IntegerVector>(oneObsAffectation_F(cur_obs, 
                                                               cur_prototypes, 
                                                               type, 
                                                               affectation,
                                                               radius_type, 
                                                               NumericVector::create(radius), 
                                                               the_grid) );
    
    //Step 8 : Representation Step
    NumericVector the_nei;
    if(dist_type.compare("letremy")==0) the_nei = selectNei_f(winner[0]-1, ind, radius, coordDists, coordDists2, radius_type);
    else the_nei = selectNei_f(winner[0]-1, ind, radius, coordDists, radius_type);
    double epsilon = 0.3*eps0[0]/(1+0.2*ind/(dim[0]*dim[1]));

    // Update
    for(int i=0; i<prototypes.rows(); i++) {
      prototypes(i,_) = (1-epsilon*the_nei[i])*prototypes(i,_) + epsilon*the_nei[i]*sel_obs(0,_);
    }
    
    //create backup
    if(backup_steps.find(ind)!=backup_steps.end()){
      NumericMatrix out_proto = scaleMatrix(prototypes, scaling);
      colnames(out_proto) = colnames(norm_x_data);
      CharacterVector rnames(out_proto.rows());
      for(int i=0; i<rnames.size(); i++) {
        rnames(i) = int_to_str(i+1);
      }
      rownames(out_proto) = rnames;
      int ind_s = std::distance(backup_steps.begin(),backup_steps.find(ind));
      // todo: check why clone is needed
      backup_proto[int_to_str(ind_s+1)] = clone(out_proto);
      backup_clustering(_,ind_s) = obsAffectation_F(prototypes, x_data);
      IntegerVector f_x = obsAffectation_F(prototypes, norm_x_data);
      for(int i=0; i<norm_x_data.rows(); i++) {
        if(dist_type.compare("letremy")==0) the_nei = selectNei_f(f_x[i]-1, ind, radius, coordDists, coordDists2, radius_type);
        else the_nei = selectNei_f(f_x[i]-1, ind, radius, coordDists, radius_type);
        backup_energy[ind_s] += sum( the_nei*obsDistance(prototypes, norm_x_data, i) );
      }
      backup_energy[ind_s] /= (norm_x_data.rows() * prototypes.rows());
    }
    
    //  Logging completion to Rcout
    if(verbose.compare("true")==0) {
      if(completion_steps.find(ind)!=completion_steps.end()) {
        int p_completion = std::distance(completion_steps.begin(),completion_steps.find(ind))+1;
        Rcout << 10*p_completion << "% done" << std::endl;
      } 
    }
    
  }
  backup["prototypes"] = backup_proto;
  backup["clustering"] = backup_clustering;
  backup["energy"] = backup_energy;
}



// [[Rcpp::export]]
void trainKorresp(NumericMatrix prototypes, List parameters, NumericMatrix x_data, NumericMatrix norm_x_data, List backup, List the_dist){
  
  // Converting arguments to Rcpp
  std::string type = as< std::vector<std::string> >(parameters["type"])[0];
  std::string affectation = as< std::vector<std::string> >(parameters["affectation"])[0];
  std::string radius_type = as< std::vector<std::string> >(parameters["radius.type"])[0];
  List the_grid = as<List>(parameters["the.grid"]);
  std::string dist_type = as< std::vector<std::string> >(the_grid["dist.type"])[0];
  NumericMatrix coords = as<NumericMatrix>(the_grid["coord"]);
  std::string scaling = as< std::vector<std::string> >(parameters["scaling"])[0];
  IntegerVector nb_save = as<IntegerVector>(parameters["nb.save"]);
  IntegerVector dim = as<IntegerVector>(the_grid["dim"]);
  NumericVector eps0 = as<NumericVector>(parameters["eps0"]);
  int maxit = as<IntegerVector>(parameters["maxit"])[0];
  std::string verbose = as< std::vector<std::string> >(parameters["verbose"])[0];
  
  // Setting up iterations for backup and progess
  std::set<int> backup_steps = generateEquiSteps(1,maxit,nb_save[0]) ;
  std::set<int> completion_steps = generateEquiSteps(1,maxit,10);
  
  // Backup Initialization
  List backup_proto;
  NumericMatrix backup_clustering(x_data.rows()+x_data.cols(), nb_save[0]);
  NumericVector backup_energy(nb_save[0],0.0);
  
  //NumericMatrix coordDists = calculateDistanceMatrix(coords);
  NumericMatrix coordDists, coordDists2;
  coordDists = as<NumericMatrix>(the_dist["one"]);
  if(dist_type.compare("letremy")==0) {
    coordDists2 = as<NumericMatrix>(the_dist["two"]);
  }
  
  for(int ind=1; ind<= maxit; ind++) {
    
    // Step 6 : Randomly choose an observation (indexed from 1)
    IntegerVector rand_ind = selectObs_F(ind, IntegerVector::create(x_data.rows(), x_data.cols()), type );
    NumericMatrix sel_obs(1, norm_x_data.cols());
    sel_obs(0, Rcpp::_) = norm_x_data(rand_ind[0]-1, Rcpp::_ );
    
    // Step 7 : Assignment step
    NumericMatrix cur_obs;
    NumericMatrix cur_prototypes;
    if (ind%2==0) {
      cur_obs = subsetNumMat(sel_obs, Range(1, x_data.cols()), 2);
      cur_prototypes = subsetNumMat(prototypes, Range(1, x_data.cols()), 2);
    } else {
      cur_obs = subsetNumMat(sel_obs, Range(x_data.cols()+1, norm_x_data.cols()), 2);
      cur_prototypes = subsetNumMat(prototypes, Range(x_data.cols()+1, norm_x_data.cols()), 2);
    }
    double radius = calculateRadius_f(dim, radius_type, ind, maxit);
    // winner indexed at 1
    IntegerVector winner = as<IntegerVector>( oneObsAffectation_F(cur_obs, 
                                                                 cur_prototypes, 
                                                                 type, 
                                                                 affectation,
                                                                 radius_type, 
                                                                 NumericVector::create(radius), 
                                                                 the_grid) );
    //Step 8 : Representation Step
    NumericVector the_nei;
    if(dist_type.compare("letremy")==0) the_nei = selectNei_f(winner[0]-1, ind, radius, coordDists, coordDists2, radius_type);
    else the_nei = selectNei_f(winner[0]-1, ind, radius, coordDists, radius_type);
    double epsilon = 0.3*eps0[0]/(1+0.2*ind/(dim[0]*dim[1]));
    
    // Update
    for(int i=0; i<prototypes.rows(); i++) {
      prototypes(i,_) = (1-epsilon*the_nei[i])*prototypes(i,_) + epsilon*the_nei[i]*sel_obs(0,_);
    }
    
    //create backup
    if(backup_steps.find(ind)!=backup_steps.end()){
      NumericMatrix out_proto = scaleMatrix(prototypes, scaling);
      colnames(out_proto) = colnames(norm_x_data);
      CharacterVector rnames(out_proto.rows());
      for(int i=0; i<rnames.size(); i++) {
        rnames(i) = int_to_str(i+1);
      }
      rownames(out_proto) = rnames;
      int ind_s = std::distance(backup_steps.begin(),backup_steps.find(ind));
      // todo: check why clone is needed
      backup_proto[int_to_str(ind_s+1)] = clone(out_proto);
      NumericMatrix input_rows = subsetNumMat(subsetNumMat(norm_x_data, Range(1,x_data.cols()), 2), Range(1,x_data.rows()), 1);
      NumericMatrix input_prototypes_rows = subsetNumMat(prototypes, Range(1,x_data.cols()), 2);
      IntegerVector winner_rows = obsAffectation_F(input_prototypes_rows, input_rows );
      NumericMatrix input_cols = subsetNumMat(subsetNumMat(norm_x_data, Range(x_data.cols()+1, norm_x_data.cols()), 2), Range(x_data.rows()+1,norm_x_data.cols()), 1);
      NumericMatrix input_prototypes_cols = subsetNumMat(prototypes, Range(x_data.cols()+1,norm_x_data.cols()), 2);
      IntegerVector winner_cols = obsAffectation_F(input_prototypes_cols, input_cols );
      backup_clustering(_,ind_s) = concatenate(winner_cols, winner_rows);
      IntegerVector f_x = obsAffectation_F(prototypes, norm_x_data);
      
      for(int i=0; i<norm_x_data.rows(); i++) {
        if(dist_type.compare("letremy")==0) the_nei = selectNei_f(f_x[i]-1, ind, radius, coordDists, coordDists2, radius_type);
        else the_nei = selectNei_f(f_x[i]-1, ind, radius, coordDists, radius_type);
        backup_energy[ind_s] += sum( the_nei*obsDistance(prototypes, norm_x_data, i) );
      }
      backup_energy[ind_s] /= (norm_x_data.rows() * prototypes.rows());
    }
    
    //  Logging completion to Rcout
    if(verbose.compare("true")==0) {
      if(completion_steps.find(ind)!=completion_steps.end()) {
        int p_completion = std::distance(completion_steps.begin(),completion_steps.find(ind))+1;
        Rcout << 10*p_completion << "% done" << std::endl;
      } 
    }
    
  }
  backup["prototypes"] = backup_proto;
  backup["clustering"] = backup_clustering;
  backup["energy"] = backup_energy;
}