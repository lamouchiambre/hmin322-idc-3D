#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <iostream>
#include <math.h> 
#include <algorithm>
using namespace std;

Eigen::Vector3d barycentre(Eigen::MatrixXd V ){
  Eigen::Vector3d b;
  std::cout << V.size() << std::endl;
  int size = V.size()/3;

  double x = 0; 
  double y = 0;
  double z = 0;

  for (int i = 0; i < size; i++){
    // std::cout << V(i)<<" "<< V(size+i)<<" "<< V(2*size+i)<< " " << std::endl;
    x +=V(i)/size;
    y +=V(size + i)/size;
    z +=V(2 * size + i)/size;
  }

  b(0) = x;
  b(1) = y;
  b(2) = z;

  // std::cout << b(0) << " " << b(1) << " " << b(2) << std::endl;
  return b;
}

Eigen::Vector3d coordSphere(Eigen::Vector3d p, Eigen::Vector3d g){
  Eigen::Vector3d o;

  o(0) = std::sqrt(std::pow(p(0) - g(0), 2) + std::pow(p(1) - g(1), 2) + std::pow(p(2) - g(2), 2));//pi
  o(1) = acos((p(2)-g(2))/o(0));//theta
  o(2) = atan2((p(1)-g(1)),(p(0)-g(0)));//phi

  // std::cout<<"coord shere " << o(0) << " " << o(1) << " "<<o(2) << " " <<std::endl;

  return o;
}

Eigen::Vector3d coordCart(Eigen::Vector3d o, Eigen::Vector3d g){
  Eigen::Vector3d p;

  p(0) = o(0)*sin(o(1))*cos(o(2))+g(0);
  p(1) = o(0)*sin(o(1))*sin(o(2))+g(1);
  p(2) = o(0)*cos(o(1))+g(2);
  // std::cout<<"coord cart " << p(0) << " " << p(1) << " "<<p(2) << " " <<std::endl;
  return p;
}
template <class T>
vector<T> normalize(vector<T> v){
    // Find the min element
  // cout << "\nMin Element = "
  //   << *min_element(v.begin(), v.end());

  T min = *min_element(v.begin(), v.end());
  T max = *max_element(v.begin(), v.end());

  for (int i = 0; i < v.size(); i++)
  {
    /* code */
    v[i] = (v[i]-min)/(max-min);
  }
  return v;
}

template <class T>
vector<T> normalizeInv(vector<T> v, T min, T max){
  for (int i = 0; i < v.size(); i++)
  {
    /* code */
    v[i] = v[i]*(max-min)+min;
  }
  return v;
}

double mean(vector<double> tab){
  double sum = 0.0;
  for (int i = 0; i < tab.size(); i++)
  {
    sum+=tab[i];
  }
  sum/=tab.size();
  
  return sum;
}

bool get_bit(int num, int position)
{
	bool bit = num & (1 << position);
	return bit;
}


// template <class T>
vector<bool> messageToVectBit(string m){
  int size = m.size();
  vector<bool> bytes;

  for (int i = 0; i < size; i++)
  {
    for(int j = 128; j > 0; j = j/2){
      if(m[i] & j){
        bytes.push_back(1);
        std::cout << "1 ";
      }else{
        bytes.push_back(0);
        std::cout << "0 ";
      }
  }
    cout<<"-";
  }
  cout<<endl;
  return bytes;
}

void binary(unsigned int num){
  for(int i = 128; i > 0; i = i/2){
    if(num & i){
      std::cout << "1 ";
    }else{
      std::cout << "0 ";
    }
  }
  std::cout << std::endl;
}



Eigen::Vector2d minMax(Eigen::MatrixXd V, Eigen::Vector3d g){
  Eigen::Vector3d tmp;
  Eigen::Vector3d tmpO;
  Eigen::Vector2d mM;

  double min, max;

  int size = V.size()/3;
  tmp = {V(0), V(size), V(2*size)};
  tmpO = coordSphere(tmp, g); 

  min = tmpO(0); 
  max = tmpO(0);


  for (int i = 1; i < size; i++){
    tmp = {V(i), V(size+i), V(2*size+i)};
    tmpO = coordSphere(tmp, g); 
    // std::cout<<"----- " << i <<"----- " <<std::endl;
    // std::cout<<"coord tmp  " << tmp(0) << " " << tmp(1) << " "<<tmp(2) << " " <<std::endl;
    coordCart(tmpO, g);
    // std::cout<<"tmp " << tmp <<std::endl;
    // std::cout<<"tmpO " << tmpO <<std::endl;
    // std::cout << "card " << coordCart(tmpO, g) <<std::endl;
    min = std::min(min, tmpO(0));
    max = std::max(max, tmpO(0));
    // std::cout<<"min " << min <<std::endl;
    // std::cout<<"max " << max <<std::endl;

  }

  mM(0) = min;
  mM(1) = max;


  return mM;
}

Eigen::MatrixXd tatouage(Eigen::MatrixXd V, vector<int> m, int bn_bin, double a, double esp = 0.1){

  vector<Eigen::Vector3d> p_vect;
  vector<Eigen::Vector3d> p_vect_sphere;
  vector<Eigen::Vector3d> p_vect_sphere_;

  int size = V.size()/3;

  Eigen::MatrixXd V_out(size,3);


  // calculer le centre de graviter
  Eigen::Vector3d bary = barycentre( V );
  Eigen::Vector2d min_max = minMax( V, bary);
  
  //conversion coord Sphere
  // cout<<"p_vect "<<endl;

  for (int i = 0; i < size; i++)
  {
    Eigen::Vector3d tmp(V(i,0),V(i,1),V(i,2));
    
    p_vect.push_back(tmp);
    p_vect_sphere.push_back(coordSphere(tmp, bary));
    // cout<< tmp <<endl;
    // cout<< "--" <<endl;
    // cout<< p_vect_sphere[i](0) <<endl;
    // cout<< "------" <<endl;


  }


  int N = int(size/bn_bin);
  printf("N %i\n",N);
  printf("bn_bin %i\n",bn_bin);

  vector<vector<double>> bins;
  vector<double> p_poinSphere;
  vector<double> p_poinSphere_message ;

  //isolation de la norme des point
  for (int i = 0; i < size; i++){
    p_poinSphere.push_back(p_vect_sphere[i](0));
    // cout<<p_poinSphere[i]<<endl;
  }
  cout<<"size"<<size <<" "<< "p_poinSphere.size "<<p_poinSphere.size()<<endl;
  //normalization
  p_poinSphere = normalize(p_poinSphere);

  // cout<<"p_poinSphereNormelize ";
  // for(auto &bi : p_poinSphere){
  //   cout<<bi<<" ";
  // }
  //   cout<<endl;


  //division en block de N bins
  for (int i = 0; i < bn_bin; i++)
  {
    vector<double> bin;
    for (int j = 0; j < N; j++)
    {
      bin.push_back(p_poinSphere[i*N + j]);
      // cout<<p_poinSphere[i*N + j]<<endl;
    }
    bins.push_back(bin);   
  }

  cout<<"bins "<<endl;
  for (int i = 0; i < bins.size(); i++)
  {
      cout<<"__";
    for(auto &bi : bins[i]){
      cout<<bi<<" ";
    }

  }
  cout<<endl;
  
  cout<<"##############################"<<endl;

  for (int i = 0; i < bn_bin; i++)
  {
    double kn = 1;
    double mn = mean(bins[i]);
    // int j = 0;
    double pp;
    vector<double> binTmp = bins[i];
    // cout<<"debut"<<endl;
    // cout<<"bins i " <<bins[i][0]<<" "<<bins[i][1]<<endl;
    // cout<<"binTmp " <<binTmp[0]<<" "<<binTmp[1]<<endl;

    
    cout<<"mi - "<<m[i]<<endl;
    kn = 1;
    if (m[i] == 0)
      {
      cout<<"mn - "<<mn<<endl;  
        while (mn > 0.5-a)
        { 
          // cout<<"mn0 - "<<mn<<" kn - "<<kn<<endl;  

          for (int j = 0; j < N; j++)
          { 
            binTmp[j] = pow(binTmp[j],kn);            
          }
          // cout<<"binTmp "<<binTmp[0]<<" "<<binTmp[1]<<endl;
          mn = mean(binTmp);
          kn+=esp;

        }
        for (int j = 0; j < N; j++)
          { 
            p_poinSphere_message.push_back(binTmp[j]);
          }

      }else
      {
        cout<<"mn - "<<mn<<endl;  

        while (mn < 0.5+a)
        {
          // cout<<"mn1 - "<<mn<<" kn - "<<kn<<endl;  
          for (int j = 0; j < N; j++)
          { 
            binTmp[j] = pow(binTmp[j],kn);            
          }
          // cout<<"binTmp "<<binTmp[0]<<" "<<binTmp[1]<<endl;
          mn = mean(binTmp);
          kn-=esp;
        }
        for (int j = 0; j < N; j++)
          { 
            p_poinSphere_message.push_back(binTmp[j]);
          }
      }
    
  }
  //denormalization premiere etape
  p_poinSphere_message = normalizeInv(p_poinSphere_message, min_max(0), min_max(1));

  //denormalization deuxieme etape
  for (int i = 0; i < size; i++)
  {
    //p_poinSphere_message[i]
    //p_vect_sphere[i](0)
    p_vect_sphere_.push_back(Eigen::Vector3d(p_poinSphere_message[i], p_vect_sphere[i](1), p_vect_sphere[i](2)));
  }
  // vector<Eigen::Vector3d
    for (int i = 0; i < size; i++)
  {
    Eigen::Vector3d tmpCard = coordCart(p_vect_sphere_[i], bary);
    V_out(i,0) = tmpCard(0);
    V_out(i,1) = tmpCard(1);
    V_out(i,2) = tmpCard(2);
  }

  return V_out;
}


int main(int argc, char *argv[])
{

  vector<int> tmpBit = {0,1,1,0};
 

  // Inline mesh of a cube
  const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
    0.0,0.0,0.1,//0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
    1,7,5,
    1,3,7,
    1,4,3,
    1,2,4,
    3,8,7,
    3,4,8,
    5,7,8,
    5,8,6,
    1,5,6,
    1,6,2,
    2,6,8,
    2,8,4).finished().array()-1;

  Eigen::MatrixXd V2 = tatouage(V,  tmpBit,  tmpBit.size(), 0.0);

  for (int i = 0; i < (V2.size()/3); i++)
  {
    /* code */
    cout<<"V2 "<<V2(i,0)<<" " <<V2(i,1)<<" "<<V2(i,2)<<endl;
  }

  // Plot the mesh
  // viewer.data().set_mesh(V, F);

  // viewer.data(0).clear();
  // viewer.data(0).set_mesh(V1, F1);
  Eigen::MatrixXd V1;
  Eigen::MatrixXi F1;
  igl::readOFF("socket.off", V1, F1);

  igl::opengl::glfw::Viewer viewer;

  // viewer.data().set_mesh(V2, F);
  viewer.data().set_mesh(V1, F1);

  // viewer.data().set_face_based(true);
  viewer.launch();
}
