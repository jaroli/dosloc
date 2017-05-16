#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>

using namespace std;


double a        = 15.35440;	//cell dimensions
double b	= 13.29732;
double c        = 45.40261;
double zsi	= 2*9.40260;	//position of interface
   
const int nband	= 1200;		//number of bands (eigenvalues)
const int nk    = 2;		//number of k-points
const int ne	= 2000;		//number of points in the DOS
const int nx	= 120;		//number of grid point of charge along the x-axis
const int ny	= 108;		//number of grid point of charge along the y-axis
const int nz	= 360;		//number of grid point of charge along the z-axis
  

double sigma	= 0.036;	//standard deviation of Gaussians
double Emin	 = -8.0;
double Emax	 = 10.0;
double dE;			//spacing of the DOS points 

double Eigen[nband][nk];	//array of eigenvalues
double qslice[nz][nband][nk];	//amount of charge in a slice
double gauss[ne][nband][nk];	//matrix of gaussians(E - Eigenval)
double E[ne];			//energy grid for DOS
double dos[nz][ne];		//output, local dos


string s_eigen	= "EIGENVAL";
string s_parchg	= "parchg_orig";//name of directory where PARCHG files are stored
string s_out	= "dosloc-out";	//main output file

void fmt(fstream& file, char c)
{
  switch (c)
  {
     case 'f':
       file.setf(ios::fixed, ios::floatfield);
       file.setf(ios::right);
       file.precision(3);
       file.width(11);
       break;
     case 's':
       file.setf(ios::scientific, ios::floatfield);
       file.setf(ios::right);
       file.precision(4);
       file.width(12);
       break;
     default:
       ;
  }

}


//fill the matrix holding Gaussians(x-mi)
//x = E[i] - dE*0.5, the x points are shifted so that their position is identical as in VASP
//mi = Eigen[j] - potcore
//this is done only to speedup the program
//the matrix is used for every y-point
void calc_gauss()
{
  int i, j, k;
  double norm;


  norm = 1.0/(sigma * pow(2*M_PI, 0.5));
  for (i = 0; i < ne; i++)
    for (j = 0; j < nband; j++)
      for (k = 0; k < nk; k++)
        gauss[i][j][k] = norm * exp( -pow(E[i] - dE*0.5 - (Eigen[j][k]), 2) / (2 * sigma * sigma) );

}


void read_eigenval()
{
  int i, k;
  string index;
  string line;
  fstream file;
 

  file.open(s_eigen.c_str(), ios::in);
  for (i = 0; i < 6; i++)
  {
    getline(file, line);
    //cout << line << endl;
  }

  for (k = 0; k < nk; k++)
  {
    getline(file, line);
    //cout << line << endl;
    getline(file, line);
    //cout << line << endl;
    
    for (i = 0; i < nband; i++)
    {
      file >> index;
      file >> Eigen[i][k];
      getline(file, line); //jump to next line
    }
  }
  file.close();
}

void calc_qslice()
{
  int i, j, k;
  double q;	//charge at a certaint point y (it is not charge acctually but the number in PARCHG file)
  string line;
  string fname;
  fstream file;
  //stringstream oss;
  stringstream ssi, ssk;

  for (i = 0; i < nband; i++)
  for (k = 0; k < nk; k++)
  {

    //form the filename gamma version
    //oss.str(""); oss.width(4); oss.fill('0'); oss << i + 1;
    //fname = s_parchg + "/PARCHG." + oss.str() + ".ALLK";

    //form the filename
    ssi.str(""); ssi.width(4); ssi.fill('0'); ssi << i + 1;
    ssk.str(""); ssk.width(4); ssk.fill('0'); ssk << k + 1;
    fname = s_parchg + "/PARCHG." + ssi.str() + "." + ssk.str();

    //read in the charge along y-axis
    file.open(fname.c_str(), ios::in);
    getline(file, line); //skip the header line
    for (j = 0; j < nz; j++)
    {
      file >> q;
      qslice[j][i][k] = q / (nx*ny*nz);
    }
    file.close();

    if (i == (nband-1) && k == (nk-1)) 
      cout << "parchg band=" << i + 1 << " kpoint=" << k + 1;
  }


/*
  check if we are not missing any electrons, the result should be equal to NELECT (2722)
  q = 0.0;
  for (i = 0; i < 1361; i++)
    for (j = 0; j < nz; j++)
      q = q + qslice[j][i];
  cout << "total charge in all occupied bands = " << q << endl;
*/

}


int main()
{
  int i, j, k, l;
  double sum;
  double Vslice;
  fstream file;

  //cout << "sigma = " << sigma << endl;
  //Vslice = a * *b * (c/nz) * 1.0e-24; //in cm^3
  Vslice = a * b * (c/nz) * 1.0e-3; //in 10^21 cm^3

  //generate energy grid
  dE = (Emax - Emin) / double(ne - 1);
  for (i = 0; i < ne; i++)
    E[i] = Emin + i*dE;

  read_eigenval();
  calc_qslice();
  calc_gauss();

  //smear
  for (l = 0; l < nz; l++)
  {
    for (i = 0; i < ne; i++)
    {
      dos[l][i] = 0.0;
      for (j = 0; j < nband; j++)
        for (k = 0; k < nk; k++)
          dos[l][i] = dos[l][i] + qslice[l][j][k]*gauss[i][j][k] / nk; //two states per eigenvalue !
      dos[l][i] = dos[l][i] / Vslice;
    }
    if (l == (nz-1)) cout << " z grid point =" << l + 1 << endl;
  }


  //////////////////////////////////////////////////write output file 
  //cout << "writing output file" << endl;
  file.open(s_out.c_str(), ios::out);
  for (i = 0; i < nz; i=i+2)
  {
    for (j = 0; j < ne; j=j+2)
    {
      fmt(file, 's'); file << i * (c/nz);
      fmt(file, 's'); file << E[j];
      fmt(file, 's'); file << float(dos[i][j]); //what about scaling ?
      file << endl;
    }
    file << endl;
  }
  //plot the silicon part again
  for (i = 0; i < ceil((zsi*nz)/c) + 1; i=i+2)
  {
    for (j = 0; j < ne; j=j+2)
    {
      
      fmt(file, 's'); file << (i + nz) * (c/nz);
      fmt(file, 's'); file << E[j];
      fmt(file, 's'); file << float(dos[i][j]);
      file << endl;
    }
    file << endl;
  }
  file.close();

/*
  //for testing purposes: check if we get the total dos back
  file.open("total", ios::out);
  for (i = 0; i < ne; i++)
  {
    sum = 0.0;
    for (j = 0; j < nz; j++)
      sum = sum + dos[j][i];

    fmt(file, 's'); file << E[i];
    fmt(file, 's'); file << float(sum/nz); //we are doing an average, so divide by number of points
    file << endl;
  }
  file.close();
*/

  return 0;
}


