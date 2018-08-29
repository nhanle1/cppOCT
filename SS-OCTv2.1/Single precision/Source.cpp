// octProc.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <numeric> //for iota interation
#include <vector>
#include <fstream>
#include <stdio.h>

#include <algorithm> //provide function 'lower_bound' for interpolation
#include <complex>
#include <mkl.h>
#include <mkl_dfti.h>
#include <ctime>
#include <windows.h>

using namespace std;

const char* octfilename("2017.10.20_14.41.31_s03F_G1_8x8mm_r8.oct");
const char* refname("ref_gumdata.bin");
//Function declaration
float polyval(vector<float>, int, float);
//interpolation
void dinterp1(float* data, int nrows, float* x, int N, float* result);
//void fft_complex(vector<complex<float>>& in,size_t NZ);
void sumstats(vector<INT64>& time);
void printlogo();

DFTI_DESCRIPTOR_HANDLE descriptor;

int main() {
	printlogo();

	ofstream outputfile("ED.bin", ios::binary);

	const long int nZ = 2048;
	const long int nX = 500;
	const long int nY = 500;
	const long int nR = 8;
	const long int nXR = nX * nR;
	const long int nBF = nXR * nZ;
	const long int nXZ2 = nX * nZ / 2;
	size_t status; //check fread

	LARGE_INTEGER begintime, endtime1, endtime2, frequency;
	vector<INT64> time(8);

	vector<float> refdata(nZ), KES(nZ);		//data reference&K-estimate
	vector<complex<float>> specC(nBF/2);
	vector<float> specS(nBF / 2,0), EDdataS(nBF / 2, 0);
	vector<float> buffer(nBF), tbuf(nBF);	//temporary buffer
	vector<float> idenv(nXR, 1); //Identity matrix for PCA
	vector<INT16> buf(nBF);		//data buffer for binary read
	FILE* octfile = NULL;


	//Linear interpolation parameter initialization
	//'interp1' equivalent to matlab (spline order 0 interp)
	//MKL_INT nscoeff;                    // number of spline coefficients
	MKL_INT dorder[1] = { 1 };
	/***** Parameters describing function *****/
	DFTaskPtr task;
	/* Limits of interpolation interval are provided in case
	of uniform partition */
	float K[2]; //Linear K-clock
	K[0] = 1; K[1] = nZ;

	/***** Parameters describing boundary conditions type *****/
	/* No boundary conditions are provided for linear spline */
	float *bc = 0;
	// pointer array of spline coefficients
	vector<float> scoeff(nX*nR*(nZ - 1)*DF_PP_LINEAR);
	float* rptr = new float[nBF]; //return pointer of interpolation results

	MKL_INT xhint = DF_UNIFORM_PARTITION;	//Uniform K-clock partition from K[0] to K[1]
	MKL_INT yhint = DF_MATRIX_STORAGE_ROWS; //data stored in rows, i.e. data(i+j*nZ);
	dfsNewTask1D(&task, nZ, K, xhint, nX*nR, buffer.data(), yhint);
	dfsEditPPSpline1D(task, DF_PP_LINEAR, DF_PP_DEFAULT, DF_NO_BC,
		0, DF_NO_IC, 0, scoeff.data(), DF_NO_HINT);
	//***** Construct linear spline using STD method ****

	//=============FFT initialization===================//


	/*===========================PCA initialization======================*/
	MKL_INT         m;				//Null value 
	MKL_INT ifail[nR], info;
	MKL_Complex8   alpha = { 1,0 }, beta = { 0,0 }, inverse = { -1,0 };
	MKL_Complex8  *matrixc = new MKL_Complex8[nR*nR];
	float w[nR];					//array of eigen values
	const long int RFilt = 4;
	MKL_Complex8 z[nR*nR] = { 0 }, ImEV[nR*nR] = { 0 }, complexone = { 1,0 };
	vector<float> Cov(nR*nR, 0), eigvec(nR*nR, 0), ImCov(nR*nR, 0);

	//========== Initialize CLOCK performance & read reference data===============/

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&begintime);
	time[7] = frequency.QuadPart;

	//OPEN REF file
	ifstream ref, kesbin;
	ref.open(refname, ios::in | ios::binary);
	kesbin.open("KES.bin", ios::in | ios::binary);
	//check to see that the file was opened correctly:
	if (ref.is_open() || kesbin.is_open()) {
		//SAVE REFERENCE DATA TO ARRAY
		ref.read((char*)refdata.data(), nZ * sizeof(float));
		ref.close();
		//Save KES data
		kesbin.read((char*)KES.data(), nZ * sizeof(float));
		kesbin.close();
	}
	else {
		std::cerr << "There was a problem"
			"opening the ref-data & KES file!\n";
		exit(1);//exit or do additional error checking
	}

	//========== MAIN PROCESS BEGIN===============/

	for (int iY = 0; iY < nY; iY++) {

		QueryPerformanceCounter(&endtime1);
		//OPEN OCT file
		if (octfile == NULL) {
			if (fopen_s(&octfile, octfilename, "rb") == 0) {
				fseek(octfile, 20580, SEEK_SET);
				cout << "\t\t Size of DATA: " << nZ << 'X' << nX
					<< 'X' << nR << 'X' << nY <<
					"\tProgress" << "\n\t\t\t\t\t";
			}
			else {
				std::cerr << "There was a problem opening the oct file!\n";
				exit(2);//exit or do additional error checking
			}
		}

		for (int i = 0; i < nR; i++) {
			fseek(octfile, 4, SEEK_CUR);
			fread(buf.data() + i * nX*nZ, sizeof(short), nX*nZ, octfile);
		};
		printf("\t%.1f %%", ((float)(iY + 1) / (float)nY * 100));


		QueryPerformanceCounter(&endtime2);
		time[0] += endtime2.QuadPart - endtime1.QuadPart;



		/*-----Checked 2018.07.16 functioning correctly----------*/
		//Assigning buffer typedef of usignedshort* -> float*

		std::copy(buf.begin(), buf.end(), buffer.begin());
		for (int i = 0; i < nR; i++) {
			for (int j = 0; j < nX; j++) {
				cblas_saxpy(nZ, -1, refdata.data(), 1, 
					&buffer[j*nZ + i * nZ*nX], 1);
			}
		}


		QueryPerformanceCounter(&endtime1);
		time[1] += endtime1.QuadPart - endtime2.QuadPart;



		/************** Create Data Fitting task ***********************/
		/***************************************************************/
		dfsConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
		//Interpolate, careful with STORAGE_COLS and STORAGE_ROWS
		dfsInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
			nZ, KES.data(), DF_NON_UNIFORM_PARTITION, 1,
			dorder, NULL, tbuf.data(), DF_MATRIX_STORAGE_ROWS, NULL);

		QueryPerformanceCounter(&endtime2);
		time[2] += endtime2.QuadPart - endtime1.QuadPart;


		/*************************** FOURIER TRANSFORM  *****************/
		/***************************************************************/
		DftiComputeForward(descriptor, buffer.data(), specC.data());
		QueryPerformanceCounter(&endtime1);
		time[3] += endtime1.QuadPart - endtime2.QuadPart;




		/*--------------------EIGEN DECOMPOSITION: PCA METHOD--------------*/
		//Construct covariance matrix nRxnR
		vcAbs(nBF / 2,(MKL_Complex8*) specC.data(), specS.data());
		cblas_sgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, nR, nR, nXZ2, 1,
			specS.data(), nXZ2, specS.data(), nXZ2, 0, Cov.data(), nR);

		/*cblas_cgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, nR, nR, nXZ2, &alpha,
			specC.data(), nXZ2, specC.data(), nXZ2, &beta, matrixc, nR);*/

			//Eigen value of first 'RFilt' eigen value, the eigen value is ascending
		info = LAPACKE_ssyevx(LAPACK_COL_MAJOR, 'V', 'I', 'L', nR, Cov.data(), nR,
			NULL, NULL, nR - RFilt + 1, nR, -1, &m, w, eigvec.data() + (nR - RFilt)*nR, nR, ifail);

		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, nR, nR, nR, -1,
			eigvec.data(), nR, eigvec.data(), nR, 1, ImCov.data(), nR);

		//projection of filter eigenvector to original data
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nXZ2, nR, nR, 1,
			specS.data(), nXZ2, ImCov.data(), nR, 0, EDdataS.data(), nXZ2);
		//Absolute data & summation, take 9% processing time
		vsAbs(nBF / 2, EDdataS.data(), tbuf.data());
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nXZ2, 1, nR, 1,
			tbuf.data(), nXZ2, idenv.data(), nR, NULL, buffer.data(), nXZ2);

		QueryPerformanceCounter(&endtime2);
		time[4] += endtime2.QuadPart - endtime1.QuadPart;

		printf("\b\b\b\b\b\b\b\b");

		outputfile.write((char*)(buffer.data()), sizeof(float)*nXZ2);
		
	};

	time[6] = (endtime2.QuadPart - begintime.QuadPart);
	sumstats(time);



	std::printf("\t --------------------------------------------\n"
		"\t Total time: \t\t %.1f ms \t \n\n\t\t ",
		float(time[6]) * 1000 / time[7]);


	std::fclose(octfile);

	outputfile.close();

	std::system("PAUSE");


	return 0;
}




//POLYVAL function
float polyval(vector<float> p, int n, float x) {
	float px;
	px = 0;

	for (int k = 0; k < n; k++) {
		px += p[k] * pow(x, k);
	}

	return px;
}





void dinterp1(float* data, int nrows, float* x, int N, float* result) {
	for (int i = 0; i < N; i++) {
		// get coordinates of bounding grid locations
		long long x_1 = (long long)std::floor(x[i] - 1);

		// handle special case where x is the last element
		if ((x[i] - 1) == (nrows - 1)) { ; x_1 -= 1; }
		// return 0 for target values that are out of bounds
		if ((x_1 < 0) || (x_1 + 1) > (nrows - 1)) {
			result[i] = 0;
		}
		else {
			// get the array values
			const float& f_1 = data[x_1];
			const float& f_2 = data[x_1 + 1];
			// compute weights
			float w_x1 = x_1 + 1 - (x[i] - 1);
			result[i] = f_1 * w_x1 + f_2 - f_2 * w_x1;
		}
	}
}

void sumstats(vector<INT64>& time) {
	time[5] = (float)(time[6] - time[4] - time[3] - time[2] - time[1] - time[0]);
	std::printf("\n\t ===============================================\n"
		"\t =======        Summary statistics       ======\n"
		"\t Bin time: \t\t %.1f ms \t %.1f %% \n"
		"\t INT to Float time: \t %.1f ms \t %.1f %% \n"
		"\t Interpolation time: \t %.1f ms \t %.1f %% \n"
		"\t FFT time: \t\t %.1f ms \t %.1f %% \n"
		"\t ED time: \t\t %.1f ms \t %.1f %% \n"
		"\t Write time: \t\t %.1f ms \t %.1f %% \n",
		float(time[0]) * 1000 / time[7],
		float(time[0]) / time[6] * 100,
		float(time[1]) * 1000 / time[7],
		float(time[1]) / time[6] * 100,
		float(time[2]) * 1000 / time[7],
		float(time[2]) / time[6] * 100,
		float(time[3]) * 1000 / time[7],
		float(time[3]) / time[6] * 100,
		float(time[4]) * 1000 / time[7],
		float(time[4]) / time[6] * 100,
		float(time[5]) * 1000 / time[7],
		float(time[5]) / time[6] * 100);
}

void printlogo() {
	std::cout << R"(     
                         .wg$$$$$$$$$$$$&gww,
                   ,ag;    *&$&$$$$$$$$$$$$$$$L
                a$$$$$$$W    *T$$$$$$$$$$$$$$$[   $&w
             y$$$$$$$$$$$$&w    *R$$$$$$$$$$$$1   }$$$&g
           a$$$$$$$$$$$$$$$$$b;    7&$$$$$$$$$&   ]$$$$$$w
         z$$$$$$$$$$$$$$$$$$$$$$W    *T$$$$$$$$   ]$$$$$$$$y
        $$$$$$$$$$$$$$$$&&RPMMM*TT*     *R&&$$$r   &$$$$$$$$&
      .PPMfT*TT*!*"                        7&$$L   $$$$$$$$$$&w
                   ,;g                       *T[   8$$$$$$$$$$$L
    ;ggA&&$$$$$$$$$$P*                             $$$$$$$$$$$*
    $$$$$$$$$$$$$&PL                               ]$$$$$$$$C   ,
   ]$$$$$$$$$$$$B*         UW OCTA program         ]$$$$$$C    4$L
   $$$$$$$$$$&$*                                    $$$$F    4$$$&
   $$$$$$$$$$C                 Nhan Le              B&F    g$$$$$$
  j$$$$$&$$*    g      E-mail: nle0601@gmail.com    F`   a$$$$$$$$
   $$$&&$F    4$$         2018.08.23 ver 2.11          a$$$$$$$$$$
   $$$$F    g$$$$U                                   a$$$$$$$$$$$$
   ]$F    g$$$$$$L                                 y&$$$$$$$$$$$$F
    `   g$$$$$$$$[                               ;&$$$$$$$$$$$$$$
      a$&$$$$$$$$1                             ,$$$$$$$$$&&RRPMM*
     ]$$$$$$$$$$$B   ]&w                      'T*|*`
      ?$$$$$$$$$$$   j$$&b;                        ,,;;wwwwgg;`
        &$$$$$$$$$r   $$$$&$w      wwwgg&@$$$$$$$$$$$$$$$$$$P
         ?&$$$$$$$L   &$$$$$$$&w    ?R&$$$$$$$$$$$$$$$$$$$$[
           *&$$$$$[   E$$$$$$$$$$W;    *&$$$$$$$$$$$$$$$BC
             ]T$$$$   $$$$$$$$$$$$&$w    *5$$$$$$$$&$$P*
                ?M$   ]$$$$$$$$$$$$$$$&g    ?R$$$$$M[
                      j$&$$$$$$$$$$$$$$$$W;    TT=
                        `|7*MRR&&$$$&&RRMMT*
     
            )" << "\n\n\n";
//	std::cout << R"(
// .=.                                           
//|                         &@@@@@@@@@@@@@@@@@(                                  
//|                  ,@@@(    ,@@@@@@@@@@@@@@@@@   .                             
//|               ,@@@@@@@@%     @@@@@@@@@@@@@@@   %@@,                          
//|             @@@@@@@@@@@@@@,    ,@@@@@@@@@@@@   *@@@@@                        
//|           &@@@@@@@@@@@@@@@@@%     @@@@@@@@@@,   @@@@@@&                      
//|         #@@@@@@@@@@@@@@@@@@@@@@,    /@@@@@@@(   @@@@@@@@#                    
//|        @@@@@@@@@@@@@@@@@@@&%(,.        @@@@@&   @@@@@@@@@@                   
//|       @@@%(*.                            @@@@   @@@@@@@@@@@                  
//|                                            ,@   &@@@@@@@@@@@                 
//|     /%&@@@@@@@@@@@@%                            #@@@@@@@@@@                  
//|    (@@@@@@@@@@@@@@                              ,@@@@@@@/.                   
//|    @@@@@@@@@@@@@*                                @@@@@@/    @@               
//|   ,@@@@@@@@@@@/              Nhan Le             @@@@&    @@@@,              
//|   (@@@@@@@@@%      E-mail: nle0601@gmail.com     @@@    *@@@@@#              
//|   %@@@@@@@@    %                                 %.   ,@@@@@@@%              
//|   %@@@@@@,   .@@                                     @@@@@@@@@#              
//|   /@@@@#    @@@@       2018.08.23 ver 2.11         %@@@@@@@@@@@              
//|    @@%    &@@@@@,                                *@@@@@@@@@@@/               
//|    @    (@@@@@@@/                               @@@@@@@@@@@@@&               
//|       .@@@@@@@@@&                             @@@@@@@@@@@@@@@                
//|      @@@@@@@@@@@@   %                       %@@@@@&%/.                      
//|      (@@@@@@@@@@@   &@@                                                      
//|       .@@@@@@@@@@   #@@@@/                  .*(%@@@@@@@@@@,                  
//|         @@@@@@@@@   *@@@@@@@     &@@@@@@@@@@@@@@@@@@@@@@@                    
//|          *@@@@@@@,  .@@@@@@@@@,    *@@@@@@@@@@@@@@@@@@@*                     
//|            %@@@@@#   @@@@@@@@@@@%    .@@@@@@@@@@@@@@@&                       
//|              .@@@&   @@@@@@@@@@@@@@,    (@@@@@@@@@@,                         
//|                 /@   @@@@@@@@@@@@@@@@/     @@@@@(                            
//|                      &@@@@@@@@@@@@@@@@@@                                     
//|                           .(&@@@@@@@@@,  
//)";
}