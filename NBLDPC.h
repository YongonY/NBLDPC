#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include "gf.h"
#include "Simulation.h"
#include<vector>
#include<algorithm>
#include<cmath>
using namespace std;

class CNBLDPC
{
public:
	CNBLDPC(void);
	~CNBLDPC(void);
	
	// Galois Field
	int GFq;
	CGF GF;
	int maxIter;

	// Code Parameter
    int CodeLen;
	int ChkLen;
	int* PuncturePosition;
	int PunctureLen;

	int maxVarDegree;
	int maxChkDegree;
	int* VarDegree;//每个变量节点的度数
	int* ChkDegree;//每个校验节点的度数
	int** VarLink;//每个变量节点连接校验节点的序号
	int** ChkLink;//每个校验节点连接变量节点的序号
	int** VarLinkGFe;//每个变量节点连接边的权值
	int** ChkLinkGFe;//每个校验节点连接边的权值
	int** VarLinkDc;//每个变量节点是其所连接校验节点的第几度
	int** ChkLinkDv;//每个校验节点是其所连接变量节点的第几度
	

	bool Initial(CSimulation &sim);	
	//通用辅助函数
	int ExchRow(int** H, int source, int destination);
	int AddRow(int** H, int source, int destination, int multifactor);
	int ExchCol(int** H, int source, int destination);
	int CopyLLRVector(double* L_d, double* L_s);
	int ClearLLRVector(double* L);
	int AddLLRVector(double* L, double* L1, double* L2);
	int MinusLLRVector(double* L, double* L1, double* L2);

	//Encode
	int* EncodeExchangeColSource;
	int* EncodeExchangeColDest;
	int EncodeExchangeColNum;
	int** EncodeVarLink;
	int** EncodeVarLinkGFe;
	int* EncodeVarDv;
	
	int InitialEncode(void);
	int Encode(int* msg_sym, int* code_sym);
	int GaussEliminate(int** H, int* ExchangeColSource, int* ExchangeColDest, int& ExchangeColNum);

	//******************************************Decode***********************************************
	double*** L_v2c;
	double*** L_c2v;
	double* LLRV;
	double** L_post;
	int DecodeMethod;

	int Decoding(double** L_ch, int* DecodeOutput,int* RelySeri_symbol,int* RelySeri_bit);//选择译码模式
	int DecideLLRVector(double* LLR);//译码判决

	//BP
	double* L_sigma;
	double* L_rho;

	int Decoding_BP(double** L_ch, int* DecodeOutput);
	int LLR_BoxPlus(double* L, double* L1, double* L2, int A1, int A2);
	int L_Back(int row, int l);
	int L_Forward(int row, int l);

	
	// EMS
	int EMS_Nm;
	int EMS_Nc;

	double EMS_Correction_Factor;
	double EMS_Correction_Offset;
	double** EMS_sort_L_v2c;
	int** EMS_sort_Entr_v2c;
	double* EMS_L_c2v;
	int Decoding_EMS(double** L_ch, int* DecodeOutput,int*RelySeri_bit);
	int SortLLRVector(double* L_sorted, int* Enre_sorted, int Len_sorted, double* L_sorting, int Len_sorting);
	int ConstructConf(int Nm, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int except, int end, int row);
	

	// Trellis EMS
	double TEMS_Factor;
	double TEMS_Offset;
	int TEMS_Nr;
	int TEMS_Nc;

	double** TEMS_deltaU;
	int* TEMS_Beta;
	int TEMS_Syndrome;
	int** TEMS_Min;
	bool** TEMS_inConf_Nr;
	double* TEMS_deltaW;
	int** TEMS_Eta;
	int* TEMS_Eta_Candidate;
	bool* TEMS_isSelected;
	double* TEMS_Lc2v;
	bool* TEMS_isUpdate;

	int Decoding_TEMS(double** L_ch, int* DecodeOutput,int* RelySeri_symbol,int* RelySeri_bit);
	int TEMS_Get_deltaU(int row);
	int TEMS_Get_Beta(int row);
	int TEMS_Get_Min(int row);
	int TEMS_ConstructConf(int Nr, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int end);

	//OSD	
	//Parameter
	double OSD_Alpha;//累积参数
	int OSD_Order;//OSD阶数,若译码算法为其他算法时，OSD_order非负就表示用对应的阶数去加强译码算法。
	int OSD_flag;//如果OSD_flag为正，增强其他译码算法时采用解调器输出软量；反之采用迭代译码的累计软量。
	int min_distance = 1000000;//初始欧式距
	int CodeLen_bit;
	int ChkLen_bit;

	double** LLR_Sym_sum;//计算每个（符号）变量节点的累积可信度
	double*** LLR_Sym;//保存初始软量以及每次迭代的输出软量
	double *L_ch_bit;

	int* infor_record;//参考信息，选取解调软输出的硬判结果或者迭代算法的累积软输出
	int* infor_record_bit;
	int *ReliableSeri_BIT;//参考信息依据的软量
	int* OSD_temp_out_bit;//每一次翻转后的OSD输出
	int* near_optimal_bit;//保存欧式距离最小的OSD输出
	int* near_optimal;

	int** H_origin_bit;//存储原始的比特级校验矩阵
	//int** H_OSD_bit;//高斯消去后的校验矩阵
	int** G_OSD_bit;//比特级生成矩阵的标准式


	int Decoding_OSD_bit(double** L_ch, int* DecodeOutput, int*ReliableSeri_BIT, int flag);//比特级OSD算法
    void MatrixH2G(int** H_matrix,int** G_matrix, int H_ChkLen, int H_CodeLen);
	void OSD_permute(int* ReliableSeri_BIT);
	int paritycheck(int* DecodeOutput);
	void G_GaussEliminate_bit(int**Matrix_in,int** MAtrix_out, int col_length, int row_length, int*order_seri);
	void H_GaussEliminate_bit(int**Matrix_in, int** MAtrix_out, int col_length, int row_length, int*order_seri);
	void OSD_Encode_bit(int**Matrix_in,int*seq_in,int* seq_out,int*relyseri_bit,int col_length,int row_length);
	void compute_min_distance_bit();
	
	//CRC
	int CRCLen;//CRC长度
	int CRCrate;//CRC中用作纠错的位数
	int **CRCMatrix_G;//信源生成的比特序列的CRC编码矩阵
	int **CRCMatrix_G_IP;//CRCMatrix_G的标准形式
	int **CRCMatrix_partH_PI;//
	int **H_OSD_CRC_bit;//添加了CRC辅助后的校验矩阵
	void CRCMatrixGen();
	


	//Basic_Set Conf
	int length_basicSet;
	int basic_set_Nc;
	int basic_set_Nm;
	double basic_set_factor;
	double basic_set_offset;
	int* Set_q;
	int* Set_col;
	int* basicSet_q;
	int* basicSet_col;
	double* basicSet_L;
	bool* basicSet_available;
	bool* col_isSelected;
	int* basicSet_Candidate;

	int Decoding_BS_TEMS(double** L_ch, int* DecodeOutput, int* RelySeri_symbol, int* RelySeri_bit);
	int BasicSet_TEMS_ConstructConf(int begin, int end, int &diff,int& sumNonele, double& sumNonLLR);
	void Construct_BasicSet(int Nm);

	//test_Conf
	double** TEMS_sort_deltaU;
	int **TEMS_sort_q;
	int test_ConstructConf(int Nm, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int end);
	void test_sort();

	//////////////////////////////////////////////////////////////////////
	//上层需要的一些数据
	double tems_route_num = 0;
	int TX_code[10000];	
};
