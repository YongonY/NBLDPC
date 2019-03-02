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
	int* VarDegree;//ÿ�������ڵ�Ķ���
	int* ChkDegree;//ÿ��У��ڵ�Ķ���
	int** VarLink;//ÿ�������ڵ�����У��ڵ�����
	int** ChkLink;//ÿ��У��ڵ����ӱ����ڵ�����
	int** VarLinkGFe;//ÿ�������ڵ����ӱߵ�Ȩֵ
	int** ChkLinkGFe;//ÿ��У��ڵ����ӱߵ�Ȩֵ
	int** VarLinkDc;//ÿ�������ڵ�����������У��ڵ�ĵڼ���
	int** ChkLinkDv;//ÿ��У��ڵ����������ӱ����ڵ�ĵڼ���
	

	bool Initial(CSimulation &sim);	
	//ͨ�ø�������
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

	int Decoding(double** L_ch, int* DecodeOutput,int* RelySeri_symbol,int* RelySeri_bit);//ѡ������ģʽ
	int DecideLLRVector(double* LLR);//�����о�

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
	double OSD_Alpha;//�ۻ�����
	int OSD_Order;//OSD����,�������㷨Ϊ�����㷨ʱ��OSD_order�Ǹ��ͱ�ʾ�ö�Ӧ�Ľ���ȥ��ǿ�����㷨��
	int OSD_flag;//���OSD_flagΪ������ǿ���������㷨ʱ���ý���������������֮���õ���������ۼ�������
	int min_distance = 1000000;//��ʼŷʽ��
	int CodeLen_bit;
	int ChkLen_bit;

	double** LLR_Sym_sum;//����ÿ�������ţ������ڵ���ۻ����Ŷ�
	double*** LLR_Sym;//�����ʼ�����Լ�ÿ�ε������������
	double *L_ch_bit;

	int* infor_record;//�ο���Ϣ��ѡȡ����������Ӳ�н�����ߵ����㷨���ۻ������
	int* infor_record_bit;
	int *ReliableSeri_BIT;//�ο���Ϣ���ݵ�����
	int* OSD_temp_out_bit;//ÿһ�η�ת���OSD���
	int* near_optimal_bit;//����ŷʽ������С��OSD���
	int* near_optimal;

	int** H_origin_bit;//�洢ԭʼ�ı��ؼ�У�����
	//int** H_OSD_bit;//��˹��ȥ���У�����
	int** G_OSD_bit;//���ؼ����ɾ���ı�׼ʽ


	int Decoding_OSD_bit(double** L_ch, int* DecodeOutput, int*ReliableSeri_BIT, int flag);//���ؼ�OSD�㷨
    void MatrixH2G(int** H_matrix,int** G_matrix, int H_ChkLen, int H_CodeLen);
	void OSD_permute(int* ReliableSeri_BIT);
	int paritycheck(int* DecodeOutput);
	void G_GaussEliminate_bit(int**Matrix_in,int** MAtrix_out, int col_length, int row_length, int*order_seri);
	void H_GaussEliminate_bit(int**Matrix_in, int** MAtrix_out, int col_length, int row_length, int*order_seri);
	void OSD_Encode_bit(int**Matrix_in,int*seq_in,int* seq_out,int*relyseri_bit,int col_length,int row_length);
	void compute_min_distance_bit();
	
	//CRC
	int CRCLen;//CRC����
	int CRCrate;//CRC�����������λ��
	int **CRCMatrix_G;//��Դ���ɵı������е�CRC�������
	int **CRCMatrix_G_IP;//CRCMatrix_G�ı�׼��ʽ
	int **CRCMatrix_partH_PI;//
	int **H_OSD_CRC_bit;//�����CRC�������У�����
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
	//�ϲ���Ҫ��һЩ����
	double tems_route_num = 0;
	int TX_code[10000];	
};
