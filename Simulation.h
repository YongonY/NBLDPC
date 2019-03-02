#pragma once
#define _CRT_SECURE_NO_WARNINGS
#define BP_DECODE 1
#define EMS_DECODE 2
#define MinMax_DECODE 3
#define T_EMS_DECODE 4
#define T_MinMax_DECODE 5
#define OSD_DECODE 6
#define BS_TEMS_DECODE 7

#define Rand_Seq 1
#define All_0_Seq 0

#define Screen_Logo 0
#define Screen_Conf	1
#define File_Conf	2
#define Screen_Head	3
#define File_Head	4
#define File_Begin_Time	5
#define File_End_Time	6
#define Screen_Sim_Data	7
#define Screen_Sim_End_Data	8
#define File_Sim_Data	9
#define File_Sim_End_Data	10

#include <string>
#include <fstream>
#include <ctime>
using namespace std;

class CSimulation
{
public:
	CSimulation(void);
	~CSimulation(void);
	string ProfileFileName;
	string NonBinaryFileName;
	
	int GFq;
	int decodeMethod;
	int parallel;//���ж�
	int maxIter;//�����㷨������������
	int randomMsg;//��Դ���ɵ���ȫ0�������ϵ��

	//Puncture
	int PuntureVarDegree;
	
	//modulation
	int nQAM;
	string ConstellationFileName;

	//channel
	int randomseed;
	double snrBegin;
	double snrStep;
	double snrStop;
	double EbN0;

	//crc
	int crcLen;//��Ϣλ��CRC�ĳ���
	int crc_correctLen;//CRC�����ڸ��������λ��

	//OSD
	int OSD_order;//OSD����,�������㷨Ϊ�����㷨ʱ��OSD_order�Ǹ��ͱ�ʾ�ö�Ӧ�Ľ���ȥ��ǿ�����㷨��
	double OSD_factor;//OSD�ۻ�����
	int OSD_flag;//���OSD_flagΪ������ǿ���������㷨ʱ���ý���������������֮���õ���������ۼ�������

	//ems
	int ems_nm;
	int ems_nc;
	double ems_factor;
	double ems_offset;

	//t-ems
	int tems_nr;
	int tems_nc;
	double tems_factor;
	double tems_offset;

	//bs_tems
	int bs_tems_nm;
	int bs_tems_nc;
	double bs_tems_factor;
	double bs_tems_offset;
	
	//Err
	double errFrame;//��ǰ���ԣ�������֡��
	double errBit;
	double errSym;
	double FER;
	double BER;
	double SER;
	double U_errFrame;//���֡��
	double U_errBit;
	double U_errSym;
	double U_FER;
	double U_BER;
	double U_SER;

	//��������
	double simCycle;//��ǰ���Ե���֡��
	int minSimCycle;//ÿ�����������Ҫ���Ե����ٲ���֡��
	int minErrFrame;//ÿ�����������Ҫ���Եõ����ٴ���֡��
	int U_minErrFrame;//ÿ�����������Ҫ���Եõ��������֡����ʵ�ʴ���ȴ���жԵģ�
	int showSimFrameStep;//����̨�����ˢ��Ƶ��
	
	ofstream fout;
	time_t t;
	clock_t start;
	clock_t stop;
	clock_t status;
	struct tm* current_time;

	int Initial(string profilename);
	int Show(int mode);//����̨����־����ʾ
	int ClearSimuCount(void);
	bool NextSNR(void);
	bool SimulateThisSNR(void);

};

