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
	int parallel;//并行度
	int maxIter;//迭代算法的最大迭代次数
	int randomMsg;//信源生成的是全0还是随机系列

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
	int crcLen;//信息位中CRC的长度
	int crc_correctLen;//CRC中用于辅助纠错的位数

	//OSD
	int OSD_order;//OSD阶数,若译码算法为其他算法时，OSD_order非负就表示用对应的阶数去加强译码算法。
	double OSD_factor;//OSD累积参数
	int OSD_flag;//如果OSD_flag为正，增强其他译码算法时采用解调器输出软量；反之采用迭代译码的累计软量。

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
	double errFrame;//当前测试，错误总帧数
	double errBit;
	double errSym;
	double FER;
	double BER;
	double SER;
	double U_errFrame;//误检帧数
	double U_errBit;
	double U_errSym;
	double U_FER;
	double U_BER;
	double U_SER;

	//测试条件
	double simCycle;//当前测试的总帧数
	int minSimCycle;//每个信噪比所需要测试的最少测试帧数
	int minErrFrame;//每个信噪比所需要测试得到最少错误帧数
	int U_minErrFrame;//每个信噪比所需要测试得到最少误检帧数（实际错误却被判对的）
	int showSimFrameStep;//控制台结果的刷新频率
	
	ofstream fout;
	time_t t;
	clock_t start;
	clock_t stop;
	clock_t status;
	struct tm* current_time;

	int Initial(string profilename);
	int Show(int mode);//控制台与日志的显示
	int ClearSimuCount(void);
	bool NextSNR(void);
	bool SimulateThisSNR(void);

};

