#pragma once
#include "Simulation.h"
#include "Complex.h"
#include "GF.h"
#include "Rand.h"
#include <string>
#include <iostream>
#include <fstream>
#include "nbldpc.h"
#include <cmath>
using namespace std;

class CComm
{
public:
	CComm(void);
	~CComm(void);
	

	int GFq;
	int parallel_num;
	int Bit_Len_PerSYM;
	int tems_route_num = 0;

	CNBLDPC NBLDPC;
	CRand Rand;
	bool Initial(CSimulation &sim, int parallel_order);
	double SetEbN0(CSimulation &sim, int parallel_order);//新的信噪比测试时的初始化
//****************************************************************************************************
	int Transmission(void);//以下是一个完整的通信过程：
	//GenerateMessage、插入CRC、Encode、modulate、经历channel、
	//demodulate、Decode、Analysis Err

	// GenerateMessage
	int randomMsg;//产生随机序列还是全0码字
	int* TX_MSG_BIT_beforeCRC;//信源产生的随机码字序列或者全0序列
	int regPN[11];
	int GenerateMessage(void);	
	int GenPN(void);
	
	//CRC
	int crcLen;//CRC长度
	int crcLen_correct;//CRC用作纠错的位数
	int MSG_SYM_LEN;
	int MSG_BIT_LEN;
	int* TX_MSG_SYM;
	int* TX_MSG_BIT;//CRC编码后的序列
	
	void CRCEncode(int *seqOut, int *seqIn, int seqInLen, int crcLen, int crc24Type);
	int CrcCheck(int *seqIn, int seqInLen, int crcLen, int crc24Type);

	//Encode
	int* TX_CODE_SYM;//编码后的序列
	int* TX_CODE_BIT;
	int CODE_SYM_LEN;
	int CODE_BIT_LEN;
	int Encode(void);//具体的编码函数定义在NBLDPC中

	//puncture
	int* PUN_SYM; int PUN_SYM_LEN;
	int* PUN_BIT; int PUN_BIT_LEN;
	double CodeRate;
	int Puncture(void);

	//modulate
	CComplex* CONSTELLATION;
	CComplex* TX_MOD_SYM;
	int modOrder;
	int MOD_BIT_PER_SYM;
	int MOD_SYM_LEN;
	int* TX_MOD_BIT; int MOD_BIT_LEN;
	int Modulate(void);
	
	//channel
	CComplex* RX_MOD_SYM;//注意这是经历过噪声信道的比特序列
	double sigma_n;
	int Channel_AWGN(void);

	//demodulate
	double* RX_LLR_BIT;//解调后输出的软量
	double** RX_LLR_SYM;
	double **RX_DEMOD_LLR_SYM;
	int Demodulate(void);

	//Decode
	int* RX_DECODE_SYM;//译码输出结果
	int* RX_DECODE_BIT;
	int* RX_MSG_SYM; 
	int* RX_MSG_BIT; 
	int* ReliableSeri_Symbol;//每个符号的可信度(序号序列)
	int* ReliableSeri_BIT;
	bool DecodeCorrect;//迭代译码书否输出标准的码字
	int Decode(void);//具体的译码函数定义在NBLDPC中

	//Err
	int Err(CSimulation &sim);
		
};

