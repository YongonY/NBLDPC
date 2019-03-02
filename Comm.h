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
	double SetEbN0(CSimulation &sim, int parallel_order);//�µ�����Ȳ���ʱ�ĳ�ʼ��
//****************************************************************************************************
	int Transmission(void);//������һ��������ͨ�Ź��̣�
	//GenerateMessage������CRC��Encode��modulate������channel��
	//demodulate��Decode��Analysis Err

	// GenerateMessage
	int randomMsg;//����������л���ȫ0����
	int* TX_MSG_BIT_beforeCRC;//��Դ����������������л���ȫ0����
	int regPN[11];
	int GenerateMessage(void);	
	int GenPN(void);
	
	//CRC
	int crcLen;//CRC����
	int crcLen_correct;//CRC���������λ��
	int MSG_SYM_LEN;
	int MSG_BIT_LEN;
	int* TX_MSG_SYM;
	int* TX_MSG_BIT;//CRC����������
	
	void CRCEncode(int *seqOut, int *seqIn, int seqInLen, int crcLen, int crc24Type);
	int CrcCheck(int *seqIn, int seqInLen, int crcLen, int crc24Type);

	//Encode
	int* TX_CODE_SYM;//����������
	int* TX_CODE_BIT;
	int CODE_SYM_LEN;
	int CODE_BIT_LEN;
	int Encode(void);//����ı��뺯��������NBLDPC��

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
	CComplex* RX_MOD_SYM;//ע�����Ǿ����������ŵ��ı�������
	double sigma_n;
	int Channel_AWGN(void);

	//demodulate
	double* RX_LLR_BIT;//��������������
	double** RX_LLR_SYM;
	double **RX_DEMOD_LLR_SYM;
	int Demodulate(void);

	//Decode
	int* RX_DECODE_SYM;//����������
	int* RX_DECODE_BIT;
	int* RX_MSG_SYM; 
	int* RX_MSG_BIT; 
	int* ReliableSeri_Symbol;//ÿ�����ŵĿ��Ŷ�(�������)
	int* ReliableSeri_BIT;
	bool DecodeCorrect;//����������������׼������
	int Decode(void);//��������뺯��������NBLDPC��

	//Err
	int Err(CSimulation &sim);
		
};

