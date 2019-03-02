#include "Comm.h"


CComm::CComm(void)
	: modOrder(0)
	, TX_MSG_SYM(nullptr)
	, TX_MSG_BIT(nullptr)
	, TX_CODE_SYM(nullptr)
	, TX_CODE_BIT(nullptr)
	, TX_MOD_SYM(nullptr)
	, RX_MOD_SYM(nullptr)
	, RX_LLR_BIT(nullptr)
	, RX_LLR_SYM(nullptr)
	, ReliableSeri_Symbol(nullptr)
	, ReliableSeri_BIT(nullptr)
	, RX_DECODE_SYM(nullptr)
	, RX_DECODE_BIT(nullptr)
	, MOD_BIT_PER_SYM(0)
	, MSG_SYM_LEN(0)
	, MSG_BIT_LEN(0)
	, CODE_SYM_LEN(0)
	, CODE_BIT_LEN(0)
	, MOD_SYM_LEN(0)
	, PUN_SYM(nullptr)
	, PUN_SYM_LEN(0)
	, RX_MSG_SYM(nullptr)
	, RX_MSG_BIT(nullptr)
	, PUN_BIT(nullptr)
	, PUN_BIT_LEN(0)
	, TX_MOD_BIT(nullptr)
	, MOD_BIT_LEN(0)
	, sigma_n(0)
	, CodeRate(0)
	, randomMsg(0)
	, CONSTELLATION(NULL)
	, RX_DEMOD_LLR_SYM(NULL)
	, DecodeCorrect(false)
{
}


CComm::~CComm(void)
{
}


//bool CComm::Initial(int GFq, int nQAM, string NBLDPCFileName, int PunctureVarDegree, int maxIter, int randMsg, int deocdemethod, int ems_nm, int ems_nc, double ems_factor, double ems_offset, int randseed)
bool CComm::Initial(CSimulation &sim,int parallel_order)
{	
	GFq = sim.GFq;
	parallel_num = sim.parallel;
	Bit_Len_PerSYM = log(GFq) / log(2);
	
	
	// intial NBLDPC
	NBLDPC.Initial(sim);
	// intial random generator
	Rand.IX = Rand.IY = Rand.IZ = sim.randomseed;
	if (parallel_order) {
		Rand.IX = Rand.IX + parallel_order;
		Rand.IY = Rand.IY + parallel_order;
		Rand.IZ = Rand.IZ + parallel_order;
	}
	//**********************************GenerateMessage&CRC&Encode******************************************

	// intial PN register
	regPN[0] = 1;	regPN[1] = 0;	regPN[2] = 1;	regPN[3] = 0;
	regPN[4] = 0;	regPN[5] = 0;	regPN[6] = 1;	regPN[7] = 1;
	regPN[8] = 0;	regPN[9] = 0;	regPN[10] = 1;
	if (parallel_order) {
		for (int num_yx = 0; num_yx < parallel_order; num_yx++) {
			GenPN();
		}
	}
	randomMsg = sim.randomMsg;
	
	//CRC
	crcLen = sim.crcLen;
	crcLen_correct = sim.crc_correctLen;

	MSG_SYM_LEN = NBLDPC.CodeLen - NBLDPC.ChkLen;
	MSG_BIT_LEN = MSG_SYM_LEN * Bit_Len_PerSYM;
	TX_MSG_BIT = new int[MSG_BIT_LEN];//CRC编码后的序列
	TX_MSG_SYM = new int[MSG_SYM_LEN]; 
	TX_MSG_BIT_beforeCRC = new int[MSG_BIT_LEN - crcLen];//信源产生的随机码字序列或者全0序列
	
   //Encode
	CODE_SYM_LEN = NBLDPC.CodeLen;
	CODE_BIT_LEN = CODE_SYM_LEN * Bit_Len_PerSYM;
	TX_CODE_SYM = new int[CODE_SYM_LEN];//编码后的序列
	TX_CODE_BIT = new int[CODE_BIT_LEN];

	//***************************************** puncture***********************************
	PUN_SYM_LEN = NBLDPC.PunctureLen;
	PUN_BIT_LEN = PUN_SYM_LEN * Bit_Len_PerSYM;
	PUN_SYM = new int[PUN_SYM_LEN]; 
	PUN_BIT = new int[PUN_BIT_LEN];
	CodeRate = double(MSG_SYM_LEN) / double(CODE_SYM_LEN - PUN_SYM_LEN);

	// initial Puncture Position
	for (int p_s_len = 0; p_s_len < PUN_SYM_LEN; p_s_len++)
	{
		PUN_SYM[p_s_len] = NBLDPC.PuncturePosition[p_s_len];
		for (int b_p_s = 0; b_p_s < Bit_Len_PerSYM; b_p_s++)
		{
			PUN_BIT[p_s_len * Bit_Len_PerSYM + b_p_s]
				= Bit_Len_PerSYM * NBLDPC.PuncturePosition[p_s_len] + b_p_s;
		}
	}
	//****************************************modulate&demodulate*****************************************
	modOrder = sim.nQAM;
	MOD_BIT_PER_SYM = log(sim.nQAM * 1.0) / log(2.0);
	MOD_SYM_LEN = (CODE_SYM_LEN - PUN_SYM_LEN) * Bit_Len_PerSYM / MOD_BIT_PER_SYM;
	MOD_BIT_LEN = MOD_SYM_LEN * MOD_BIT_PER_SYM;

	// initial modulation and demodulation constellation
	CONSTELLATION = new CComplex[modOrder];
	ifstream fcons(sim.ConstellationFileName);
	string consrub;
	int index;
	for (int k = 0; k < modOrder; k++)
	{
		fcons >> consrub >> index;
		fcons >> consrub >> CONSTELLATION[index].Real >> consrub >> CONSTELLATION[index].Image;
	}
	fcons.close();
		
	//modulate
	TX_MOD_BIT = new int[MOD_BIT_LEN];
	TX_MOD_SYM = new CComplex[MOD_SYM_LEN];
	
	//demodulate
	RX_MOD_SYM = new CComplex[MOD_SYM_LEN];//注意这是经历过噪声信道的比特序列
	RX_LLR_SYM = new double*[CODE_SYM_LEN];//解调后输出的软量
	RX_LLR_BIT = new double[CODE_BIT_LEN];

	RX_DEMOD_LLR_SYM = new double *[MOD_SYM_LEN];
	for(int s = 0; s < MOD_SYM_LEN; s ++)
	{
		RX_DEMOD_LLR_SYM[s] = new double[modOrder - 1];
	}

	//**************************************Decode&&Err****************************************************
	RX_DECODE_SYM = new int[CODE_SYM_LEN];//译码输出结果
	RX_DECODE_BIT = new int[CODE_BIT_LEN];
	RX_MSG_SYM = new int[MSG_SYM_LEN];
	RX_MSG_BIT = new int[MSG_BIT_LEN];
	ReliableSeri_Symbol = new int[CODE_SYM_LEN];
	ReliableSeri_BIT = new int[CODE_BIT_LEN];
	for (int s = 0; s < CODE_SYM_LEN; s++)
	{
		RX_LLR_SYM[s] = new double[GFq - 1];
	}
	return false;
}

double CComm::SetEbN0(CSimulation &sim, int parallel_order)
{
	// reset the seed
	Rand.IX = Rand.IY = Rand.IZ = sim.randomseed;
	if (parallel_order) {
		Rand.IX = Rand.IX + parallel_order;
		Rand.IY = Rand.IY + parallel_order;
		Rand.IZ = Rand.IZ + parallel_order;
	}
	// reset the PN register
	regPN[0] = 1;	regPN[1] = 0;	regPN[2] = 1;	regPN[3] = 0;
	regPN[4] = 0;	regPN[5] = 0;	regPN[6] = 1;	regPN[7] = 1;
	regPN[8] = 0;	regPN[9] = 0;	regPN[10] = 1;
	if (parallel_order) {
		for (int num_yx = 0; num_yx < parallel_order; num_yx++) {
			GenPN();
		}
	}
	// set sigma
	double EbN0 = pow(10.0, sim.EbN0 / 10.0);
	sigma_n = 1.0 / sqrt(2 * MOD_BIT_PER_SYM * CodeRate * EbN0);
	return sigma_n;
}

int CComm::Transmission()
{
	GenerateMessage();
	Encode();
	Puncture();
	Modulate();
	Channel_AWGN();
	Demodulate();
	Decode();
	return 0;
}

//**********************************************************GenerateMessage&Puncture&Encode************************************************
int CComm::GenerateMessage()
{
	// generate TX_MSG_BIT & TX_MSG_SYM
	for(int b = 0; b < MSG_BIT_LEN-crcLen; b ++)
	{
		if(randomMsg)
		{
			for (int i = 0; i < parallel_num - 1; i++) { GenPN(); }
			TX_MSG_BIT_beforeCRC[b] = GenPN();
		}
		else
		{
			TX_MSG_BIT_beforeCRC[b] = 0;//GenPN();
		}
	}

	//CRC编码
	CRCEncode(TX_MSG_BIT, TX_MSG_BIT_beforeCRC, MSG_BIT_LEN- crcLen,crcLen, 0);
	/////////////////////////////////////////////////////////////////////////////////////////
	/*
	cout << "TX_MSG_bit_beforeCRC" << endl;
	for (int i = 0; i < MSG_BIT_LEN - crcLen; i++) {
		cout << TX_MSG_BIT_beforeCRC[i];
	}
	cout << endl;

	cout << "TX_MSG_bit" << endl;
	for (int i = 0; i < MSG_BIT_LEN; i++) {
		cout << TX_MSG_BIT[i];
	}
	cout << endl;
	*/
	/////////////////////////////////////////////////////////////////////////////////////////


	for(int s = 0; s < MSG_SYM_LEN; s ++)
	{
		TX_MSG_SYM[s] = 0;
		for (int i = 0; i < Bit_Len_PerSYM; i++) {
			TX_MSG_SYM[s] = 2 * TX_MSG_SYM[s] + TX_MSG_BIT[Bit_Len_PerSYM*s+ Bit_Len_PerSYM -1-i];
 		}
			
	}

	return 0;
	
}
int CComm::GenPN(void)
{
	/*shift the shift register*/
	for (int i = 10; i >= 1; i--)
	{
		regPN[i] = regPN[i - 1];
	}
	/*calculate the output*/
	regPN[0] = regPN[10] ^ regPN[3];
	/*output the result*/
	return regPN[10];
}


int CComm::Encode(void)
{
	// TX_MSG_SYM --> TX_CODE_SYM
	if(randomMsg)
	{
		NBLDPC.Encode(TX_MSG_SYM, TX_CODE_SYM);
	}
	else
	{
		for(int s = 0; s < CODE_SYM_LEN; s ++)
		{
			TX_CODE_SYM[s] = 0;
		}
	}
	

	// TX_CODE_SYM --> TX_CODE_BIT
	for(int s = 0; s < CODE_SYM_LEN; s ++)
	{
		for(int b_p_s = 0; b_p_s < Bit_Len_PerSYM; b_p_s ++)
		{
			TX_CODE_BIT[s * Bit_Len_PerSYM + b_p_s]
			= ((TX_CODE_SYM[s] & (1 << b_p_s)) == 0)? 0 : 1;
		}
	}
	
	// TX_CODE_BIT --> TX_MSG_BIT
	for(int b = 0; b < MSG_BIT_LEN; b ++)
	{
		TX_MSG_BIT[b] = TX_CODE_BIT[b];
	}
	
	return 0;

}
int CComm::Puncture(void)
{
	// TX_CODE_BIT --> TX_MOD_BIT
	int p_i = 0;
	int m_b = 0;
	for (int b = 0; b < CODE_BIT_LEN; b++)
	{
		if (PUN_BIT[p_i] == b && PUN_BIT_LEN != 0) // this bit is punctured
		{
			p_i++;
		}
		else
		{
			TX_MOD_BIT[m_b] = TX_CODE_BIT[b];
			m_b++;
		}
	}
	return 0;
}
//***********************************************Modulate&channel&Demodulate********************************************************
int CComm::Modulate(void)
{
	// TX_MOD_BIT --> TX_MOD_SYM
	int constellation_index;
	for(int s = 0; s < MOD_SYM_LEN; s ++)
	{
		constellation_index = 0;
		for(int b_p_s = 0; b_p_s < MOD_BIT_PER_SYM; b_p_s ++)
		{
			constellation_index += (TX_MOD_BIT[s * MOD_BIT_PER_SYM + b_p_s] << ( MOD_BIT_PER_SYM - 1 - b_p_s ));
		}
		TX_MOD_SYM[s].Real = CONSTELLATION[constellation_index].Real;
		TX_MOD_SYM[s].Image = CONSTELLATION[constellation_index].Image;
	}
	return 0;
}


int CComm::Channel_AWGN()
{
	// only support bpsk now
	for(int s = 0; s < MOD_SYM_LEN; s ++)
	{
		RX_MOD_SYM[s].Real = TX_MOD_SYM[s].Real + Rand.Rand_Norm(0, sigma_n);
		RX_MOD_SYM[s].Image = TX_MOD_SYM[s].Image + Rand.Rand_Norm(0, sigma_n);
	}
	return 0;
}


int CComm::Demodulate(void)
{
	if(modOrder == 2)
	{
		// RX_MOD_SYM --> RX_LLR_BIT --> RX_LLR_SYM
		// only support bpsk now
		int p_i = 0;
		for(int b = 0; b < CODE_BIT_LEN; b ++)
		{
			if(PUN_BIT[p_i] == b && PUN_BIT_LEN != 0)
			{
				p_i ++;
				RX_LLR_BIT[b] = 0;
			}
			else
			{
				RX_LLR_BIT[b] = -2 * RX_MOD_SYM[b - p_i].Real / (sigma_n * sigma_n);
			}
		}
		// RX_LLE_BIT --> RX_LLR_SYM
		for(int s = 0; s < CODE_SYM_LEN; s ++)
		{
			for(int q = 1; q < GFq; q ++)
			{
				RX_LLR_SYM[s][q - 1] = 0;
				for(int b_p_s = 0; b_p_s < Bit_Len_PerSYM; b_p_s ++)
				{
					if((q & (1 << b_p_s)) != 0)
					{
						RX_LLR_SYM[s][q - 1] += RX_LLR_BIT[s * Bit_Len_PerSYM + b_p_s];
					}
				}
			}
		}
	}
	
	else if(modOrder == GFq)
	{
		int p_i = 0;
		for(int s = 0; s < CODE_SYM_LEN; s ++)
		{
			if(PUN_SYM[p_i] == s && PUN_SYM_LEN != 0)
			{
				p_i ++;
				for(int q = 1; q < GFq; q ++)
				{
					RX_LLR_SYM[s][q - 1] = 0;
				}
				
			}
			else
			{
				for(int q = 1; q < GFq; q ++)
				{
					RX_LLR_SYM[s][q - 1] = ( (2 * RX_MOD_SYM[s - p_i].Real - CONSTELLATION[0].Real - CONSTELLATION[q].Real ) * (CONSTELLATION[q].Real - CONSTELLATION[0].Real) 
						+ (2 * RX_MOD_SYM[s - p_i].Image - CONSTELLATION[0].Image - CONSTELLATION[q].Image ) * (CONSTELLATION[q].Image - CONSTELLATION[0].Image) ) / (2 * sigma_n * sigma_n);
				}
			}
		}
	}
	else
	{
		cout << "This module ( code order ~= modulation order ) haven't been developed!" << endl;
		exit(-1);

	}
	return 0;
}

//****************************************************************Decode&Err********************************************************
int CComm::Decode(void)
{
	// NBLDPC deocding
	for (int i = 0; i < CODE_BIT_LEN; i++) {
		NBLDPC.TX_code[i] = TX_CODE_BIT[i];
	}

	DecodeCorrect = NBLDPC.Decoding(RX_LLR_SYM, RX_DECODE_SYM, ReliableSeri_Symbol,ReliableSeri_BIT);
	
	tems_route_num = NBLDPC.tems_route_num;
	 
	// RX_DECODE_SYM --> RX_DECODE_BIT
	for(int s = 0; s < CODE_SYM_LEN; s ++)
	{
		for(int b_p_s = 0; b_p_s < Bit_Len_PerSYM; b_p_s ++)
		{
			RX_DECODE_BIT[s * Bit_Len_PerSYM + b_p_s]
			= ((RX_DECODE_SYM[s] & (1 << b_p_s)) == 0)? 0 : 1;
		}
	}

	// RX_DECODE_SYM --> RX_MSG_SYM
	for(int s = 0; s < MSG_SYM_LEN; s ++)
	{
		RX_MSG_SYM[s] = RX_DECODE_SYM[s];
		
	}
	
	// RX_DECODE_BIT --> RX_MSG_BIT
	for(int b = 0; b < MSG_BIT_LEN; b ++)
	{
		RX_MSG_BIT[b] = RX_DECODE_BIT[b];
	}
	return 0;
}

int CComm::Err(CSimulation &sim)
{
	// BER, SER, FER
	double errSym, errBit, errFrame;
	errSym = errBit = errFrame = 0;

	for (int s = 0; s < MSG_SYM_LEN; s++)
	{
		if (TX_MSG_SYM[s] != RX_MSG_SYM[s])
		{
			errSym++;
		}
	}
	int pos;
	for (int b = 0; b < MSG_BIT_LEN; b++)
	{
		//pos = ReliableSeri_BIT[b];
		pos = b;
		if (TX_MSG_BIT[pos] != RX_MSG_BIT[pos])
		{
			errBit++;
		}
	}

	int CRC_ckeck_flag;
	CRC_ckeck_flag = CrcCheck(RX_MSG_BIT, MSG_BIT_LEN, crcLen, 1);


	/*
	cout << endl;
	cout << "通信出错  "<< errBit <<endl;
	cout << "发射端  "  << endl;
	for (int i = 0; i < 128; i++) { cout << TX_CODE_BIT[i]; }
	cout << endl;
	cout << "接收端  " << endl;
	for (int i = 0; i < 128; i++) { cout << RX_DECODE_BIT[i]; }
	cout << endl;
	*/
	sim.errSym += errSym;
	sim.errBit += errBit;
	sim.errFrame += (errSym != 0) ? 1 : 0;

	if (CRC_ckeck_flag && (errSym != 0))
	{
		sim.U_errSym += errSym;
		sim.U_errBit += errBit;
		sim.U_errFrame += 1;
	}

	sim.SER = sim.errSym / (sim.simCycle * MSG_SYM_LEN*sim.parallel);
	sim.BER = sim.errBit / (sim.simCycle * MSG_BIT_LEN*sim.parallel);
	sim.FER = sim.errFrame / (sim.simCycle*sim.parallel);

	sim.U_SER = sim.U_errSym / (sim.simCycle * MSG_SYM_LEN*sim.parallel);
	sim.U_BER = sim.U_errBit / (sim.simCycle * MSG_BIT_LEN*sim.parallel);
	sim.U_FER = sim.U_errFrame / (sim.simCycle*sim.parallel);
	return 0;
}

//**************************************************CRC*********************************************
void CComm::CRCEncode(int *seqOut, int *seqIn, int seqInLen, int crcLen, int crc24Type)
{
	int i, j;
	int temp;
	int reg[24] = { 0 };
	int G[25] = { 0 };

	//-	gCRC24A(D) = [D24 + D23 + D18 + D17 + D14 + D11 + D10 + D7 + D6 + D5 + D4 + D3 + D + 1] and;
	//-	gCRC24B(D) = [D24 + D23 + D6 + D5 + D + 1] for a CRC length L = 24 and;
	//-	gCRC16(D) = [D16 + D12 + D5 + 1] for a CRC length L = 16.
	//-	gCRC8(D) = [D8 + D7 + D4 + D3 + D + 1] for a CRC length of L = 8.
	switch (crcLen)
	{
	case 8: G[0] = G[1] = G[3] = G[4] = G[7] = G[8] = 1; break;
	case 16:G[0] = G[5] = G[12] = G[16] = 1; break;
	case 24:
		if (!crc24Type)
		{
			G[0] = G[1] = G[3] = G[4] = G[5] = G[6] = G[7] = G[10] = G[11] = G[14] = G[17] = G[18] = G[23] = G[24] = 1;
		}
		else
		{
			G[0] = G[1] = G[5] = G[6] = G[23] = G[24] = 1;
		}
		break;
	default: break;
	}

	if (crcLen != 0)
	{
		for (i = 0; i < seqInLen; i++)
		{
			temp = reg[crcLen - 1] ^ seqIn[i];
			for (j = crcLen - 1; j > 0; j--)
			{
				reg[j] = reg[j - 1] ^ (G[j] && temp);
			}
			reg[0] = temp;
		}
		if (seqIn != seqOut)
		{
			memcpy(seqOut, seqIn, seqInLen * sizeof(int));
		}
		for (i = 0; i < crcLen; i++)
		{
			seqOut[i + seqInLen] = reg[crcLen - 1 - i];
		}
	}
	else//if CRC length is 0, just copy input to output
	{
		if (seqIn != seqOut)
		{
			memcpy(seqOut, seqIn, seqInLen * sizeof(int));
		}
	}
}


int CComm::CrcCheck(int *seqIn, int seqInLen, int crcLen, int crc24Type)
{
	int i = 0;
	int j = 0;
	int temp = 0;
	int sum = 0;
	int reg[24] = { 0 };
	int G[25] = { 0 };
	int sumZero = 0;		//used to find out whether input sequence is a all 0 sequence or not

							//-	gCRC24A(D) = [D24 + D23 + D18 + D17 + D14 + D11 + D10 + D7 + D6 + D5 + D4 + D3 + D + 1] and;
							//-	gCRC24B(D) = [D24 + D23 + D6 + D5 + D + 1] for a CRC length L = 24 and;
							//-	gCRC16(D) = [D16 + D12 + D5 + 1] for a CRC length L = 16.
							//-	gCRC8(D) = [D8 + D7 + D4 + D3 + D + 1] for a CRC length of L = 8.
	switch (crcLen)
	{
	case 8: G[0] = G[1] = G[3] = G[4] = G[7] = G[8] = 1; break;
	case 16:G[0] = G[5] = G[12] = G[16] = 1; break;
	case 24:
		if (!crc24Type)
		{
			G[0] = G[1] = G[3] = G[4] = G[5] = G[6] = G[7] = G[10] = G[11] = G[14] = G[17] = G[18] = G[23] = G[24] = 1;
		}
		else
		{
			G[0] = G[1] = G[5] = G[6] = G[23] = G[24] = 1;
		}
		break;
	default: break;
	}


	if (crcLen != 0)
	{
		for (i = 0; i < seqInLen; i++)
		{
			temp = reg[crcLen - 1];
			for (j = crcLen - 1; j >0; j--)
			{
				reg[j] = reg[j - 1] ^ (G[j] && temp);
			}
			reg[0] = temp ^ seqIn[i];
			sumZero += seqIn[i];
		}

		sum = 0;
		for (i = 0; i < crcLen; i++)
		{
			sum += reg[i];
		}

		if (sum != 0)
		{
			return 0;
		}
		else
		{
			if (0 != sumZero)// can't be all zero sqence
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
	}
	else
	{
		return 1;
	}

}




