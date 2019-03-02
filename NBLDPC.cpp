#include "NBLDPC.h"
#include "OSD.h"


CNBLDPC::CNBLDPC(void)
	: GFq(0)
	, CodeLen(0)
	, ChkLen(0)
	, PuncturePosition(nullptr)
	, PunctureLen(0)
	, maxVarDegree(0)
	, maxChkDegree(0)
	, VarDegree(nullptr)
	, ChkDegree(nullptr)
	, VarLink(nullptr)
	, ChkLink(nullptr)
	, VarLinkGFe(NULL)
	, ChkLinkGFe(NULL)
	, maxIter(0)
	, L_v2c(NULL)
	, L_c2v(NULL)
	, LLRV(NULL)
	, L_post(NULL)
	, VarLinkDc(NULL)
	, ChkLinkDv(NULL)
	, L_sigma(NULL)
	, L_rho(NULL)
	, EncodeExchangeColSource(NULL)
	, EncodeExchangeColDest(NULL)
	, EncodeExchangeColNum(0)
	, EncodeVarLink(NULL)
	, EncodeVarLinkGFe(NULL)
	, EncodeVarDv(NULL)
	, EMS_Nm(0)
	, EMS_Nc(0)
	, DecodeMethod(0)
	, EMS_sort_L_v2c(NULL)
	, EMS_sort_Entr_v2c(NULL)
	, EMS_L_c2v(NULL)
	, EMS_Correction_Factor(0)
	, EMS_Correction_Offset(0)
	, TEMS_Nr(0)
	, TEMS_Nc(0)
	, TEMS_deltaU(NULL)
	, TEMS_deltaW(NULL)
	, TEMS_Beta(NULL)
	, TEMS_Eta(NULL)
	, TEMS_Min(NULL)
	, TEMS_inConf_Nr(NULL)
	, TEMS_Eta_Candidate(NULL)
	, TEMS_Lc2v(NULL)
	, TEMS_isUpdate(NULL)
	, TEMS_Syndrome(0)
	, TEMS_Factor(0)
	, TEMS_Offset(0)
	, TEMS_isSelected(NULL)
{
}


CNBLDPC::~CNBLDPC(void)
{
}

//外界辅助代码
void Symol_LLR2P(double* LLR_S, double* Prob_S, int Q)
{
	int q;
	double temp_sum = 0;

	for (q = 0; q < Q - 1; q++)
	{
		temp_sum += exp(LLR_S[q]);
	}

	Prob_S[0] = 1 / (1 + temp_sum);

	for (q = 1; q < Q; q++)
	{
		Prob_S[q] = exp(LLR_S[q - 1])*Prob_S[0];
	}

}

void Prob_S2LLR_b(double* Prob_S, int Q, double* LLR_b, int L)
{
	int q, l;

	double* Bit_Prob0 = new double[L];
	double* Bit_Prob1 = new double[L];

	for (l = 0; l < L; l++)
	{
		Bit_Prob0[l] = 0;
		Bit_Prob1[l] = 0;
	}

	for (q = 0; q < Q; q++)
	{
		for (l = 0; l < L; l++)
		{
			if ((q & (1 << (L - 1 - l))) == 0)
			{
				Bit_Prob0[l] += Prob_S[q];
			}
			else
			{
				Bit_Prob1[l] += Prob_S[q];
			}
		}
	}
	for (l = 0; l < L; l++)
	{
		LLR_b[l] = log(Bit_Prob1[l] / Bit_Prob0[l]);
	}

	delete[] Bit_Prob0;
	delete[] Bit_Prob1;
}


void LLR_b2S(double* LLR_b, int L, double* LLR_S, int Q)
{
	int l;

	for (int q = 1; q < Q; q++)
	{
		LLR_S[q - 1] = 0;
		for (l = 0; l < L; l++)
		{
			if ((q & (1 << (L - 1 - l))) != 0)
			{
				LLR_S[q - 1] += LLR_b[l];
			}
		}
	}
}


bool CNBLDPC::Initial(CSimulation &sim)
{
	GFq = sim.GFq;
	maxIter = sim.maxIter;
	// initial GF
	GF.Initial(GFq);

	ifstream fnb(sim.NonBinaryFileName);
	fnb >> CodeLen >> ChkLen >> GFq;
	fnb >> maxVarDegree >> maxChkDegree;
	VarDegree = new int[CodeLen];
	ChkDegree = new int[ChkLen];

	PunctureLen = 0;
	VarLink = new int*[CodeLen];
	VarLinkGFe = new int*[CodeLen];
	VarLinkDc = new int*[CodeLen];
	for(int col = 0; col < CodeLen; col ++)
	{
		fnb >> VarDegree[col];
		VarLink[col] = new int[VarDegree[col]];
		VarLinkGFe[col] = new int[VarDegree[col]];
		VarLinkDc[col] = new int[VarDegree[col]];
		if(VarDegree[col] == sim.PuntureVarDegree)
		{
			PunctureLen ++;
		}
	}
	PuncturePosition = new int[PunctureLen];
	int PunctureIndex = 0;
	for(int col = 0; col < CodeLen; col ++)
	{
		if(VarDegree[col] == sim.PuntureVarDegree)
		{
			PuncturePosition[PunctureIndex] = col;
			PunctureIndex ++;
		}
	}
	ChkLink = new int*[ChkLen];
	ChkLinkGFe = new int*[ChkLen];
	ChkLinkDv = new int*[ChkLen];
	for(int row = 0; row < ChkLen; row ++)
	{
		fnb >> ChkDegree[row];
		ChkLink[row] = new int[ChkDegree[row]];
		ChkLinkGFe[row] = new int[ChkDegree[row]];
		ChkLinkDv[row] = new int[ChkDegree[row]];
	}

	for(int col = 0; col < CodeLen; col ++)
	{
		for(int d = 0; d < VarDegree[col]; d ++)
		{
			fnb >> VarLink[col][d] >> VarLinkGFe[col][d];
			VarLink[col][d] --;
		}
	}
	for(int row = 0; row < ChkLen; row ++)
	{
		for(int d = 0; d < ChkDegree[row]; d ++)
		{
			fnb >> ChkLink[row][d] >> ChkLinkGFe[row][d];
			ChkLink[row][d] --;
		}
	}
	fnb.close();
	
	
	
	// allocate memory
	LLRV = new double[GFq - 1];
	L_post = new double *[CodeLen];
	for(int col = 0; col < CodeLen; col ++)
	{
		L_post[col] = new double[GFq - 1];
	}
	L_v2c = new double **[CodeLen];
	for(int col = 0; col < CodeLen; col ++)
	{
		L_v2c[col] = new double *[VarDegree[col]];
		for(int d = 0; d < VarDegree[col]; d ++)
		{
			L_v2c[col][d] = new double[GFq - 1];
		}
	}
	L_c2v = new double **[ChkLen];
	for(int row = 0; row < ChkLen; row ++)
	{
		L_c2v[row] = new double* [ChkDegree[row]];
		for(int d = 0; d < ChkDegree[row]; d ++)
		{
			L_c2v[row][d] = new double[GFq - 1];
		}
	}

	// create the link of VarLinkDc & ChkLinkDv
	for(int col = 0; col < CodeLen; col ++)
	{
		for(int dv = 0; dv < VarDegree[col]; dv ++)
		{
			int row = VarLink[col][dv];
			for(int dc = 0; dc < ChkDegree[row]; dc++)
			{
				if(ChkLink[row][dc] == col)
				{
					VarLinkDc[col][dv] = dc;
				}
			}
		}
	}
	for(int row = 0; row < ChkLen; row ++)
	{
		for(int dc = 0; dc < ChkDegree[row]; dc ++)
		{
			int col = ChkLink[row][dc];
			for(int dv = 0; dv < VarDegree[col]; dv ++)
			{
				if(VarLink[col][dv] == row)
				{
					ChkLinkDv[row][dc] = dv;
				}
			}
		}
	}


	// allocate memory for iterative decoding
	DecodeMethod = sim.decodeMethod;
	if(DecodeMethod == BP_DECODE)
	{
		L_sigma = new double[GFq - 1];
		L_rho = new double[GFq - 1];
	}

	//EMS：变量赋值与内存分配
	else if(DecodeMethod == EMS_DECODE)
	{
		EMS_Nm = sim.ems_nm;
		EMS_Nc = sim.ems_nc;
		EMS_Correction_Factor = sim.ems_factor;
		EMS_Correction_Offset = sim.ems_offset;

		if(EMS_Nm > GFq)
		{
			cerr << "EMS configuration error! EMS_Nm is too large!" << endl;
			exit(-1);
		}

		EMS_sort_L_v2c = new double *[maxChkDegree]; 
		EMS_sort_Entr_v2c = new int*[maxChkDegree];
		for(int d = 0; d < maxChkDegree; d ++)
		{
			EMS_sort_L_v2c[d] = new double[GFq];
			EMS_sort_Entr_v2c[d] = new int[GFq];
		}
		EMS_L_c2v = new double[GFq];
	}

	//T-EMS：变量赋值与分配内存
	else if(DecodeMethod == T_EMS_DECODE)
	{
		TEMS_Nr = sim.tems_nr;
		TEMS_Nc = sim.tems_nc;
		TEMS_Factor = sim.tems_factor;
		TEMS_Offset = sim.tems_offset;

		TEMS_deltaU = new double*[maxChkDegree];
		for(int d = 0; d < maxChkDegree; d ++)
		{
			TEMS_deltaU[d] = new double[GFq];
		}
		TEMS_Beta = new int[maxChkDegree];
		TEMS_Min = new int*[GFq];
		TEMS_inConf_Nr = new bool*[GFq];
		for(int q = 0; q < GFq; q ++)
		{
			TEMS_Min[q] = new int[maxChkDegree];
			TEMS_inConf_Nr[q] = new bool[maxChkDegree];
		}
		TEMS_deltaW = new double[GFq];
		TEMS_Eta = new int*[GFq];
		for(int q = 0; q < GFq; q ++)
		{
			TEMS_Eta[q] = new int[maxChkDegree];
		}
		TEMS_Eta_Candidate = new int[maxChkDegree];
		TEMS_Lc2v = new double[GFq];
		TEMS_isUpdate = new bool[GFq];
		TEMS_isSelected = new bool[GFq];
	}
	
	//****************************************Basic_Set*****************************************************************
	else if (DecodeMethod == BS_TEMS_DECODE) {
		length_basicSet = GFq;
		basic_set_Nc = sim.bs_tems_nc;
		basic_set_Nm = sim.bs_tems_nm;
		basic_set_factor = sim.bs_tems_factor;
		basic_set_offset = sim.bs_tems_offset;

		TEMS_deltaU = new double*[maxChkDegree];
		for (int d = 0; d < maxChkDegree; d++)
		{
			TEMS_deltaU[d] = new double[GFq];
		}
		TEMS_Beta = new int[maxChkDegree];
		TEMS_Min = new int*[GFq];
		TEMS_inConf_Nr = new bool*[GFq];
		for (int q = 0; q < GFq; q++)
		{
			TEMS_Min[q] = new int[maxChkDegree];
			TEMS_inConf_Nr[q] = new bool[maxChkDegree];
		}
		TEMS_deltaW = new double[GFq];
		TEMS_Eta = new int*[GFq];
		for (int q = 0; q < GFq; q++)
		{
			TEMS_Eta[q] = new int[maxChkDegree];
		}
		TEMS_Eta_Candidate = new int[maxChkDegree];
		TEMS_Lc2v = new double[GFq];
		TEMS_isUpdate = new bool[GFq];
		TEMS_isSelected = new bool[GFq];

		basicSet_col = new int[length_basicSet];
		basicSet_q = new int[length_basicSet];
		basicSet_L = new double[length_basicSet];

		Set_col = new int[GFq];
		Set_q = new int[GFq];
		basicSet_available = new bool[GFq];
		col_isSelected = new bool[maxChkDegree];
		TEMS_sort_deltaU = new double*[maxChkDegree];
		TEMS_sort_q = new int*[maxChkDegree];
		for (int i = 0; i < maxChkDegree; i++) {
			TEMS_sort_deltaU[i] = new double[GFq];
			TEMS_sort_q[i] = new int[GFq];
		}
	}
	

	//***********************************************OSD******************************************
	OSD_Alpha = sim.OSD_factor;
	OSD_Order = sim.OSD_order;
	OSD_flag = sim.OSD_flag;
	CodeLen_bit = CodeLen*(log(GFq) / log(2));
	ChkLen_bit = ChkLen*(log(GFq) / log(2));

	L_ch_bit = new double[CodeLen_bit];//解调软输出的比特结
	LLR_Sym = new double**[maxIter];//迭代算法每一次迭代的软输出
	for (int i = 0; i < maxIter; i++) {
		LLR_Sym[i] = new double*[CodeLen];
		for (int j = 0; j < CodeLen; j++) {
			LLR_Sym[i][j] = new double[GFq - 1];
		}
	}
	LLR_Sym_sum = new double*[CodeLen];//迭代算法的累积软输出
	for (int i = 0; i < CodeLen; i++) {
		LLR_Sym_sum[i] = new double[GFq - 1];
	}

	infor_record = new int[CodeLen];//参考信息，选取解调软输出的硬判结果或者迭代算法的累积软输出
	infor_record_bit = new int[CodeLen_bit];
	ReliableSeri_BIT = new int[CodeLen_bit];//参考信息依据的软量
	OSD_temp_out_bit = new int[CodeLen_bit];//每一次翻转后的OSD输出
	near_optimal = new int[CodeLen];//保存欧式距离最小的OSD输出
	near_optimal_bit = new int[CodeLen_bit];//

		
	
	H_origin_bit = new int*[ChkLen_bit];
	//H_OSD_bit = new int*[ChkLen_bit];
	for (int i = 0; i < ChkLen_bit; i++) {
		H_origin_bit[i] = new int[CodeLen_bit];
		for (int j = 0; j < CodeLen_bit; j++) { H_origin_bit[i][j] = 0; }
		//H_OSD_bit[i] = new int[CodeLen_bit];
	}
	
	int P = log(GFq) / log(2);
	int element_temp;
	int chk_pos_temp;
	for (int i = 0; i < CodeLen; i++) {
		for (int j = 0; j < VarDegree[i]; j++) {
			chk_pos_temp = VarLink[i][j];
			element_temp = VarLinkGFe[i][j];
			for (int k = 0; k < P; k++) {
				for (int l = 0; l < P; l++) {
					H_origin_bit[P*chk_pos_temp + k][P*i + l] = GF.GFElement[element_temp].ValueMatric[l][k];
				}
			}
		}	
	}	

	int MsgLen_bit = CodeLen_bit - ChkLen_bit;
	G_OSD_bit = new int*[MsgLen_bit];
	for (int i = 0; i < MsgLen_bit; i++) {
		G_OSD_bit[i] = new int[CodeLen_bit];
	}
	MatrixH2G(H_origin_bit, G_OSD_bit,ChkLen_bit,CodeLen_bit);
	

	//**************************************CRC****************************************
	CRCLen = sim.crcLen;
	CRCrate = sim.crc_correctLen;
	CRCMatrix_G = new int*[MsgLen_bit-CRCLen];
	CRCMatrix_G_IP = new int*[MsgLen_bit - CRCLen];
	CRCMatrix_partH_PI = new int*[CRCrate];
	H_OSD_CRC_bit = new int*[ChkLen_bit + CRCrate];
	for (int i = 0; i < MsgLen_bit - CRCLen; i++) {
		CRCMatrix_G[i] = new int[MsgLen_bit]();
		CRCMatrix_G_IP[i] = new int[MsgLen_bit]();
	}
	for (int i = 0; i < CRCrate; i++) {
		CRCMatrix_partH_PI[i] = new int[MsgLen_bit]();
	}
	for (int i = 0; i < ChkLen_bit + CRCrate; i++) {
		H_OSD_CRC_bit[i] = new int[CodeLen_bit]();
	}
	CRCMatrixGen();

	//**************************************test_Conf********************************************************************
	TEMS_sort_deltaU = new double*[maxChkDegree];
	TEMS_sort_q = new int*[maxChkDegree];
	for (int i = 0; i < maxChkDegree; i++) {
		TEMS_sort_deltaU[i] = new double[GFq];
		TEMS_sort_q[i] = new int[GFq];
	}

	
	// initial encode
	if (sim.randomMsg)
	{
		InitialEncode();
	}
	
	return true;
}

int CNBLDPC::InitialEncode(void)
{

	// allocate memory for encode exchange
	EncodeExchangeColSource = new int[CodeLen];
	EncodeExchangeColDest = new int[CodeLen];
	EncodeExchangeColNum = 0;
	// allocate matrix H
	int ** H = new int*[ChkLen];
	for (int row = 0; row < ChkLen; row++)
	{
		H[row] = new int[CodeLen];
	}
	// reset value of H
	for (int row = 0; row < ChkLen; row++)
	{
		for (int col = 0; col < CodeLen; col++)
		{
			H[row][col] = 0;
		}
	}
	for (int row = 0; row < ChkLen; row++)
	{
		for (int d = 0; d < ChkDegree[row]; d++)
		{
			H[row][ChkLink[row][d]] = ChkLinkGFe[row][d];
		}
	}
	/*
	cout << "H" << endl;
	for (int i = 0; i < ChkLen; i++) {
		for (int j = 0; j < CodeLen; j++) {
			cout << H[i][j] << "\t";
		}
		cout << endl;
	}
	*/

	GaussEliminate(H, EncodeExchangeColSource, EncodeExchangeColDest, EncodeExchangeColNum);
	
	// allocate memory for encode var node
	EncodeVarDv = new int[ChkLen];
	EncodeVarLink = new int*[ChkLen];
	EncodeVarLinkGFe = new int*[ChkLen];

	int all = 0;
	for (int p = 0; p < ChkLen; p++)
	{
		EncodeVarDv[p] = 0;
		for (int col = 0; col < CodeLen - ChkLen + p; col++)
		{
			if (H[p][col] != 0)
			{
				EncodeVarDv[p] ++;
			}
		}
		all += EncodeVarDv[p];
		EncodeVarLink[p] = new int[EncodeVarDv[p]];
		EncodeVarLinkGFe[p] = new int[EncodeVarDv[p]];
	}
	for (int p = 0; p < ChkLen; p++)
	{
		int d = 0;
		for (int col = 0; col < CodeLen - ChkLen + p; col++)
		{
			if (H[p][col] != 0)
			{
				EncodeVarLink[p][d] = col;
				EncodeVarLinkGFe[p][d] = H[p][col];
				d++;
			}
		}
	}

	// release memory of H
	
	for (int row = 0; row < ChkLen; row++)
	{
		delete[] H[row];
	}
	delete[] H;
	
	return 0;
}

int CNBLDPC::Encode(int* msg_sym, int* code_sym)
{
	// put information symbol
	for (int col = 0; col < CodeLen - ChkLen; col++)
	{
		code_sym[col] = msg_sym[col];
	}
	for (int col = CodeLen - ChkLen; col < CodeLen; col++)
	{
		code_sym[col] = 0;
	}

	// encode
	for (int p = 0; p < ChkLen; p++)
	{
		int col = CodeLen - ChkLen + p;
		int v, h, x;
		for (int d = 0; d < EncodeVarDv[p]; d++)
		{
			v = EncodeVarLink[p][d];
			h = EncodeVarLinkGFe[p][d];
			x = code_sym[v];
			code_sym[col] = GF.GFAdd(code_sym[col], GF.GFMultiply(h, x));
		}
	}

	// permute
	int temp, col_s, col_d;
	for (int num = EncodeExchangeColNum - 1; num >= 0; num--)
	{
		col_s = EncodeExchangeColSource[num];
		col_d = EncodeExchangeColDest[num];
		temp = code_sym[col_s];
		code_sym[col_s] = code_sym[col_d];
		code_sym[col_d] = temp;
	}
	for (int col = 0; col < CodeLen - ChkLen; col++)
	{
		msg_sym[col] = code_sym[col];
	}

	return 0;
}


int CNBLDPC::Decoding(double** L_ch, int* DecodeOutput,int* RelySeri_symbol,int* RelySeri_bit)
{
	int decoderesult;
	switch (DecodeMethod)
	{
	case BP_DECODE:
		decoderesult = Decoding_BP(L_ch, DecodeOutput);
		break;
	case EMS_DECODE:
		decoderesult = Decoding_EMS(L_ch, DecodeOutput, RelySeri_bit);
		break;
	case MinMax_DECODE:
		cout << "Min Max algorithm has not been developed" << endl;
		exit(-1);
		break;
	case T_EMS_DECODE:	
		decoderesult =Decoding_TEMS(L_ch, DecodeOutput,RelySeri_symbol, RelySeri_bit);
		break;
	case OSD_DECODE:
		decoderesult = Decoding_OSD_bit(L_ch, DecodeOutput, RelySeri_bit,1);
		break;
	case BS_TEMS_DECODE:
		decoderesult = Decoding_BS_TEMS(L_ch, DecodeOutput, RelySeri_symbol, RelySeri_bit);
		break;
	case T_MinMax_DECODE:
		cout << "Trellis Min Max algorithm has not been developed" << endl;
		exit(-1);
		break;

	default:
		cout << "Decode Method is not specified!" << endl;
		exit(-1);
	}
	return decoderesult;
}

int CNBLDPC::Decoding_BP(double** L_ch, int* DecodeOutput)
{
//	ofstream decoderecord("decode.txt");
	// initial
	double* old_L_v2c = new double[GFq - 1];

	for(int col = 0; col < CodeLen; col ++)
	{
		for(int d = 0; d < VarDegree[col]; d ++)
		{
			CopyLLRVector(L_v2c[col][d], L_ch[col]);
		}
	}
	for(int row = 0; row < ChkLen; row ++)
	{
		for(int d = 0; d < ChkDegree[row]; d ++)
		{
			ClearLLRVector(L_c2v[row][d]);
		}
	}
	for (int i = 0; i < CodeLen; i++) {
		for (int j = 0; j < GFq - 1; j++) {
			LLR_Sym_sum[i][j] = 0;
		}
	}


	// iterative decoding
	int iter_number = 0;
	bool decode_correct = true;
	while(iter_number ++ < maxIter)
	{
		// tentative decoding
		for(int col = 0; col < CodeLen; col ++)
		{
			CopyLLRVector(L_post[col], L_ch[col]);
			int row, dc;
			for(int d = 0; d < VarDegree[col]; d ++)
			{
				row = VarLink[col][d];
				dc = VarLinkDc[col][d];
				AddLLRVector(L_post[col], L_post[col], L_c2v[row][dc]);
			}
			for (int i = 0; i < GFq - 1; i++) {
				LLR_Sym_sum[col][i] = OSD_Alpha*LLR_Sym_sum[col][i] + L_post[col][i];
			}
			DecodeOutput[col] = DecideLLRVector(L_post[col]);
	//		decoderecord << DecodeOutput[col] << " " << L_post[col][DecodeOutput[col] - 1] << " ";
		}
		// check whether satisfy the parity matrix
		decode_correct = true;
		for(int row = 0; row < ChkLen; row ++)
		{
			int col, multip, syndrome;
			syndrome = 0;
			for(int d = 0; d < ChkDegree[row]; d ++)
			{
				col = ChkLink[row][d];
				multip = GF.GFMultiply(ChkLinkGFe[row][d], DecodeOutput[col]);
				syndrome = GF.GFAdd(syndrome, multip);
			}
			if(syndrome != 0)
			{
				decode_correct = false;
				break;
			}
		}
		if(decode_correct)
		{
			delete[] old_L_v2c;

			return decode_correct;
		}

		// message from var to check
		for(int col = 0; col < CodeLen; col ++)
		{
			int row, dc;

			int update_before = -1;
			int update_after = -1;

			for(int dv = 0; dv < VarDegree[col]; dv ++)
			{
				row = VarLink[col][dv];
				dc = VarLinkDc[col][dv];

				update_before = DecideLLRVector(L_v2c[col][dv]);
				CopyLLRVector(old_L_v2c, L_v2c[col][dv]);
				MinusLLRVector(L_v2c[col][dv], L_post[col], L_c2v[row][dc]);

				update_after = DecideLLRVector(L_v2c[col][dv]);
				if (update_after != update_before)
				{
					for (int i = 0; i < GFq - 1; i++)
					{
						L_v2c[col][dv][i] = 0.5*old_L_v2c[i] + 0.5*L_v2c[col][dv][i];
					}
				}

			}
		}

		// message from check to var
		for(int row = 0; row < ChkLen; row ++)
		{
			for(int d = 0; d < ChkDegree[row]; d ++)
			{
				ClearLLRVector(L_sigma);
				ClearLLRVector(L_rho);
				L_Back(row, d - 1);
				L_Forward(row, d + 1);
				int A1, A2;
				A1 = A2 = GF.GFInverse(ChkLinkGFe[row][d]);
				if(d == 0)
				{
					A1 = 0;
				}
				else if(d == ChkDegree[row] - 1)
				{
					A2 = 0;
				}
				LLR_BoxPlus(L_c2v[row][d], L_sigma, L_rho, A1, A2);
			}
		}
	}
	//OSD部分,
	if (OSD_Order >= 0) {
		if (decode_correct == false) {
			Decoding_OSD_bit(L_ch, DecodeOutput, ReliableSeri_BIT, OSD_flag);
		}

	}
}

int CNBLDPC::Decoding_EMS(double** L_ch, int* DecodeOutput,int*RelySeri_bit)
{
	// initial
	for(int col = 0; col < CodeLen; col ++)
	{
		for(int d = 0; d < VarDegree[col]; d ++)
		{
			CopyLLRVector(L_v2c[col][d], L_ch[col]);
		}
	}
	for(int row = 0; row < ChkLen; row ++)
	{
		for(int d = 0; d < ChkDegree[row]; d ++)
		{
			ClearLLRVector(L_c2v[row][d]);
		}
	}

	for (int i = 0; i < CodeLen; i++) {
		for (int j = 0; j < GFq - 1; j++) {
			LLR_Sym_sum[i][j] = 0;
		}
	}

	// iterative decoding
	int iter_number = 0;
	bool decode_correct = true;
	while(iter_number ++ < maxIter)
	{
		// tentative decoding
		for(int col = 0; col < CodeLen; col ++)
		{
			CopyLLRVector(L_post[col], L_ch[col]);
			int row, dc;
			for(int d = 0; d < VarDegree[col]; d ++)
			{
				row = VarLink[col][d];
				dc = VarLinkDc[col][d];
				AddLLRVector(L_post[col], L_post[col], L_c2v[row][dc]);
			}
			//每个变量节点的累积可信度计算
			for (int i = 0; i < GFq - 1; i++) {
				LLR_Sym_sum[col][i] = OSD_Alpha*LLR_Sym_sum[col][i] + L_post[col][i];
			}
			DecodeOutput[col] = DecideLLRVector(L_post[col]);
		}

		// check whether satisfy the parity matrix
		decode_correct = true;
		for(int row = 0; row < ChkLen; row ++)
		{
			int col, multip, syndrome;
			syndrome = 0;
			for(int d = 0; d < ChkDegree[row]; d ++)
			{
				col = ChkLink[row][d];
				multip = GF.GFMultiply(ChkLinkGFe[row][d], DecodeOutput[col]);
				syndrome = GF.GFAdd(syndrome, multip);
			}
			if(syndrome != 0)
			{
				decode_correct = false;
				break;
			}
		}
		if(decode_correct)
		{
			return 1;
		}
		// message from var to check
		for(int col = 0; col < CodeLen; col ++)
		{
			int row, dc;
			for(int dv = 0; dv < VarDegree[col]; dv ++)
			{
				row = VarLink[col][dv];
				dc = VarLinkDc[col][dv];
				MinusLLRVector(L_v2c[col][dv], L_post[col], L_c2v[row][dc]);
			}
		}
		// message from check to var
		for(int row = 0; row < ChkLen; row ++)
		{
			// sort to get the Nm maximum LLR
			int col, dv;
			for(int dc = 0; dc < ChkDegree[row]; dc ++)
			{
				col = ChkLink[row][dc];
				dv = ChkLinkDv[row][dc];
				SortLLRVector(EMS_sort_L_v2c[dc], EMS_sort_Entr_v2c[dc], GFq, L_v2c[col][dv], GFq - 1);
			}
			// print
			/*			
			ofstream fsort("sort.txt");
			for(int dc = 0; dc < ChkDegree[row]; dc ++)
			{
				for(int i = 0; i < GFq; i ++)
				{
					fsort << EMS_sort_L_v2c[dc][i] << " " << EMS_sort_Entr_v2c[dc][i] << " ";
				}
				fsort << endl;
			}
			fsort.close();
			*/
			for(int dc = 0; dc < ChkDegree[row]; dc ++)
			{
				// reset the sum store vector to the munimum
				for(int q = 0; q < GFq; q ++)
				{
					EMS_L_c2v[q] = -DBL_MAX;
				}
				// recursly exhaustly
				int sumNonele, diff;
				double sumNonLLR;
				// conf(q, 1)
				sumNonele = 0; sumNonLLR = 0; diff = 0;
				ConstructConf(GFq, 1, sumNonele, sumNonLLR, diff, 0, dc, ChkDegree[row] - 1, row);
				// conf(nm, nc)
				sumNonele = 0; sumNonLLR = 0; diff = 0;
				ConstructConf(EMS_Nm, EMS_Nc, sumNonele, sumNonLLR, diff, 0, dc, ChkDegree[row] - 1, row);
				// calculate each c2v LLR
				int v = 0;
				for(int k = 1; k < GFq; k ++)
				{
					v = GF.GFMultiply(k, ChkLinkGFe[row][dc]);
					L_c2v[row][dc][k - 1] = (EMS_L_c2v[v] - EMS_L_c2v[0]) / EMS_Correction_Factor;
					if(L_c2v[row][dc][k - 1] < -1 * EMS_Correction_Offset)
					{
						L_c2v[row][dc][k - 1] = L_c2v[row][dc][k - 1] + EMS_Correction_Offset;
					}
					else if(L_c2v[row][dc][k - 1] >  EMS_Correction_Offset)
					{
						L_c2v[row][dc][k - 1] = L_c2v[row][dc][k - 1] - EMS_Correction_Offset;
					}
					else
					{
						L_c2v[row][dc][k - 1] = 0;
					}
				}
			}
		}
	}
	//OSD部分,
	if (OSD_Order >= 0) {
		if (decode_correct == false) {
			Decoding_OSD_bit(L_ch, DecodeOutput, RelySeri_bit, OSD_flag);
		}
	}
	return 0;
}

int CNBLDPC::Decoding_TEMS(double** L_ch, int* DecodeOutput,int* RelySeri_symbol,int* RelySeri_bit)
{

	// initial
	tems_route_num = 0;
	for (int col = 0; col < CodeLen; col++)
	{
		for (int d = 0; d < VarDegree[col]; d++)
		{
			CopyLLRVector(L_v2c[col][d], L_ch[col]);
		}
	}
	for (int i = 0; i < maxIter; i++) {
		for (int j = 0; j < CodeLen; j++) {
			for (int k = 0; k < GFq - 1; k++) {
				LLR_Sym[i][j][k] = 0;
			}
		}
	}

	int Bit_Len = log(GFq * 1.0) / log(2.0);

	double* old_L_v2c = new double[GFq - 1];
	double* Prob_Symbl = new double[GFq];
	double* LLR_Bit = new double[Bit_Len];
	double* old_Prob_Symbl = new double[GFq];
	double* old_LLR_Bit = new double[Bit_Len];
	
	for (int i = 0; i < CodeLen; i++) {
		for (int j = 0; j < GFq - 1; j++) {
			LLR_Sym_sum[i][j] = 0;
		}	
	}

	for (int row = 0; row < ChkLen; row++)
	{
		for (int d = 0; d < ChkDegree[row]; d++)
		{
			ClearLLRVector(L_c2v[row][d]);
		}
	}
	// iterative decoding
	int iter_number = 0;
	bool decode_correct = true;
	while (iter_number++ < maxIter)
	{
		//		cout << iter_number << " ";
		// tentative decoding
		for (int col = 0; col < CodeLen; col++)
		{
			CopyLLRVector(L_post[col], L_ch[col]);
			int row, dc;
			for (int d = 0; d < VarDegree[col]; d++)
			{
				row = VarLink[col][d];
				dc = VarLinkDc[col][d];
				AddLLRVector(L_post[col], L_post[col], L_c2v[row][dc]);
			}
			//OSD 每个符号的可信度计算
			for (int i = 0; i < GFq-1; i++) {
				LLR_Sym_sum[col][i] = OSD_Alpha*LLR_Sym_sum[col][i] + L_post[col][i];
			}
			DecodeOutput[col] = DecideLLRVector(L_post[col]);
		}
		//保存每次迭代计算后的输出软量
		for (int d = 0; d < CodeLen; d++) {
			CopyLLRVector(LLR_Sym[iter_number-1][d], L_post[d]);
		}

		// check whether satisfy the parity matrix
		decode_correct = true;
		for (int row = 0; row < ChkLen; row++)
		{
			int col, multip, syndrome;
			syndrome = 0;
			for (int d = 0; d < ChkDegree[row]; d++)
			{
				col = ChkLink[row][d];
				multip = GF.GFMultiply(ChkLinkGFe[row][d], DecodeOutput[col]);
				syndrome = GF.GFAdd(syndrome, multip);
			}
			if (syndrome != 0)
			{
				decode_correct = false;
				break;
			}
		}
		if (decode_correct)
		{
			delete[] old_L_v2c;
			delete[] Prob_Symbl;
			delete[] LLR_Bit;
			delete[] old_Prob_Symbl;
			delete[] old_LLR_Bit;			
			OSD_permute(RelySeri_bit);
			tems_route_num = tems_route_num / iter_number;
			return 1;
		}
		// message from var to check

		for (int col = 0; col < CodeLen; col++)
		{
			int row, dc;
			int update_before = -1;
			int update_after = -1;
			for (int dv = 0; dv < VarDegree[col]; dv++)
			{
				row = VarLink[col][dv];
				dc = VarLinkDc[col][dv];
				update_before = DecideLLRVector(L_v2c[col][dv]);
				CopyLLRVector(old_L_v2c, L_v2c[col][dv]);
				MinusLLRVector(L_v2c[col][dv], L_post[col], L_c2v[row][dc]);
				update_after = DecideLLRVector(L_v2c[col][dv]);
				if (update_after != update_before)
				{
					for (int i = 0; i < GFq - 1; i++)
					{
						L_v2c[col][dv][i] = 0.25*old_L_v2c[i] + 0.75*L_v2c[col][dv][i];
					}


				}
			}
		}
		//cout << iter_number << " ";
		// message from check to var
		for (int row = 0; row < ChkLen; row++)
		{
			// get the most reliable symbol beta
			TEMS_Get_Beta(row);
			// transform to the delta domain
			TEMS_Get_deltaU(row);
			// get the two minimum value of each row and mark the Nr minimum pos
			TEMS_Get_Min(row);
			// mark the conf
			for (int q = 0; q < GFq; q++)
			{
				TEMS_deltaW[q] = DBL_MAX;
				TEMS_isSelected[q] = false;
			}
			int sumNonele, diff;
			double sumNonLLR;
			sumNonele = diff = 0; sumNonLLR = 0.0;
			TEMS_ConstructConf(TEMS_Nr, TEMS_Nc, sumNonele, sumNonLLR, diff, 0, ChkDegree[row] - 1);//通过遍历的方式获得TEMS_deltaW[q]以及其对应的路径
			
			//calculate the output
			for (int dc = 0; dc < ChkDegree[row]; dc++)
			{
				// clear the temp vector TEMS_Lc2v
				for (int q = 0; q < GFq; q++)
				{
					TEMS_Lc2v[q] = DBL_MAX;
					TEMS_isUpdate[q] = false;
				}
				// calculate the dc output
				for (int eta = 0; eta < GFq; eta++)
				{
					int deviation = TEMS_Eta[eta][dc];
					int eta_plus_deviation = GF.GFAdd(eta, deviation);
					double minus = TEMS_deltaW[eta] - TEMS_deltaU[dc][deviation];
					if (TEMS_Lc2v[eta_plus_deviation] > minus)
					{
						TEMS_Lc2v[eta_plus_deviation] = minus;
						TEMS_isUpdate[eta_plus_deviation] = true;
					}
				}
				for (int q = 0; q < GFq; q++)
				{
					if (!TEMS_isUpdate[q])
					{
						TEMS_Lc2v[q] = (dc == TEMS_Min[q][0]) ?
							TEMS_deltaU[TEMS_Min[q][1]][q] : TEMS_deltaU[TEMS_Min[q][0]][q];//补足的这一句很重要
					}
				}

				// delta domain --> normal domain & repermutation
				int h_inverse = GF.GFInverse(ChkLinkGFe[row][dc]);
				int beta_syn = GF.GFAdd(TEMS_Syndrome, TEMS_Beta[dc]);
				double L0 = -1.0 * TEMS_Lc2v[beta_syn];
				for (int eta = 0; eta < GFq; eta++)
				{
					if (eta != beta_syn)
					{
						int beta = GF.GFMultiply(h_inverse, GF.GFAdd(eta, beta_syn));
						L_c2v[row][dc][beta - 1] = -1.0 * TEMS_Lc2v[eta] - L0;
						L_c2v[row][dc][beta - 1] = L_c2v[row][dc][beta - 1] / TEMS_Factor;
						if (L_c2v[row][dc][beta - 1] < -1 * TEMS_Offset)
						{
							L_c2v[row][dc][beta - 1] = L_c2v[row][dc][beta - 1] + TEMS_Offset;
						}
						else if (L_c2v[row][dc][beta - 1] >  TEMS_Offset)
						{
							L_c2v[row][dc][beta - 1] = L_c2v[row][dc][beta - 1] - TEMS_Offset;
						}
						else
						{
							L_c2v[row][dc][beta - 1] = 0;
						}
					}
				}
			}
		}
	}
	tems_route_num = tems_route_num / iter_number;
	
	//OSD部分,
	if (OSD_Order >= 0) {
		if (decode_correct == false) {
			Decoding_OSD_bit(L_ch, DecodeOutput, RelySeri_bit, OSD_flag);
		}
	}
	
	delete[] old_L_v2c;
	return 0;
}

int CNBLDPC::Decoding_BS_TEMS(double** L_ch, int* DecodeOutput, int* RelySeri_symbol, int* RelySeri_bit)
{

	// initial
	tems_route_num = 0;
	for (int col = 0; col < CodeLen; col++)
	{
		for (int d = 0; d < VarDegree[col]; d++)
		{
			CopyLLRVector(L_v2c[col][d], L_ch[col]);
		}
	}
	for (int i = 0; i < maxIter; i++) {
		for (int j = 0; j < CodeLen; j++) {
			for (int k = 0; k < GFq - 1; k++) {
				LLR_Sym[i][j][k] = 0;
			}
		}
	}

	int Bit_Len = log(GFq * 1.0) / log(2.0);

	double* old_L_v2c = new double[GFq - 1];
	double* Prob_Symbl = new double[GFq];
	double* LLR_Bit = new double[Bit_Len];
	double* old_Prob_Symbl = new double[GFq];
	double* old_LLR_Bit = new double[Bit_Len];

	for (int i = 0; i < CodeLen; i++) {
		for (int j = 0; j < GFq - 1; j++) {
			LLR_Sym_sum[i][j] = 0;
		}
	}

	for (int row = 0; row < ChkLen; row++)
	{
		for (int d = 0; d < ChkDegree[row]; d++)
		{
			ClearLLRVector(L_c2v[row][d]);
		}
	}
	// iterative decoding
	int iter_number = 0;
	bool decode_correct = true;
	while (iter_number++ < maxIter)
	{
		//		cout << iter_number << " ";
		// tentative decoding
		for (int col = 0; col < CodeLen; col++)
		{
			CopyLLRVector(L_post[col], L_ch[col]);
			int row, dc;
			for (int d = 0; d < VarDegree[col]; d++)
			{
				row = VarLink[col][d];
				dc = VarLinkDc[col][d];
				AddLLRVector(L_post[col], L_post[col], L_c2v[row][dc]);
			}
			//OSD 每个符号的可信度计算
			for (int i = 0; i < GFq - 1; i++) {
				LLR_Sym_sum[col][i] = OSD_Alpha*LLR_Sym_sum[col][i] + L_post[col][i];
			}
			DecodeOutput[col] = DecideLLRVector(L_post[col]);
		}
		//保存每次迭代计算后的输出软量
		for (int d = 0; d < CodeLen; d++) {
			CopyLLRVector(LLR_Sym[iter_number - 1][d], L_post[d]);
		}

		// check whether satisfy the parity matrix
		decode_correct = true;
		for (int row = 0; row < ChkLen; row++)
		{
			int col, multip, syndrome;
			syndrome = 0;
			for (int d = 0; d < ChkDegree[row]; d++)
			{
				col = ChkLink[row][d];
				multip = GF.GFMultiply(ChkLinkGFe[row][d], DecodeOutput[col]);
				syndrome = GF.GFAdd(syndrome, multip);
			}
			if (syndrome != 0)
			{
				decode_correct = false;
				break;
			}
		}
		if (decode_correct)
		{
			delete[] old_L_v2c;
			delete[] Prob_Symbl;
			delete[] LLR_Bit;
			delete[] old_Prob_Symbl;
			delete[] old_LLR_Bit;
			OSD_permute(RelySeri_bit);
			tems_route_num = tems_route_num / iter_number;
			return 1;
		}
		// message from var to check

		for (int col = 0; col < CodeLen; col++)
		{
			int row, dc;
			int update_before = -1;
			int update_after = -1;
			for (int dv = 0; dv < VarDegree[col]; dv++)
			{
				row = VarLink[col][dv];
				dc = VarLinkDc[col][dv];
				update_before = DecideLLRVector(L_v2c[col][dv]);
				CopyLLRVector(old_L_v2c, L_v2c[col][dv]);
				MinusLLRVector(L_v2c[col][dv], L_post[col], L_c2v[row][dc]);
				update_after = DecideLLRVector(L_v2c[col][dv]);
				if (update_after != update_before)
				{
					for (int i = 0; i < GFq - 1; i++)
					{
						L_v2c[col][dv][i] = 0.25*old_L_v2c[i] + 0.75*L_v2c[col][dv][i];
					}


				}
			}
		}
		//cout << iter_number << " ";
		// message from check to var
		for (int row = 0; row < ChkLen; row++)
		{
			// get the most reliable symbol beta
			TEMS_Get_Beta(row);
			// transform to the delta domain
			TEMS_Get_deltaU(row);
			// get the two minimum value of each row and mark the Nr minimum pos
			TEMS_Get_Min(row);
			// mark the conf
			for (int q = 0; q < GFq; q++)
			{
				TEMS_deltaW[q] = DBL_MAX;
				TEMS_isSelected[q] = false;
			}
			////////////////////////////////////////////////////////////////////////////////
			for (int i = 0; i < maxChkDegree; i++) {
				col_isSelected[i] = false;
			}
			////////////////////////////////////////////////////////////////////////////////
			int sumNonele, diff;
			double sumNonLLR;
			sumNonele = diff = 0; sumNonLLR = 0.0;
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			//BS_TEMS与TEMS的核心区别
			TEMS_ConstructConf(1, 1, sumNonele, sumNonLLR, diff, 0, ChkDegree[row] - 1);
			Construct_BasicSet(basic_set_Nm);
			sumNonele = diff = 0; sumNonLLR = 0.0;
			BasicSet_TEMS_ConstructConf(0, basic_set_Nm - 1, diff, sumNonele, sumNonLLR);

			//////////////////////////////////////////////////////////////////////////////////////////////
			/*
			cout << "TEMS(1,1)" << endl;
			for (int i = 0; i < GFq; i++) {
			cout << i << "-" << TEMS_deltaW[i] << ":";
			for (int d = 0; d < ChkDegree[row]; d++) {
			cout << TEMS_Eta[i][d] << " ";
			}
			cout << endl;
			}
			cout << endl;
			*/
			/*
			cout << "basicSet" << endl;
			for (int i = 0; i < GFq; i++) {
			cout << i << "-" << TEMS_deltaW[i] << ":";
			for (int d = 0; d < ChkDegree[row]; d++) {
			cout << TEMS_Eta[i][d] << " ";
			}
			cout << endl;
			}
			cout << endl;
			*/
			/////////////////////////////////////////////////////////////////////////////////////////////////////////
			//test_sort();
			//test_ConstructConf(8, 3, sumNonele, sumNonLLR, diff, 0, ChkDegree[row] - 1);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//calculate the output
			for (int dc = 0; dc < ChkDegree[row]; dc++)
			{
				// clear the temp vector TEMS_Lc2v
				for (int q = 0; q < GFq; q++)
				{
					TEMS_Lc2v[q] = DBL_MAX;
					TEMS_isUpdate[q] = false;
				}
				// calculate the dc output
				for (int eta = 0; eta < GFq; eta++)
				{
					int deviation = TEMS_Eta[eta][dc];
					int eta_plus_deviation = GF.GFAdd(eta, deviation);
					double minus = TEMS_deltaW[eta] - TEMS_deltaU[dc][deviation];
					if (TEMS_Lc2v[eta_plus_deviation] > minus)
					{
						TEMS_Lc2v[eta_plus_deviation] = minus;
						TEMS_isUpdate[eta_plus_deviation] = true;
					}
				}
				for (int q = 0; q < GFq; q++)
				{
					if (!TEMS_isUpdate[q])
					{
						TEMS_Lc2v[q] = (dc == TEMS_Min[q][0]) ?
							TEMS_deltaU[TEMS_Min[q][1]][q] : TEMS_deltaU[TEMS_Min[q][0]][q];//补足的这一句很重要
					}
				}

				// delta domain --> normal domain & repermutation
				int h_inverse = GF.GFInverse(ChkLinkGFe[row][dc]);
				int beta_syn = GF.GFAdd(TEMS_Syndrome, TEMS_Beta[dc]);
				double L0 = -1.0 * TEMS_Lc2v[beta_syn];
				for (int eta = 0; eta < GFq; eta++)
				{
					if (eta != beta_syn)
					{
						int beta = GF.GFMultiply(h_inverse, GF.GFAdd(eta, beta_syn));
						L_c2v[row][dc][beta - 1] = -1.0 * TEMS_Lc2v[eta] - L0;
						L_c2v[row][dc][beta - 1] = L_c2v[row][dc][beta - 1] / basic_set_factor;
						if (L_c2v[row][dc][beta - 1] < -1 * basic_set_offset)
						{
							L_c2v[row][dc][beta - 1] = L_c2v[row][dc][beta - 1] + basic_set_offset;
						}
						else if (L_c2v[row][dc][beta - 1] >  basic_set_offset)
						{
							L_c2v[row][dc][beta - 1] = L_c2v[row][dc][beta - 1] - basic_set_offset;
						}
						else
						{
							L_c2v[row][dc][beta - 1] = 0;
						}
					}
				}
			}
		}
	}
	tems_route_num = tems_route_num / iter_number;

	//OSD部分,
	if (OSD_Order >= 0) {
		if (decode_correct == false) {
			Decoding_OSD_bit(L_ch, DecodeOutput, RelySeri_bit, OSD_flag);
		}
	}

	delete[] old_L_v2c;
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//通用辅助函数
int CNBLDPC::CopyLLRVector(double* L_d, double* L_s)
{
	for (int q = 0; q < GFq - 1; q++)
	{
		L_d[q] = L_s[q];
	}
	return 0;
}

int CNBLDPC::ClearLLRVector(double* L)
{
	for (int q = 0; q < GFq - 1; q++)
	{
		L[q] = 0;
	}
	return 0;
}

int CNBLDPC::AddLLRVector(double* L, double* L1, double* L2)
{
	for (int q = 0; q < GFq - 1; q++)
	{
		L[q] = L1[q] + L2[q];
	}
	return 0;
}

int CNBLDPC::MinusLLRVector(double* L, double* L1, double* L2)
{
	for (int q = 0; q < GFq - 1; q++)
	{
		L[q] = L1[q] - L2[q];
	}
	return 0;
}

int CNBLDPC::ExchRow(int** H, int source, int destination)
{
	int temp;
	for (int col = 0; col < CodeLen; col++)
	{
		temp = H[source][col];
		H[source][col] = H[destination][col];
		H[destination][col] = temp;
		//		H[source][col] =  H[source][col] + H[destination][col];
		//		H[destination][col] = H[source][col] - H[destination][col];
		//		H[source][col] = H[source][col] - H[destination][col];
	}
	return 0;
}

int CNBLDPC::AddRow(int** H, int source, int destination, int multifactor)
{
	for (int col = 0; col < CodeLen; col++)
	{
		H[destination][col] = GF.GFAdd(H[destination][col], GF.GFMultiply(multifactor, H[source][col]));
	}
	return 0;
}

int CNBLDPC::ExchCol(int** H, int source, int destination)
{
	int temp;
	for (int row = 0; row < ChkLen; row++)
	{
		temp = H[row][source];
		H[row][source] = H[row][destination];
		H[row][destination] = temp;
	}
	return 0;
}


//编码辅助函数
int CNBLDPC::GaussEliminate(int** H, int* ExchangeColSource, int* ExchangeColDest, int& ExchangeColNum)
{
	//用于encode的高斯消去
	for (int row = ChkLen - 1; row >= 0; row--)
	{
		// eliminate row 'row'
		int col = row + CodeLen - ChkLen;
		// make the H[row][col] != 0
		if (H[row][col] == 0)
		{
			bool exchanged = false;
			// search up
			if (!exchanged)
			{
				for (int row_up = row - 1; row_up >= 0; row_up--)
				{
					if (H[row_up][col] != 0)
					{
						ExchRow(H, row, row_up);
						exchanged = true;
						break;
					}
				}
			}
			// if cannot find nonzero ele, search left
			if (!exchanged)
			{
				for (int col_left = col - 1; col_left >= 0; col_left--)
				{
					if (H[row][col_left] != 0)
					{
						ExchCol(H, col, col_left);
						EncodeExchangeColSource[EncodeExchangeColNum] = col;
						EncodeExchangeColDest[EncodeExchangeColNum] = col_left;
						EncodeExchangeColNum++;
						exchanged = true;
						break;
					}
				}
			}
			if (!exchanged)
			{
				cerr << "NB matrix is not full rank" << endl;
				exit(-1);
			}
		}
		// make the H[0:row-1][col] = 0
		int h_inverse = GF.GFInverse(H[row][col]);
		for (int row_up = row - 1; row_up >= 0; row_up--)
		{
			int x;
			if (H[row_up][col] != 0)
			{
				x = GF.GFMultiply(h_inverse, H[row_up][col]);
				AddRow(H, row, row_up, x);
			}
		}
		// make the H[row][col] = 1
		for (int col_left = 0; col_left <= col; col_left++)
		{
			H[row][col_left] = GF.GFMultiply(H[row][col_left], h_inverse);
		}
	}
	return 0;
}


//译码辅助函数
int CNBLDPC::DecideLLRVector(double* LLR)
{
	double max = 0;
	int alpha_i;
	for (int q = 0; q < GFq - 1; q++)
	{
		if (LLR[q] > max)
		{
			max = LLR[q];
			alpha_i = q + 1;
		}
	}
	if (max <= 0)
	{
		return 0;
	}
	else
	{
		return alpha_i;
	}
}

//BP
int CNBLDPC::L_Back(int row, int l)
{
	int A1, A2;
	int col, dv;
	if (l < 0)
	{
		return 0;
	}
	else if (l == 0)
	{
		A1 = 0;
		A2 = ChkLinkGFe[row][l];
		col = ChkLink[row][l];
		dv = ChkLinkDv[row][l];
		LLR_BoxPlus(L_sigma, L_sigma, L_v2c[col][dv], A1, A2);
	}
	else
	{
		L_Back(row, l - 1);
		A1 = 1;
		A2 = ChkLinkGFe[row][l];
		col = ChkLink[row][l];
		dv = ChkLinkDv[row][l];
		LLR_BoxPlus(L_sigma, L_sigma, L_v2c[col][dv], A1, A2);
	}
	return 1;
}

int CNBLDPC::L_Forward(int row, int l)
{
	int A1, A2;
	int col, dv;
	if (l >= ChkDegree[row])
	{
		return 0;
	}
	else if (l == ChkDegree[row] - 1)
	{
		A1 = 0;
		A2 = ChkLinkGFe[row][l];
		col = ChkLink[row][l];
		dv = ChkLinkDv[row][l];
		LLR_BoxPlus(L_rho, L_rho, L_v2c[col][dv], A1, A2);
	}
	else
	{
		L_Forward(row, l + 1);
		A1 = 1;
		A2 = ChkLinkGFe[row][l];
		col = ChkLink[row][l];
		dv = ChkLinkDv[row][l];
		LLR_BoxPlus(L_rho, L_rho, L_v2c[col][dv], A1, A2);
	}
	return 1;
}

int CNBLDPC::LLR_BoxPlus(double* L, double* L1, double* L2, int A1, int A2)
{
	if (A1 == 0)
	{
		int v2;
		for (int alpha_i = 1; alpha_i < GFq; alpha_i++)
		{
			v2 = GF.GFMultiply(GF.GFInverse(A2), alpha_i);
			LLRV[alpha_i - 1] = L2[v2 - 1];
		}
		CopyLLRVector(L, LLRV);
	}
	else if (A2 == 0)
	{
		int v1;
		for (int alpha_i = 1; alpha_i < GFq; alpha_i++)
		{
			v1 = GF.GFMultiply(GF.GFInverse(A1), alpha_i);
			LLRV[alpha_i - 1] = L1[v1 - 1];
		}
		CopyLLRVector(L, LLRV);
	}
	else
	{
		long double sum1, sum2;
		//		long double sum1_prime, sum2_prime;
		int v1, v2;
		// down
		sum2 = 0;
		//		sum2_prime = 0;
		for (int x = 0; x < GFq; x++)
		{
			v1 = x;
			v2 = GF.GFMultiply(GF.GFInverse(A2), GF.GFMultiply(A1, x));
			//			if(v1 == 0 || v2 == 0)
			//			{
			//				sum2_prime = 0;
			//				sum2 += 1;
			//			}
			//			else
			if (v1 != 0 && v2 != 0)
			{
				//				sum2_prime = log(sum2);
				if (sum2 > L1[v1 - 1] + L2[v2 - 1])
				{
					sum2 = sum2 + log(1 + exp(-1 * (sum2 - (L1[v1 - 1] + L2[v2 - 1]))));
				}
				else
				{
					sum2 = L1[v1 - 1] + L2[v2 - 1] + log(1 + exp(-1 * (L1[v1 - 1] + L2[v2 - 1] - sum2)));
				}
			}
		}
		for (int alpha_i = 1; alpha_i < GFq; alpha_i++)
		{
			// up
			sum1 = 0;
			//			sum1_prime = 0;
			v1 = GF.GFMultiply(GF.GFInverse(A1), alpha_i);
			v2 = GF.GFMultiply(GF.GFInverse(A2), alpha_i);
			if (L1[v1 - 1] > L2[v2 - 1])
			{
				sum1 = L1[v1 - 1] + log(1 + exp(-1 * (L1[v1 - 1] - L2[v2 - 1])));
			}
			else
			{
				sum1 = L2[v2 - 1] + log(1 + exp(-1 * (L2[v2 - 1] - L1[v1 - 1])));
			}
			for (int x = 0; x < GFq; x++)
			{
				// up
				v1 = x;
				v2 = GF.GFMultiply(GF.GFInverse(A2), GF.GFAdd(alpha_i, GF.GFMultiply(x, A1)));
				if (v1 != 0 && v2 != 0)
				{
					//					sum1_prime = log(sum1);
					if (sum1 > L1[v1 - 1] + L2[v2 - 1])
					{
						sum1 = sum1 + log(1 + exp(-1 * (sum1 - (L1[v1 - 1] + L2[v2 - 1]))));
					}
					else
					{
						sum1 = L1[v1 - 1] + L2[v2 - 1] + log(1 + exp(-1 * (L1[v1 - 1] + L2[v2 - 1] - sum1)));
					}
				}
			}
			LLRV[alpha_i - 1] = sum1 - sum2;
		}
		CopyLLRVector(L, LLRV);
	}
	return 0;
}

//EMS
int CNBLDPC::SortLLRVector(double* L_sorted, int* Enre_sorted, int Len_sorted, double* L_sorting, int Len_sorting)
{
	int sortedCount = 0;
	L_sorted[0] = 0;
	Enre_sorted[0] = 0;
	double cand;
	while (sortedCount < Len_sorting)
	{
		// the candidate value
		cand = L_sorting[sortedCount];
		// put the cand on the tail
		L_sorted[sortedCount + 1] = cand;
		Enre_sorted[sortedCount + 1] = sortedCount + 1;
		// search forward
		for (int i = sortedCount; i >= 0; i--)
		{
			if (cand >= L_sorted[i])
			{
				L_sorted[i + 1] = L_sorted[i];
				L_sorted[i] = cand;
				Enre_sorted[i + 1] = Enre_sorted[i];
				Enre_sorted[i] = sortedCount + 1;
			}
			else
			{
				break;
			}
		}
		sortedCount++;
	}
	return 0;
}

int CNBLDPC::ConstructConf(int Nm, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int except, int end, int row)
{
	if (begin > end)
	{
		if (sumNonLLR > EMS_L_c2v[sumNonele])
		{
			EMS_L_c2v[sumNonele] = sumNonLLR;
		}
	}
	else if (begin == except)
	{
		ConstructConf(Nm, Nc, sumNonele, sumNonLLR, diff, begin + 1, except, end, row);
		return 0;
	}
	else
	{
		for (int k = 0; k < Nm; k++)
		{
			sumNonele = GF.GFAdd(GF.GFMultiply(EMS_sort_Entr_v2c[begin][k], ChkLinkGFe[row][begin]), sumNonele);
			sumNonLLR = sumNonLLR + EMS_sort_L_v2c[begin][k];
			diff += (k != 0) ? 1 : 0;
			if (diff <= Nc)
			{
				ConstructConf(Nm, Nc, sumNonele, sumNonLLR, diff, begin + 1, except, end, row);
				sumNonele = GF.GFAdd(GF.GFMultiply(EMS_sort_Entr_v2c[begin][k], ChkLinkGFe[row][begin]), sumNonele);
				sumNonLLR = sumNonLLR - EMS_sort_L_v2c[begin][k];
				diff -= (k != 0) ? 1 : 0;
			}
			else
			{
				sumNonele = GF.GFAdd(GF.GFMultiply(EMS_sort_Entr_v2c[begin][k], ChkLinkGFe[row][begin]), sumNonele);
				sumNonLLR = sumNonLLR - EMS_sort_L_v2c[begin][k];
				diff -= (k != 0) ? 1 : 0;
				break;
			}
		}
	}
	return 0;
}

//TEMS
int CNBLDPC::TEMS_Get_Beta(int row)
{

	TEMS_Syndrome = 0;
	for (int dc = 0; dc < ChkDegree[row]; dc++)
	{
		int col = ChkLink[row][dc];
		int dv = ChkLinkDv[row][dc];
		int h = ChkLinkGFe[row][dc];
		double max = 0;
		int max_ele = 0;
		for (int q = 1; q < GFq; q++)
		{
			if (L_v2c[col][dv][q - 1] > max)
			{
				max = L_v2c[col][dv][q - 1];
				max_ele = GF.GFMultiply(q, h);
			}
		}
		TEMS_Beta[dc] = max_ele;
		TEMS_Syndrome = GF.GFAdd(TEMS_Syndrome, max_ele);
	}
	return 0;
}

int CNBLDPC::TEMS_Get_deltaU(int row)
{
	for (int dc = 0; dc < ChkDegree[row]; dc++)
	{
		int col = ChkLink[row][dc];
		int dv = ChkLinkDv[row][dc];
		int h_inverse = GF.GFInverse(ChkLinkGFe[row][dc]);

		int beta_p = GF.GFMultiply(h_inverse, TEMS_Beta[dc]);
		double max;
		max = (beta_p == 0) ? 0 : L_v2c[col][dv][beta_p - 1];////////////////////////////////////////////////////////////////////

		TEMS_deltaU[dc][TEMS_Beta[dc]] = max - 0;
		for (int x = 1; x < GFq; x++)
		{
			int eta = GF.GFAdd(x, TEMS_Beta[dc]);
			TEMS_deltaU[dc][eta] = max - L_v2c[col][dv][GF.GFMultiply(h_inverse, x) - 1];
		}
	}
	return 0;
}

int CNBLDPC::TEMS_Get_Min(int row)
{
	// sort
	for (int q = 0; q < GFq; q++)
	{
		// clear
		for (int dc = 0; dc < ChkDegree[row]; dc++)
		{
			TEMS_Min[q][dc] = dc;
		}
		// search min
		double min = TEMS_deltaU[0][q];
		int min_pos = 0;
		for (int dc = 1; dc < ChkDegree[row]; dc++)
		{
			for (int d_s = dc; d_s >= 1; d_s--)
			{
				double du = TEMS_deltaU[TEMS_Min[q][d_s]][q];
				double du_forward = TEMS_deltaU[TEMS_Min[q][d_s - 1]][q];
				if (du < du_forward)
				{
					int temp;
					temp = TEMS_Min[q][d_s];
					TEMS_Min[q][d_s] = TEMS_Min[q][d_s - 1];
					TEMS_Min[q][d_s - 1] = temp;
				}
				else
				{
					break;
				}
			}
		}
	}
	// set TEMS_inConf_Nr
	for (int dc = 0; dc < ChkDegree[row]; dc++)
	{
		TEMS_inConf_Nr[0][dc] = true;
	}
	for (int q = 1; q < GFq; q++)
	{
		for (int dc = 0; dc < ChkDegree[row]; dc++)
		{
			TEMS_inConf_Nr[q][dc] = false;
		}
		for (int index = 0; index < TEMS_Nr; index++)
		{
			if (index >= ChkDegree[row])
				break;
			int dc = TEMS_Min[q][index];
			TEMS_inConf_Nr[q][dc] = true;
		}
	}

	return 0;
}

int CNBLDPC::TEMS_ConstructConf(int Nr, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int end)
{
	if (begin > end)
	{
		tems_route_num++;
		if (sumNonLLR < TEMS_deltaW[sumNonele])
		{
			TEMS_deltaW[sumNonele] = sumNonLLR;
			for (int d = 0; d <= end; d++)
			{
				TEMS_Eta[sumNonele][d] = TEMS_Eta_Candidate[d];
			}
		}
	}
	else
	{
		for (int q = 0; q < GFq; q++)
		{
			if (TEMS_inConf_Nr[q][begin] && !TEMS_isSelected[q])
			{
				diff += (q != 0) ? 1 : 0;
				if (diff <= Nc)
				{
					TEMS_isSelected[q] = (q != 0) ? true : false;//保证一条路径中不存在两个校验节点q是相同的。
					sumNonele = GF.GFAdd(sumNonele, q);
					sumNonLLR += TEMS_deltaU[begin][q];
					TEMS_Eta_Candidate[begin] = q;
					TEMS_ConstructConf(Nr, Nc, sumNonele, sumNonLLR, diff, begin + 1, end);//以下语句是迭代回退
					sumNonele = GF.GFAdd(sumNonele, q);
					sumNonLLR -= TEMS_deltaU[begin][q];
					TEMS_isSelected[q] = false;//
					diff -= (q != 0) ? 1 : 0;
				}
				else
				{
					diff -= (q != 0) ? 1 : 0;
					break;
				}

				/*
				else
				{
				sumNonele = GF.GFAdd(sumNonele, q);
				sumNonLLR -= TEMS_deltaU[begin][q];
				diff -= (q != 0)? 1:0;
				break;
				}
				*/
			}
		}
	}
	return 0;
}
//*****************************test_conf***************************************
void CNBLDPC::test_sort() {
	for (int i = 0; i < maxChkDegree; i++) {
		vector<double> temp_sort;
		for (int j = 0; j < GFq; j++) {
			temp_sort.push_back(TEMS_deltaU[i][j]);
			TEMS_sort_deltaU[i][j] = TEMS_deltaU[i][j];
			TEMS_sort_q[i][j] = j;
		}
		sort(TEMS_sort_q[i], TEMS_sort_q[i] + GFq, [&temp_sort](int t1, int t2) {return temp_sort[t1] < temp_sort[t2]; });
		sort(TEMS_sort_deltaU[i], TEMS_sort_deltaU[i] + GFq);
	}
	/*
	double temp_L_forward;
	double temp_L;
	int temp_q;
	for (int i = 0; i < maxChkDegree; i++) {
		TEMS_sort_deltaU[i][0] = TEMS_deltaU[i][0];
		TEMS_sort_q[i][0] = 0;
		for (int j = 1; j < GFq; j++) {
			TEMS_sort_deltaU[i][j] = TEMS_deltaU[i][j];
			TEMS_sort_q[i][j] = j;
			for (int k = j; k >= 1; k--) {
				temp_L = TEMS_sort_deltaU[i][k];
				temp_L_forward = TEMS_sort_deltaU[i][k - 1];
				if (temp_L < temp_L_forward) {
					temp_q = TEMS_sort_q[i][k];
					TEMS_sort_q[i][k] = TEMS_sort_q[i][k - 1];
					TEMS_sort_q[i][k - 1] = temp_q;

					TEMS_sort_deltaU[i][k] = temp_L_forward;
					TEMS_sort_deltaU[i][k - 1] = temp_L;
				}
			}
		}
	}
	*/
}

int CNBLDPC::test_ConstructConf(int Nm, int Nc, int& sumNonele, double& sumNonLLR, int& diff, int begin, int end)
{
	if (begin > end)
	{
		if (sumNonLLR < TEMS_deltaW[sumNonele])
		{
			TEMS_deltaW[sumNonele] = sumNonLLR;
			for (int d = 0; d <= end; d++)
			{
				TEMS_Eta[sumNonele][d] = TEMS_Eta_Candidate[d];
			}
		}
	}
	else
	{
		for (int k = 0; k < Nm; k++)
		{
			diff += (k != 0) ? 1 : 0;
			if (diff <= Nc)
			{
				sumNonele = GF.GFAdd(sumNonele, TEMS_sort_q[begin][k]);
				sumNonLLR = sumNonLLR + TEMS_sort_deltaU[begin][k];
				TEMS_Eta_Candidate[begin] = TEMS_sort_q[begin][k];
				test_ConstructConf(Nm, Nc, sumNonele, sumNonLLR, diff, begin + 1, end);
				sumNonele = GF.GFAdd(sumNonele, TEMS_sort_q[begin][k]);
				sumNonLLR = sumNonLLR - TEMS_sort_deltaU[begin][k];
				diff -= (k != 0) ? 1 : 0;
			}
			else
			{
				diff -= (k != 0) ? 1 : 0;
				break;
			}
		}
	}
	return 0;
}

//**********************************************Basic_Set*******************************************************************
void CNBLDPC::Construct_BasicSet(int Nm) {

	vector<int> order;
	vector<double> LLV;
	vector<int> R;
	int R_num = 0;
	int temp_q;
	int temp_sum_q;
	//initial
	for (int i = 0; i < GFq; i++) {
		Set_q[i] = i;
		Set_col[i] = TEMS_Min[i][0];
		LLV.push_back(TEMS_deltaU[TEMS_Min[i][0]][i]);
		order.push_back(i);
		basicSet_available[i] = true;
	}
	basicSet_available[0] = false;
	sort(order.begin(), order.end(), [&LLV](int t1, int t2) {return LLV[t1] < LLV[t2]; });
	
	int bit_len = log(GFq) / log(2);
	if (Nm > bit_len) {
		for (int i = 0; i < (GFq-1); i++) {
			basicSet_q[i] = Set_q[order[i+1]];
			basicSet_col[i] = Set_col[order[i+1]];
			basicSet_L[i] = LLV[order[i+1]];
		}
		///////////////////////////////////////////////////////////////
	   /*
		for (int i = 0; i < GFq; i++) {
			cout << "q:" << basicSet_q[i] << " col:" << basicSet_col[i] << " L:" << basicSet_L[i] << endl;
		}
		*/
		///////////////////////////////////////////////////////////////
	}
	else {
		//set the first elemnet of BasicSet,in T-EMS order[0] is always q = 0 
		basicSet_q[0] = Set_q[order[1]];
		basicSet_col[0] = Set_col[order[1]];
		basicSet_L[0] = LLV[order[1]];


		R.push_back(basicSet_q[0]);
		basicSet_available[R[0]] = false;
		R_num += 1;
		int temp_num = 1;
		//find the rest elements
		for (int i = 2; i < GFq; i++) {
			if (R_num >= (GFq - 1)) { break; }//if all the elements have been found
			else {
				temp_q = Set_q[order[i]];//if this q is independent
				if (basicSet_available[temp_q]) {
					basicSet_q[temp_num] = Set_q[temp_q];
					basicSet_col[temp_num] = Set_col[temp_q];
					basicSet_L[temp_num] = LLV[temp_q];
					temp_num++;

					for (int j = 0; j <R_num; j++) {
						temp_sum_q = GF.GFAdd(temp_q, R[j]);
						R.push_back(temp_sum_q);
						basicSet_available[temp_sum_q] = false;
					}
					R.push_back(temp_q);
					basicSet_available[temp_q] = false;
					R_num = 2 * R_num + 1;
				}
			}
		}

		////////////////////////////////////////
		/*
		cout << "迭代" << endl;
		for (int i = 0; i < temp_num; i++) {
		cout << basicSet_q[i] << " " << basicSet_col[i] << " " << basicSet_L[i] << "; ";
		}
		cout << endl;
		*/
		////////////////////////////////////////////////////
	}
}
int CNBLDPC::BasicSet_TEMS_ConstructConf(int begin, int end, int &diff,int& sumNonele, double& sumNonLLR)
{
	if (begin > end)
	{
		if (sumNonLLR < TEMS_deltaW[sumNonele])
		{
			TEMS_deltaW[sumNonele] = sumNonLLR;
			//////////////////////////////////////////////////////////
			for (int d = 0; d < maxChkDegree; d++)////////////////////////////////////////////////
			{
				TEMS_Eta[sumNonele][d] = TEMS_Eta_Candidate[d];
			}
			////////////////////////////////////////////////////////////
		}
	}
	else {
		int temp_col = basicSet_col[begin];	
		if (col_isSelected[temp_col] == false) {
			diff += 1;
			if(diff<= basic_set_Nc){
				col_isSelected[temp_col] = true;
				sumNonele = GF.GFAdd(sumNonele, basicSet_q[begin]);
				sumNonLLR += basicSet_L[begin];
				//////////////////////////////////////////////////////////
				TEMS_Eta_Candidate[temp_col] = basicSet_q[begin];
				////////////////////////////////////////////////////////////
				BasicSet_TEMS_ConstructConf(begin + 1, end, diff, sumNonele, sumNonLLR);//以下语句是迭代回退
				sumNonele = GF.GFAdd(sumNonele, basicSet_q[begin]);
				sumNonLLR -= basicSet_L[begin];
				col_isSelected[temp_col] = false;
				diff -= 1;
				TEMS_Eta_Candidate[temp_col] = 0;
				BasicSet_TEMS_ConstructConf(begin + 1, end, diff, sumNonele, sumNonLLR);
			}
			else {
				diff -= 1;
				TEMS_Eta_Candidate[temp_col] = 0;
				BasicSet_TEMS_ConstructConf(begin + 1, end, diff, sumNonele, sumNonLLR);
			}
		}
		else {
			BasicSet_TEMS_ConstructConf(begin + 1, end, diff, sumNonele, sumNonLLR);
		}
		
	}
	return 0;
}