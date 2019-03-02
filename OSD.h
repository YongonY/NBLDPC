#pragma once
#include "NBLDPC.h"
#include "CRC.h"

int CNBLDPC::Decoding_OSD_bit(double** L_ch, int* DecodeOutput, int*ReliableSeri_BIT,int flag) {
	//step1: initial
	if (flag) {//flag = 1 表示前面没有级联译码算法，需要自己定义log_Sym_sum
		for (int i = 0; i < CodeLen; i++) {
			CopyLLRVector(LLR_Sym_sum[i], L_ch[i]);
			infor_record[i] = DecideLLRVector(L_ch[i]);
			near_optimal[i] = infor_record[i];//DecodeOutput是TEMS输出的译码结果
			DecodeOutput[i] = infor_record[i];
		}
	}
	else {
		for (int i = 0; i < CodeLen; i++) {
			near_optimal[i] = DecodeOutput[i];//DecodeOutput是上一级输出的译码结果
			infor_record[i] = DecodeOutput[i];
		}
	}
	min_distance = 1000000;
	//symbol2bit
	int bitlen = log(GFq) / log(2);
	int temp3;

	for (int i = 0; i < CodeLen; i++) {
		temp3 = GFq - 1;
		for (int j = 0; j < bitlen; j++) {
			//infor_record_bit[i * bitlen + j] = TX_code[i * bitlen + j];
			infor_record_bit[i * bitlen + j] = ((infor_record[i] & (1 << j)) == 0) ? 0 : 1;
			near_optimal_bit[i * bitlen + j] = ((near_optimal[i] & (1 << j)) == 0) ? 0 : 1;
			temp3 = temp3 / 2;
			L_ch_bit[i*bitlen+bitlen-1-j] = L_ch[i][temp3];
		}
	}
	
	
//step2 permute
	OSD_permute(ReliableSeri_BIT);
	
//step3: Gauss eliminate
	///////////////////////////////////////////////////////////////////////////
	//添加CRC功能
	int MsgLen_bit = CodeLen_bit - ChkLen_bit;
	for (int i = 0; i < CRCrate; i++) {
		for (int j = 0; j < MsgLen_bit; j++) {
			H_OSD_CRC_bit[i][j] = CRCMatrix_partH_PI[i][j];
		}
		for (int j = MsgLen_bit; j < CodeLen_bit; j++) {
			H_OSD_CRC_bit[i][j] = 0;
		}
	}
	for (int i = CRCrate; i < CRCrate + ChkLen_bit; i++) {
		for (int j = 0; j < CodeLen_bit; j++) {
			H_OSD_CRC_bit[i][j] = H_origin_bit[i - CRCrate][j];
		}
	}

	H_GaussEliminate_bit(H_OSD_CRC_bit, H_OSD_CRC_bit, CodeLen_bit, CRCrate + ChkLen_bit, ReliableSeri_BIT);
	//////////////////////////////////////////////////////////////////////////
	
	/*
	cout << endl;
	cout << "check" << endl;
	for (int i = 0; i < ChkLen_bit+CRCrate; i++) {
		int temp_check = 0;
		for (int j = 0; j < CodeLen_bit; j++) {
			//temp_check = H_origin_bit[i][j];
			temp_check = temp_check ^ (TX_code[j] * H_OSD_CRC_bit[i][j]);
		}
		cout << temp_check;
	}
	cout << endl;
	
	cout << "CRCMAtrixG" << endl;
	int test_CRCG[1000] = { 0 };
	for (int i = 0; i < MsgLen_bit - CRCLen; i++) {
		test_CRCG[i] = TX_code[i];
	}
	for (int i = 0; i < CRCLen; i++) {
		for (int j = 0; j < MsgLen_bit - CRCLen; j++) {
			test_CRCG[i + MsgLen_bit - CRCLen] = test_CRCG[i + MsgLen_bit - CRCLen] ^ (test_CRCG[j] * CRCMatrix_G_IP[j][i + MsgLen_bit - CRCLen]);
		}
	}
	for (int i = 0; i < 64; i++) {
		cout << test_CRCG[i];
	}
	cout << endl;
	
	cout << "反编码之前" << endl;
	for (int j = 0; j < CodeLen_bit; j++) {
		cout << TX_code[j];
	}
	cout << endl;

	cout << "check" << endl;
	for (int i = 0; i < CRCrate; i++) {
		int temp_check = 0;
		for (int j = 0; j < CodeLen_bit; j++) {
			temp_check = temp_check ^ ((TX_code[j])* H_OSD_CRC_bit[i][j]);
		}
		cout << temp_check;
	}
	cout << endl;
	
	*/
	

	/////////////////////////////////////////////////////////
	/*
	cout << "check" << endl;
	for (int i = 0; i < CRCrate + ChkLen_bit; i++) {
		for (int j = 0; j < CodeLen_bit - CRCrate - ChkLen_bit; j++) {
			cout << H_OSD_CRC_bit[i][ReliableSeri_BIT[j]];
		}
		cout << endl;
	}
	cout << "check2" << endl;
	for (int i = 0; i < CRCrate + ChkLen_bit; i++) {
	for (int j = CodeLen_bit - CRCrate - ChkLen_bit; j < CodeLen_bit; j++) {
	cout << H_OSD_CRC_bit[i][ReliableSeri_BIT[j]];
	}
	cout << endl;
	}
	*/

	//////////////////////////////////////////////////////////


	//step4:翻转比特
	//0阶OSD

	OSD_Encode_bit(H_OSD_CRC_bit, infor_record_bit,OSD_temp_out_bit,ReliableSeri_BIT,CodeLen_bit,CRCrate + ChkLen_bit);
	compute_min_distance_bit();

	
	
	//1阶OSD
	if (OSD_Order >= 1) {
		//cout << "uu";
		for (int i = 0; i < MsgLen_bit; i++) {
			infor_record_bit[ReliableSeri_BIT[i]] ^= 1;
			OSD_Encode_bit(H_OSD_CRC_bit, infor_record_bit, OSD_temp_out_bit, ReliableSeri_BIT, CodeLen_bit, CRCrate + ChkLen_bit);
			compute_min_distance_bit();
			infor_record_bit[ReliableSeri_BIT[i]] ^= 1;
		}
	}

	//2阶OSD
	if (OSD_Order >= 2) {
		//cout << "ee";
		for (int i = 0; i < MsgLen_bit; i++) {
			infor_record_bit[ReliableSeri_BIT[i]] ^= 1;
			for (int j = i + 1; j < MsgLen_bit; j++) {
				infor_record_bit[ReliableSeri_BIT[j]] ^= 1;
				OSD_Encode_bit(H_OSD_CRC_bit, infor_record_bit, OSD_temp_out_bit, ReliableSeri_BIT, CodeLen_bit, CRCrate + ChkLen_bit);
				compute_min_distance_bit();
				infor_record_bit[ReliableSeri_BIT[j]] ^= 1;
			}
			infor_record_bit[ReliableSeri_BIT[i]] ^= 1;
		}
	}

	//3阶OSD
	if (OSD_Order >= 3) {
		cout << "mm";
		for (int i = 0; i < MsgLen_bit; i++) {
			infor_record_bit[ReliableSeri_BIT[i]] ^= 1;
			for (int j = i + 1; j < MsgLen_bit; j++) {
				infor_record_bit[ReliableSeri_BIT[j]] ^= 1;
				for (int k = j + 1; k < MsgLen_bit; k++) {
					infor_record_bit[ReliableSeri_BIT[k]] ^= 1;
					OSD_Encode_bit(H_OSD_CRC_bit, infor_record_bit, OSD_temp_out_bit, ReliableSeri_BIT, CodeLen_bit, CRCrate + ChkLen_bit);
					compute_min_distance_bit();
					infor_record_bit[ReliableSeri_BIT[k]] ^= 1;
				}
				infor_record_bit[ReliableSeri_BIT[j]] ^= 1;
			}
			infor_record_bit[ReliableSeri_BIT[i]] ^= 1;
		}
	}
	

	for (int i = 0; i < CodeLen; i++) {
		DecodeOutput[i] = 0;
		for (int j = 0; j < bitlen; j++) {
			DecodeOutput[i] = 2*DecodeOutput[i]+ near_optimal_bit[i*bitlen +bitlen-1-j];
		}
	}
	/*
	cout << "OSD输出结果" << endl;
	for (int j = 0; j < CodeLen_bit; j++) {
		cout << near_optimal_bit[j];
	}
	cout << endl;
	*/
	//int ret = paritycheck(DecodeOutput);
	return 0;
}

//依据累积似然比获得符号（比特）的可信度排序
void CNBLDPC::OSD_permute(int* ReliableSeri_BIT) {
	//注意：依赖外部变量：LLR_Sym_sum

	//initial
	vector<int> col_bit;
	vector<double> Belief_bit;

	int Bit_Len = log(GFq * 1.0) / log(2.0);
	double *temp;
	temp = new double[Bit_Len];

	//bit permute
	for (int j = 0; j < CodeLen; j++) {
		for (int k = 0; k < Bit_Len; k++) {
			temp[k] = LLR_Sym_sum[j][(1 << k) -1];
		}
		for (int i = 0; i < Bit_Len; i++) {
			col_bit.push_back(j*Bit_Len + i);
			Belief_bit.push_back(abs(temp[i]));
		}
	}
	sort(col_bit.begin(), col_bit.end(), [&Belief_bit](int t1, int t2) {return Belief_bit[t1] > Belief_bit[t2]; });

	//end
	for (int i = 0; i < CodeLen*Bit_Len; i++) {
		ReliableSeri_BIT[i] = col_bit[i];//将一帧中的比特依据软量从大到小排序
	}

	delete[] temp;
}

int CNBLDPC::paritycheck(int*DecodeOutput) {
	int decode_correct = 1;
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
			decode_correct = 0;
			break;
		}
	}
	if (decode_correct)
	{
		return decode_correct;
	}
}


void CNBLDPC::G_GaussEliminate_bit(int**Matrix_in, int** MAtrix_out,int col_length,int row_length, int*order_seri) {
	//当对生成矩阵进行高斯消去。
	//step1: Initial
	for (int i = 0; i < row_length; i++) {
		for (int j = 0; j < col_length; j++) {
			MAtrix_out[i][j] = Matrix_in[i][j];
		}
	}

	//step2:gauss eliminate
	for (int row = 0; row <row_length; row++)
	{
		int col = order_seri[row];
		if (MAtrix_out[row][col] == 0)
		{
			bool exchanged = false;
			// search up
			if (!exchanged)
			{
				for (int row_up = row + 1 ; row_up<row_length; row_up++)
				{
					if (MAtrix_out[row_up][col] != 0)
					{
						for (int i = 0; i < col_length; i++) {
							MAtrix_out[row][i] = MAtrix_out[row][i]^MAtrix_out[row_up][i];
						}
						exchanged = true;
						break;
					}
				}
			}
			// if cannot find nonzero ele, search left
			if (!exchanged)
			{
				int flag = row;
				row--;
				for (int i = flag; i<col_length-1; i++)
				{
					int col_left = order_seri[i];
					order_seri[i] = order_seri[i + 1];
					order_seri[i + 1] = col_left;

				}
			}
		}
		// make the Matrix_out[row+1:row_length][col] = 0
		for (int row_up = row + 1; row_up<row_length; row_up++)
		{
			if (MAtrix_out[row_up][col] != 0)
			{
				for (int i = 0; i < col_length; i++) {
					MAtrix_out[row_up][i] = MAtrix_out[row_up][i]^MAtrix_out[row][i];
				}
			}
		}
	}
	// make the MAtrix_out彻底变成一个标准型
	for (int row = row_length - 1; row > 0; row--) {
		for (int row_up = row - 1; row_up >= 0; row_up--) {
			if (MAtrix_out[row_up][order_seri[row]] == 1) {
				for (int i = 0; i < col_length; i++) {
					MAtrix_out[row_up][i] = MAtrix_out[row_up][i] ^ MAtrix_out[row][i];
				}
			}
		}
	}
}

void CNBLDPC::H_GaussEliminate_bit(int**Matrix_in, int** MAtrix_out,int col_length,int row_length, int*order_seri) {

	//step1: Initial
	for (int i = 0; i < row_length; i++) {
		for (int j = 0; j < col_length; j++) {
			MAtrix_out[i][j] = Matrix_in[i][j];
		}
	}
	//step2:gauss eliminate
	for (int row = row_length - 1; row >= 0; row--)
	{
		//cout << "row:" << row << " ";
		int num_temp = row + col_length - row_length;
		int col = order_seri[num_temp];
		if (MAtrix_out[row][col] == 0)
		{
			bool exchanged = false;
			// search up
			if (!exchanged)
			{
				for (int row_up = row - 1; row_up >= 0; row_up--)
				{
					if (MAtrix_out[row_up][col] != 0)
					{
						for (int i = 0; i < col_length; i++) {
							MAtrix_out[row][i] = MAtrix_out[row][i]^ MAtrix_out[row_up][i];
						}
						exchanged = true;
						break;
					}
				}
			}
			// if cannot find nonzero ele, search left
			if (!exchanged)
			{
				row++;
				for (int i = num_temp - 1; i >= 0; i--)
				{
					int col_left = order_seri[i];
					order_seri[i] = order_seri[i+1];
					order_seri[i+1] = col_left;
					
				}
			}
		}
		// make the H[0:row-1][col] = 0
		for (int row_up = row - 1; row_up >= 0; row_up--)
		{
			//cout << row_up<<" ";
			if (MAtrix_out[row_up][col] != 0)
			{
				for (int i = 0; i < col_length; i++) {
					MAtrix_out[row_up][i] = MAtrix_out[row_up][i] ^ MAtrix_out[row][i];
				}
			}
		}
	}
	//
	for (int i = 0; i < row_length; i++) {
		int row_start = i + 1;
		for (int j = row_start; j < row_length; j++) {
			int pos = order_seri[col_length - row_length + i];
			if(MAtrix_out[j][pos]==1){
				for (int k = 0; k < col_length; k++) {
					MAtrix_out[j][k] = MAtrix_out[j][k] ^ MAtrix_out[i][k];
				}
			}			
		}
	}
}

void CNBLDPC::OSD_Encode_bit(int**Matrix_in, int*seq_in, int* seq_out, int*relyseri_bit, int col_length, int row_length) {
	
	int depend_bit_length = col_length - row_length;
	int temp;
	//step1 
	for (int i = 0; i < col_length; i++) {
		seq_out[i] = seq_in[i];
	}
	//step2 利用relyseri_bit与H_OSD_bit改变特定的比特位，使之满足初始的校验式。
	for (int i = 0; i < row_length; i++) {
		temp = 0;
		for (int j = 0; j <depend_bit_length; j++) {
			temp = temp ^ (Matrix_in[i][relyseri_bit[j]] * seq_out[relyseri_bit[j]]);

		}
		seq_out[relyseri_bit[depend_bit_length + i]] = temp;
	}
}




void CNBLDPC::compute_min_distance_bit() {
	int CodeLen_bit = CodeLen*log(GFq) / log(2);
	double temp = 0;
	double temp_ch = 0;
	//step1 计算OSD_output的欧式距
	
	for (int i = 0; i < CodeLen_bit; i++) {
		if (((L_ch_bit[i] < 0) && (OSD_temp_out_bit[i] == 1)) || ((L_ch_bit[i] > 0) && (OSD_temp_out_bit[i] == 0))) {
			temp_ch = temp_ch + fabs(L_ch_bit[i]);
		}
	}

	//step2 选择欧式距最小的码保存	
	if (temp_ch < min_distance) {
		for (int i = 0; i < CodeLen_bit; i++) {
			near_optimal_bit[i] = OSD_temp_out_bit[i];
		}
		min_distance = temp_ch;
	}
}

void CNBLDPC::MatrixH2G(int** H_matrix, int** G_matrix, int H_ChkLen, int H_CodeLen ) {
	//initial
	int *order_seri;
	int **H_IP;
	int MsgLen = H_CodeLen - H_ChkLen;
	order_seri = new int[H_CodeLen];
	H_IP = new int*[H_ChkLen];
	for (int i = 0; i < H_ChkLen; i++) { H_IP[i] = new int[H_CodeLen]; }
	for (int i = 0; i < MsgLen; i++) {
		for (int j = 0; j < H_CodeLen; j++) {
			G_matrix[i][j] = 0;
		}
	}
    //Generate H_IP
	for (int i = 0; i < H_CodeLen; i++) { order_seri[i] = i; }
	H_GaussEliminate_bit(H_matrix, H_IP, H_CodeLen, H_ChkLen,order_seri);
	
	//Generate G_matrix from H_IP
	for (int i = 0; i < MsgLen; i++) {
		G_matrix[i][i] = 1;
		for (int j = MsgLen; j < H_CodeLen; j++) {
			G_matrix[i][j] = H_IP[j-MsgLen][i];
		}
	}

    //delete
	delete[]order_seri;
	for (int i = 0; i < H_ChkLen; i++) { delete[]H_IP[i]; }
	delete[]H_IP;

}

void CNBLDPC::CRCMatrixGen() {
	int MsgLen_bit = CodeLen_bit - ChkLen_bit;
	switch (CRCLen)
	{
	case 8:
		for (int i = 0; i < MsgLen_bit - CRCLen; i++) {
			CRCMatrix_G[i][i] = CRCMatrix_G[i][i+1] = CRCMatrix_G[i][i+4] = CRCMatrix_G[i][i+5] = CRCMatrix_G[i][i+7] = CRCMatrix_G[i][i+8] = 1;
		}
		break;
	case 16:
		for (int i = 0; i < MsgLen_bit - CRCLen; i++) {
			CRCMatrix_G[i][i] = CRCMatrix_G[i][i + 4] = CRCMatrix_G[i][i + 11] = CRCMatrix_G[i][i + 16]  = 1;
		}
		break;
	case 24://24B
		for (int i = 0; i < MsgLen_bit - CRCLen; i++) {
			CRCMatrix_G[i][i] = CRCMatrix_G[i][i + 1] = CRCMatrix_G[i][i + 18] = CRCMatrix_G[i][i + 19] = CRCMatrix_G[i][i + 23] = CRCMatrix_G[i][i + 24] = 1;
			
		}
		break;
	default:
		break;
	}

	int* seri;
	seri = new int[MsgLen_bit];
	for (int i = 0; i < MsgLen_bit; i++) { seri[i] = i;}
	G_GaussEliminate_bit(CRCMatrix_G, CRCMatrix_G_IP, MsgLen_bit, MsgLen_bit - CRCLen, seri);

	for (int i = 0; i < CRCrate; i++) {
		CRCMatrix_partH_PI[i][i + MsgLen_bit - CRCrate] = 1;
		for (int j = 0; j < MsgLen_bit - CRCLen; j++) {
			CRCMatrix_partH_PI[i][j] = CRCMatrix_G_IP[j][MsgLen_bit - CRCrate + i];
		}		
	}

	delete[] seri;

}