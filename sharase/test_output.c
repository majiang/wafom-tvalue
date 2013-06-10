#include <stdio.h>
#include <stdlib.h>

#define S_MAX 20
#define R_MAX 28
#define R 32
#define WORD_SIZE 32

unsigned int pow_2[WORD_SIZE+1] = {0x0,
		0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000, 0x4000000, 0x2000000, 0x1000000,
		0x800000, 0x400000, 0x200000, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000,
		0x8000, 0x4000, 0x2000, 0x1000, 0x800, 0x400, 0x200, 0x100,
		0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1};

void output_gray(int S, int num, unsigned int P[], unsigned int C[][R+1])
{
	int j = 1, s;
	unsigned int bitmask = 1;
	
	while(((num - 1) & bitmask) != 0)
	{
		bitmask = (bitmask << 1);
		j++;
	}
	
	for(s = 1; s <= S; s++)
	{
		P[s] = P[s] ^ C[s][j];
	}
}		

void construct_C(int S, unsigned int C[][R+1])
{
	int i, j, k, l;
	unsigned int temp_col;
	FILE *f1;
	unsigned int hex[600];
	unsigned int bitmask1;
	char fn[256];
	
	sprintf(fn, "harase_gen_mat_%d.txt", S);
	/* ファイルからの入力 */
	f1=fopen(fn, "r");
	for(i = 1; i <= S * R_MAX; i++)
	fscanf(f1, "%x", &hex[i]);
	fclose(f1);
	
	for(i = 1; i <= S; i++)
		for(j = 1; j <= R; j++) C[i][j] = 0;
		
	
	/* 逆転させて列ベクトルを作る */
	k = 1;
	for(i = 1; i <= S; i++)
	{
		for(j = 1; j <= R_MAX; j++)
		{
			C[i][j] = hex[k];
			k++;
		}
	}
}

int main(void)
{
	int i, j;
	int S = 10;
	
	unsigned int C[S_MAX+1][R+1];
	unsigned int P[S_MAX+1];
	
	construct_C(S, C);
	
	/*生成行列の出力*/
	for(i = 1; i <= S; i++)
	{
		for(j = 1; j <= R_MAX; j++)
		{
			printf("%x ", C[i][j]);
		}
		printf("\n");
	}
	
	/*点集合の出力*/
	for(j = 0; j <= S_MAX; j++) P[j] = 0;
	for(j = 1; j <= S; j++) printf("%10.10f ", P[j]*(1.0/4294967295.0));
	printf("\n");
	
	for(i = 1; i < 100; i++)
	{
		output_gray(S, i, P, C);
		for(j = 1; j <= S; j++) printf("%10.10f ", P[j]*(1.0/4294967295.0));
		printf("\n");
	}
	
	return 0;
}