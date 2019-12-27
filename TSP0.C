#define _CRT_SECURE_NO_WARNINGS
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<stdio.h>
#include<mpi.h>

#define FILE_PATH  "d:\\pcb442.tsp","r"  //�����ļ���
#define NPROC 4
#define N_COLONY 100  // N_COLONY>=xColony
#define CITY     442  // CITY>=xCity
#define GRPSCALE 8//��Ĺ�ģ@@@@@@@@@@@@
#define MAXGEN 2000000
#define INTERVAL 10000
int     xColony = 100;     //##//  ������
int     xCity = CITY;
double  probab1 = 0.02;    //##//  �������
long    NOCHANGE = 200000;  //##//  ���ֹͣ�ı����
long    maxGen = 200000;    //##//  ͣ������
int     colony[N_COLONY * 2][CITY], colony2[N_COLONY][CITY]; //zhongqun		��200����442����

double  cityXY[CITY][2];
double  city_dis[CITY][CITY];
double  dis_p[N_COLONY * 2]; //��Ӧֵ
double  sumbest, sumTemp;
int     temp[CITY], ibest;
clock_t timeStart, timeNow, timeTemp;
long    GenNum, Ni;
void    init();
int     position(int* tmp, int C);
void    invert(int pos_start, int pos_end);
void    printBest(long GenNum);
void    tempTest(int i);
void    mapped();
void    LastCP();
double  path(int tmp[], int k1, int k2);
void select1();
void select2();
double SPAD_compute();
void initm();
FILE* fpme;

//�ҵĲ���===================================
char filepath[100], filepath2[100], filepath3[100];





//=========================================
int main(int argc, char* argv[])
{



	//�ҵĲ���=====================================
	int i = 0;
	int k1, k2, l1, l2, pos_flag, icount, sendyet;
	int sm = 0;//by me to ��־�������
	int     tempcol[CITY];//���صĸ����Ⱦɫ��
	double  tempdis;//���صĸ��������ֵ
	int sign = 1, tids[NPROC], ccs[NPROC], n, nproc,  who, msgtype, i1, j1, tempi = 0;
	int ibest2, iworst, ipass;
	int mbest_i=0, PDi = 0, PDi0 = 0, intr_best[25000];
	double timesave[25000];
	double period=0;
	int times = 0;
	int loopcounter = 0;
	char mname[30];
	int mytid, numprocs;
	int r = 0;
	//int nproc = 0;
	char nameofm[NPROC][30];
	double data[NPROC], data2[NPROC], tempdata, tempfk[N_COLONY], * ptrdata[NPROC], * ptrfk[NPROC], dt, databest[NPROC], * ptrtemp, * ptrfktemp;//ָ������������ž�Ӣ��˳��
	int chrom[NPROC][CITY], tempchr[CITY], * ptrchrom[NPROC], chrombest[NPROC][CITY], * ptrctemp;//��������أ���һά�ǽ��������ڶ�ά�ȳ�������
	double temps = 0.0, result, fk[NPROC][N_COLONY];
	double best;
	int hao[NPROC];
	clock_t timeStartm, timeNowm, timeTempm;
	time_t Start, End;
	double  speed;

	MPI_Status status;
	FILE* fppp;
	int back_data = 0;
	//=========================================
	register int C1, j, k, pos_C, pos_C1; int mem;
	register double disChange;

	// timeStart=timeNow=timeTemp=clock();  
	// init();
	 //�ҵĲ���=====================================

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mytid);/**//*�õ���ǰ���̺�*/
	MPI_Get_processor_name(mname, &r);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	if (mytid == 0) printf("number of processes: %d\n...", numprocs);
	printf("%s: Hello world from process %d \n", mname, mytid);
	srand((unsigned)time(NULL) + mytid);
	nproc = NPROC;
	//if ((fppp = fopen(FILE_PATH)) == NULL)exit(0);


	


	if (mytid == 0)
	{
			for (i = 0; i < NPROC; i++)
		{
			ptrdata[i] = &data[i];
			ptrchrom[i] = chrom[i];
			ptrfk[i] = fk[i];
		}
		//��ʼ��ָ������
		//��ʼ��x
		for (i = 0; i < NPROC; i++)
		{
			hao[i] = i;
		}
		timeStartm = timeNowm = timeTempm = clock();
		time(&Start);
		initm();
		for (i = 0; i < NPROC; i++)
		{
			MPI_Send(*cityXY, CITY * 2, MPI_DOUBLE, i + 1, 99, MPI_COMM_WORLD);////////
			printf("���֣���������%d����\n", i + 1);
		}
		//fclose(fppp);
	}
	else
	{
		MPI_Recv(*cityXY, CITY * 2, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
		init();
		printf("hahah");
	}

	for (loopcounter = 0; loopcounter < MAXGEN; loopcounter++)
	{
		if (loopcounter % INTERVAL == 0 && loopcounter != 0)
		{
			times++;
		}
		if (mytid == 0) {
			if (loopcounter % INTERVAL == 0 && loopcounter != 0)
			{
				printf("Time(s)=%d,The generation is %ld\n", times, loopcounter);
				for (i= 0; i < NPROC; i++){
					
						MPI_Recv(&data[i], 1, MPI_DOUBLE, i + 1, 99, MPI_COMM_WORLD, &status);
						MPI_Recv(chrom[i], CITY, MPI_INT, i + 1, 99, MPI_COMM_WORLD, &status);
					}	//�յ������ӽ����������Ž�
							//�ҳ���ý�
				best = data[0];
				mbest_i = 0;
			
				for (back_data = 0; back_data < nproc; back_data++)
				{
					//printf("best = %Lf\n",best);
					if (data[back_data] < best)
					{
						best = data[back_data];
						mbest_i = back_data;
					}
				}
				
		
					
				//�ҳ���ý�
				for (i = 0; i < NPROC; i++)
				{
					MPI_Send(&data[mbest_i], 1, MPI_DOUBLE, i + 1, 99, MPI_COMM_WORLD);
					printf("��������");
					MPI_Send(chrom[mbest_i], CITY, MPI_INT, i + 1, 99, MPI_COMM_WORLD);
				}
				intr_best[times] = best;
				timeNowm = clock();
				period = (double)(timeNowm - timeStartm) / CLK_TCK; //CLK_TCKΪclock()������ʱ�䵥λ����ʱ�Ӵ��
				timesave[times] = period;
				//PDi0++;
			}//��������	
		}//mytid ==0
		else {
			if (loopcounter % INTERVAL == 0 && loopcounter != 0)//��Ǩ�ƿ�ʼ
			{
				//��
				iworst = 0;
				ibest2 = 0;
				ipass = rand() % xColony;
				for (j1 = 1; j1 < xColony; j1++)
				{
					if (dis_p[j1] > dis_p[iworst])
						iworst = j1;
					if (dis_p[j1] < dis_p[ibest2])
						ibest2 = j1;
				}
				MPI_Send(&dis_p[ibest2], 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
				//printf("disp[o] =%d\n", dis_p[ibest2]);

				MPI_Send(colony[ibest2], CITY, MPI_INT, 0, 99, MPI_COMM_WORLD);
				//��������Ⱥ��ý�
				MPI_Recv(&tempdis, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);

				MPI_Recv(tempcol, CITY, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
				for (j1 = 1; j1 < xColony; j1++)
				{
					if (dis_p[j1] == tempdis)
					{
						sm = 1;

						break;
					}
				}
				if (!sm)//�������Ⱥ�����ж���ͬ��
				{
					dis_p[iworst] = tempdis;
					for (j1 = 0; j1 < xCity; j1++)
					{
						colony[iworst][j1] = tempcol[j1];

					}
					sendyet = 1;
					sm = 0;
				}

			
			}
			for (j = 0; j < xCity; j++)temp[j] = colony[i][j];
			disChange = 0; pos_flag = 0;
			pos_C = rand() % xCity;
			for (;;)
			{
				if ((rand() / 32768.0) < probab1)     //�ڱ�������
				{
					do
						pos_C1 = rand() % xCity;
					while (pos_C1 == pos_C);
					C1 = colony[i][pos_C1];
				}
				else
				{
					do
						j = rand() % xColony;
					while (j == i);
					k = position(colony[j], temp[pos_C]);
					C1 = colony[j][(k + 1) % xCity];
					pos_C1 = position(temp, C1);
				}
				if ((pos_C + 1) % xCity == pos_C1 || (pos_C - 1 + xCity) % xCity == pos_C1)break;
				k1 = temp[pos_C];
				k2 = temp[(pos_C + 1) % xCity];
				l1 = temp[pos_C1];
				l2 = temp[(pos_C1 + 1) % xCity];
				disChange += city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];
				invert(pos_C, pos_C1);
				pos_flag++;
				if (pos_flag > xCity - 1)
					break;  ////////////
				pos_C++;
				if (pos_C >= xCity)pos_C = 0;                 /**********************/
			}
			dis_p[N_COLONY + i] = dis_p[i] + disChange;
			disChange = 0;
			for (j = 0; j < xCity; j++)
				colony[N_COLONY + i][j] = temp[j];
			i++;
			if (i >= xColony)//�˴��Ǳ���+�Ӵ��� 
			{
				select1();
				Ni++; GenNum++; i = 0;
	
			}
		
		}//�ӽ��̲���

	}
	////��ʼ�����������쿪ʼ============================================
	//for (loopcounter = 0; loopcounter <= MAXGEN; loopcounter++) {
	//	if (loopcounter % INTERVAL == 0 && loopcounter != 0)
	//		times++;
	//	if (mytid == 0)//������
	//	{
	//		if (loopcounter % INTERVAL == 0 && loopcounter != 0)//��Ǩ�ƿ�ʼ
	//		{
	//			printf("Time(s)=%d,The generation is %ld", times, loopcounter);
	//			//���ո�������Ⱥ����ý�
	//			for (i = 0; i < NPROC; i++)
	//			{
	//				MPI_Recv(&data[i], 1, MPI_DOUBLE, i + 1, 99, MPI_COMM_WORLD, &status);
	//				MPI_Recv(chrom[i], CITY, MPI_INT, i + 1, 99, MPI_COMM_WORLD, &status);
	//			}
	//			//���ո�������Ⱥ����ý�
	//			//�ҳ���ý�
	//			best = data[0];
	//			for (i = 0; i < NPROC; i++)
	//			{
	//				if (data[i] < best)
	//				{
	//					best = data[i];
	//					mbest_i = i;
	//				}
	//			}
	//			intr_best[times] = (int)best;
	//			//�ҳ���ý�
	//			for (i = 0; i < NPROC; i++)
	//			{
	//				MPI_Send(&data[mbest_i], 1, MPI_DOUBLE, i + 1, 99, MPI_COMM_WORLD);
	//				MPI_Send(chrom[mbest_i], CITY, MPI_INT, i + 1, 99, MPI_COMM_WORLD);
	//			}
	//			//PDi0++;

	//		}
	//	}
	//	else
	//	{
	//		if (loopcounter % INTERVAL == 0 && loopcounter != 0)//��Ǩ�ƿ�ʼ
	//		{
	//			//��
	//			iworst = 0;
	//			ibest2 = 0;
	//			ipass = rand() % xColony;
	//			for (j1 = 1; j1 < xColony; j1++)
	//			{
	//				if (dis_p[j1] > dis_p[iworst])
	//					iworst = j1;
	//				if (dis_p[j1] < dis_p[ibest2])
	//					ibest2 = j1;
	//			}
	//			MPI_Send(&dis_p[ibest2], 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
	//			MPI_Send(colony[ibest2], CITY, MPI_INT, 0, 99, MPI_COMM_WORLD);
	//			//��������Ⱥ��ý�
	//			//����Ⱥ����ý�%%%%%%%%%%%%%%%%%%%%%%%%%%
	//			MPI_Recv(&tempdis, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
	//			MPI_Recv(tempcol, CITY, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
	//			for (j1 = 1; j1 < xColony; j1++)
	//			{
	//				if (dis_p[j1] == tempdis)
	//				{
	//					sm = 1;
	//					break;
	//				}
	//			}
	//			if (!sm)//�������Ⱥ�����ж���ͬ��
	//			{
	//				dis_p[iworst] = tempdis;
	//				for (j1 = 0; j1 < xCity; j1++)
	//				{
	//					colony[iworst][j1] = tempcol[j1];
	//				}
	//				sendyet = 1;
	//				sm = 0;
	//			}
	//		}
	//		for (j = 0; j < xCity; j++)temp[j] = colony[i][j];			//��0����Ⱥ�ֲ��ĳ���

	//		disChange = 0; pos_flag = 0;
	//		pos_C = rand() % xCity;			//�������һ�����֣��ָ��λ��
	//		for (;;)
	//		{
	//			if ((rand() / 32768.0) < probab1)     //�ڱ�������
	//			{
	//				do
	//					pos_C1 = rand() % xCity;			//�ָ��2
	//				while (pos_C1 == pos_C);
	//				C1 = colony[i][pos_C1];  //int��
	//			}
	//			else
	//			{
	//				do
	//					j = rand() % xColony;
	//				while (j == i);			//��i�����j����Ⱥ��   jΪ�������
	//				k = position(colony[j], temp[pos_C]);	//temp�ݴ�ľ��ǵ�i����Ⱥ���������ǵ�0��
	//				//��Ϊtemp�ݴ���ǵ�i����������
	//				//��������������ķָ��Ϊpos_C���ڵ�i���������е�λ�þ���temp[pos_C]
	//				//position���������þ��ǣ����ң��ڳ�������j���У��ָ���λ��
	//				// pos_C =32 temp[32] = 11
	//				// ��colony[j=52]�У� k=293�� colony[52][293] = 11
	//				C1 = colony[j][(k + 1) % xCity];			//��j����Ⱥ�ĵ�k+1�����У����Ƿָ��֮��һ��

	//				//����1   .......11 | 86........63|415........


	//				pos_C1 = position(temp, C1);		//�������ˣ��ڵ�i�����������У����k+1��λ�ã�����Ϊpos_C11111111111111
	//			}
	//			if ((pos_C + 1) % xCity == pos_C1 || (pos_C - 1 + xCity) % xCity == pos_C1)break;			//����ѡ��ͬһ����
	//			k1 = temp[pos_C];
	//			k2 = temp[(pos_C + 1) % xCity];			//������ǵ�i���������еĵ�k���ĺ���һ��
	//			l1 = temp[pos_C1];
	//			l2 = temp[(pos_C1 + 1) % xCity];
	//			disChange += city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];//�仯 ��11 �� 63 +��86��415��   �ټ���ԭ����
	//			invert(pos_C, pos_C1);			//pos ����temp���±�־
	//			pos_flag++;
	//			if (pos_flag > xCity - 1)
	//				break;  //////////400��Σ�
	//			pos_C++;
	//			if (pos_C >= xCity)pos_C = 0;                 /**********************/
	//		}
	//		dis_p[N_COLONY + i] = dis_p[i] + disChange;			//ÿ��i���У�����
	//		disChange = 0;
	//		for (j = 0; j < xCity; j++)
	//			colony[N_COLONY + i][j] = temp[j];
	//		i++;
	//		if (i >= xColony)//�˴��Ǳ���+�Ӵ��� 
	//		{
	//			Ni++;

	//			sendyet = 0;
	//			select1();
	//			//Ni++; GenNum++; i = 0;
	//			//sumbest = dis_p[0];		//��ֻ�Ƕ�Ϊ��һ��
	//			//for (j = 0; j < N_COLONY; j++)
	//			//	if (sumbest > dis_p[j])
	//			//		sumbest = dis_p[j];	//����¼·��
	//			//printf("%d:%f\n", GenNum, sumbest);
	//			//if (GenNum % 2000 == 0 && GenNum < maxGen)
	//			//	printBest(GenNum);
	//			//if (GenNum >= maxGen)
	//			//{
	//			//	timeNow = clock();
	//			//	printf("Finial solution:");
	//			//	printBest(GenNum);
	//			//	exit(1);
	//			//}
	//		}


	//	}
	//}
	if (mytid == 0)
	{
		if ((fpme = fopen("tsp0.txt", "a")) == NULL)exit(0);
		fprintf(fpme, "This is a result of %d:%lf\n", CITY, best);

		for (int i = 0; i < times; i++)
		{
			fprintf(fpme, "%d %f %d\n", i,timesave[i],intr_best[i]);
		}

		//fprintf(fpme,"PD%d:%f\n",xx,PD[xx]);

	//�ڼ��� ��ֵ ƽ����ֵ ����ֵ ��С��ֵ ��ֵ���� ��ǰ��ý�
		fclose(fpme);
	}
	MPI_Finalize();
	return 0;
}



//===================================================

void initm() {

	int i;
	double x, y;
	double d;
	FILE* fp;
	srand((unsigned)time(NULL));

	if ((fp = fopen(FILE_PATH)) == NULL)
	{
		MPI_Finalize();
		exit(0);
	}
	fscanf(fp, "%d", &xCity);
	printf("%d", xCity);
	for (i = 0; i < xCity; i++)      /*  init cityXY[][]  */
	{
		fscanf(fp, "%*d%Lf%Lf", &x, &y);
		cityXY[i][0] = x;
		cityXY[i][1] = y;
		printf(" line = %d,x =%Lf,y = %Lf\n", i, x, y);
	}
	fclose(fp);

}
//  for(;;)
//  {
//  }
// 
// return 1;
//}
void select1()
{
	int j, k;
	for (j = 0; j < N_COLONY; j++)
		if (dis_p[N_COLONY + j] < dis_p[j])	//��Ϊ�ǶԳƵģ�����ֱ�Ӽ�j����
		{
			dis_p[j] = dis_p[N_COLONY + j];
			for (k = 0; k < CITY; k++)
				colony[j][k] = colony[N_COLONY + j][k];
		}
}
void init()
{
	int i, j, t, sign, mod, array[CITY];
	double x, y;
	double d;
	//FILE *fp;
	//srand( (unsigned)time( NULL ) );
	//if((fp=fopen(FILE_PATH))==NULL)exit(0);
	//fscanf(fp,"%d",&xCity);
	//for(i=0;i<xCity;i++)      /*  init cityXY[][]  */
	//{ fscanf(fp,"%*d%Lf%Lf",&x,&y);
	//  cityXY[i][0]=x;
	//  cityXY[i][1]=y;
	//}
	//fclose(fp);			//��ȡ�ļ����city��xy

	for (i = 0; i < xCity; i++)    /*  init city_dis[] */
		for (j = 0; j < xCity; j++)
		{
			if (j > i)
			{
				d = (cityXY[i][0] - cityXY[j][0]) * (cityXY[i][0] - cityXY[j][0]) * 1.0 +
					(cityXY[i][1] - cityXY[j][1]) * (cityXY[i][1] - cityXY[j][1]) * 1.0;
				city_dis[i][j] = (int)(sqrt(d) + 0.5);
				continue;
			}
			if (j == i) { city_dis[i][j] = 0; continue; }
			if (j < i)  city_dis[i][j] = city_dis[j][i];			//���ÿ������֮��ľ���
		}

	mod = xCity;
	for (i = 0; i < xCity; i++)array[i] = i;     //    init colony[][]     
	for (i = 0; i < xColony; i++, mod = xCity)		//mod����Ϊcity������������ÿ�����У�442����ÿ����Ⱥ��100��
		for (j = 0; j < xCity; j++)
		{
			sign = rand() % mod;
			colony[i][j] = array[sign];		//array����0��442�����ֵ����i�����壬��J������
			t = array[mod - 1];				//һ��city���������һ�����壬��һ��tΪ���һ��ֵ
			array[mod - 1] = array[sign];		//array�����һλ����¼Ϊ�����signֵ
			array[sign] = t;						//������λ��
			mod--;
			if (mod == 1) colony[i][++j] = array[0];
		}						//���Ϊ��100�������ɢ��442�����еĽ��

	for (i = 0; i < xColony; i++)		    /*    init dis_p[]       */
	{
		dis_p[i] = 0;//100��������Ӧֵ,Ҳ���Ǿ���
		for (j = 0; j < xCity - 1; j++)
			dis_p[i] = dis_p[i] + city_dis[*(*(colony + i) + j)][*(*(colony + i) + j + 1)];
		dis_p[i] = dis_p[i] + city_dis[**(colony + i)][*(*(colony + i) + xCity - 1)];
	}

	ibest = 0; sumbest = dis_p[0];	    /*  init ibest & sumbest */
	sumTemp = sumbest * 5;

	Ni = 0;               /*   initialize GunNum & Ni    */
	printf("init success!!!\n");
}

void invert(int pos_start, int pos_end)			//�����ָ��֮�����λ��
{
	int j, k, t;
	if (pos_start < pos_end)
	{
		j = pos_start + 1; k = pos_end;
		for (; j <= k; j++, k--)
		{
			t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			// tΪ���һ������
		 //��j��k��temp��������


		}
	}
	else
	{
		if (xCity - 1 - pos_start <= pos_end + 1)
		{
			j = pos_end; k = pos_start + 1;
			for (; k < xCity; j--, k++)
			{
				t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			}
			k = 0;
			for (; k <= j; k++, j--)
			{
				t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			}
		}
		else
		{
			j = pos_end; k = pos_start + 1;
			for (; j >= 0; j--, k++)
			{
				t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			}
			j = xCity - 1;
			for (; k <= j; k++, j--)
			{
				t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			}
		}
	}
}
int position(int* tmp, int C)
{
	int j;
	for (j = 0; j < xCity; j++)
		if (*(tmp + j) == C)break;
	return(j);
}
void printBest(long GenNum)
{
	int i;
	if ((fpme = fopen("tsp0.txt", "a")) == NULL)exit(0);
	//fprintf(fpme,"\n   CITY      %d\t\tN_COLONY  %d",CITY,N_COLONY);
	//fprintf(fpme,"\ntime     %4.2f",(double)(timeNow-timeStart)/CLOCKS_PER_SEC);
	//fprintf(fpme,"\n   distance  %f",sumbest);
	//fprintf(fpme,"\n   GenNum    %d\n\n",GenNum);
	fprintf(fpme, "%d\t%4.2f\t%d\n", GenNum, (double)(timeNow - timeStart) / CLOCKS_PER_SEC, (int)sumbest);
	fclose(fpme);

}
double path(int tmp[], int k1, int k2)
{
	int j, t1, t2; double temp_dis = 0;
	if (k2 > k1)
		for (j = k1; j < k2; j++)
			temp_dis += city_dis[tmp[j]][tmp[j + 1]];
	else
		for (j = k1; j < k2 + xCity; j++)
		{
			t1 = j % xCity; t2 = (j + 1) % xCity;
			temp_dis += city_dis[tmp[t1]][tmp[t2]];
		}
	return temp_dis;
}
