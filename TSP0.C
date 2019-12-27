#define _CRT_SECURE_NO_WARNINGS
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<stdio.h>
#include<mpi.h>

#define FILE_PATH  "d:\\pcb442.tsp","r"  //数据文件名
#define NPROC 4
#define N_COLONY 100  // N_COLONY>=xColony
#define CITY     442  // CITY>=xCity
#define GRPSCALE 8//组的规模@@@@@@@@@@@@
#define MAXGEN 2000000
#define INTERVAL 10000
int     xColony = 100;     //##//  个体数
int     xCity = CITY;
double  probab1 = 0.02;    //##//  变异概率
long    NOCHANGE = 200000;  //##//  最大停止改变代数
long    maxGen = 200000;    //##//  停机代数
int     colony[N_COLONY * 2][CITY], colony2[N_COLONY][CITY]; //zhongqun		【200】【442】？

double  cityXY[CITY][2];
double  city_dis[CITY][CITY];
double  dis_p[N_COLONY * 2]; //适应值
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

//我的部分===================================
char filepath[100], filepath2[100], filepath3[100];





//=========================================
int main(int argc, char* argv[])
{



	//我的部分=====================================
	int i = 0;
	int k1, k2, l1, l2, pos_flag, icount, sendyet;
	int sm = 0;//by me to 标志插入与否
	int     tempcol[CITY];//传回的个体的染色体
	double  tempdis;//传回的个体的评估值
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
	double data[NPROC], data2[NPROC], tempdata, tempfk[N_COLONY], * ptrdata[NPROC], * ptrfk[NPROC], dt, databest[NPROC], * ptrtemp, * ptrfktemp;//指针数组用来存放精英的顺序
	int chrom[NPROC][CITY], tempchr[CITY], * ptrchrom[NPROC], chrombest[NPROC][CITY], * ptrctemp;//（问题相关）第一维是进程数，第二维比城市数大。
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
	 //我的部分=====================================

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mytid);/**//*得到当前进程号*/
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
		//初始化指针数组
		//初始化x
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
			printf("王兄，机器上有%d进程\n", i + 1);
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
					}	//收到各个子进程所得最优解
							//找出最好解
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
				
		
					
				//找出最好解
				for (i = 0; i < NPROC; i++)
				{
					MPI_Send(&data[mbest_i], 1, MPI_DOUBLE, i + 1, 99, MPI_COMM_WORLD);
					printf("发。。。");
					MPI_Send(chrom[mbest_i], CITY, MPI_INT, i + 1, 99, MPI_COMM_WORLD);
				}
				intr_best[times] = best;
				timeNowm = clock();
				period = (double)(timeNowm - timeStartm) / CLK_TCK; //CLK_TCK为clock()函数的时间单位，即时钟打点
				timesave[times] = period;
				//PDi0++;
			}//传送最优	
		}//mytid ==0
		else {
			if (loopcounter % INTERVAL == 0 && loopcounter != 0)//当迁移开始
			{
				//发
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
				//发送子种群最好解
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
				if (!sm)//和这个种群中所有都不同的
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
				if ((rand() / 32768.0) < probab1)     //内变异算子
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
			if (i >= xColony)//此处是报道+加代数 
			{
				select1();
				Ni++; GenNum++; i = 0;
	
			}
		
		}//子进程部分

	}
	////初始化结束，变异开始============================================
	//for (loopcounter = 0; loopcounter <= MAXGEN; loopcounter++) {
	//	if (loopcounter % INTERVAL == 0 && loopcounter != 0)
	//		times++;
	//	if (mytid == 0)//主进程
	//	{
	//		if (loopcounter % INTERVAL == 0 && loopcounter != 0)//当迁移开始
	//		{
	//			printf("Time(s)=%d,The generation is %ld", times, loopcounter);
	//			//接收各个子种群的最好解
	//			for (i = 0; i < NPROC; i++)
	//			{
	//				MPI_Recv(&data[i], 1, MPI_DOUBLE, i + 1, 99, MPI_COMM_WORLD, &status);
	//				MPI_Recv(chrom[i], CITY, MPI_INT, i + 1, 99, MPI_COMM_WORLD, &status);
	//			}
	//			//接收各个子种群的最好解
	//			//找出最好解
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
	//			//找出最好解
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
	//		if (loopcounter % INTERVAL == 0 && loopcounter != 0)//当迁移开始
	//		{
	//			//发
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
	//			//发送子种群最好解
	//			//接收群体最好解%%%%%%%%%%%%%%%%%%%%%%%%%%
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
	//			if (!sm)//和这个种群中所有都不同的
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
	//		for (j = 0; j < xCity; j++)temp[j] = colony[i][j];			//第0个族群分布的城市

	//		disChange = 0; pos_flag = 0;
	//		pos_C = rand() % xCity;			//随机产生一个数字，分割点位置
	//		for (;;)
	//		{
	//			if ((rand() / 32768.0) < probab1)     //内变异算子
	//			{
	//				do
	//					pos_C1 = rand() % xCity;			//分割点2
	//				while (pos_C1 == pos_C);
	//				C1 = colony[i][pos_C1];  //int型
	//			}
	//			else
	//			{
	//				do
	//					j = rand() % xColony;
	//				while (j == i);			//第i个与第j个种群·   j为随机产生
	//				k = position(colony[j], temp[pos_C]);	//temp暂存的就是第i个种群，在这里是第0个
	//				//因为temp暂存的是第i个城市序列
	//				//再由于随机产生的分割点为pos_C，在第i个城市序列的位置就是temp[pos_C]
	//				//position函数的作用就是，查找，在城市序列j，中，分割点的位置
	//				// pos_C =32 temp[32] = 11
	//				// 在colony[j=52]中， k=293， colony[52][293] = 11
	//				C1 = colony[j][(k + 1) % xCity];			//第j个种群的第k+1个城市，就是分割点之后一个

	//				//序列1   .......11 | 86........63|415........


	//				pos_C1 = position(temp, C1);		//反过来了，在第i个城市序列中，获得k+1的位置，命名为pos_C11111111111111
	//			}
	//			if ((pos_C + 1) % xCity == pos_C1 || (pos_C - 1 + xCity) % xCity == pos_C1)break;			//避免选到同一个点
	//			k1 = temp[pos_C];
	//			k2 = temp[(pos_C + 1) % xCity];			//这里就是第i个城市序列的第k个的后面一个
	//			l1 = temp[pos_C1];
	//			l2 = temp[(pos_C1 + 1) % xCity];
	//			disChange += city_dis[k1][l1] + city_dis[k2][l2] - city_dis[k1][k2] - city_dis[l1][l2];//变化 从11 到 63 +从86到415，   再减掉原来的
	//			invert(pos_C, pos_C1);			//pos 还是temp的下标志
	//			pos_flag++;
	//			if (pos_flag > xCity - 1)
	//				break;  //////////400多次？
	//			pos_C++;
	//			if (pos_C >= xCity)pos_C = 0;                 /**********************/
	//		}
	//		dis_p[N_COLONY + i] = dis_p[i] + disChange;			//每个i都有，不急
	//		disChange = 0;
	//		for (j = 0; j < xCity; j++)
	//			colony[N_COLONY + i][j] = temp[j];
	//		i++;
	//		if (i >= xColony)//此处是报道+加代数 
	//		{
	//			Ni++;

	//			sendyet = 0;
	//			select1();
	//			//Ni++; GenNum++; i = 0;
	//			//sumbest = dis_p[0];		//就只是定为第一个
	//			//for (j = 0; j < N_COLONY; j++)
	//			//	if (sumbest > dis_p[j])
	//			//		sumbest = dis_p[j];	//不记录路径
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

	//第几次 总值 平均分值 最大分值 最小分值 分值方差 当前最好解
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
		if (dis_p[N_COLONY + j] < dis_p[j])	//因为是对称的，所以直接加j即可
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
	//fclose(fp);			//读取文件获得city的xy

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
			if (j < i)  city_dis[i][j] = city_dis[j][i];			//求出每个城市之间的距离
		}

	mod = xCity;
	for (i = 0; i < xCity; i++)array[i] = i;     //    init colony[][]     
	for (i = 0; i < xColony; i++, mod = xCity)		//mod设置为city的数量，对于每个城市（442），每个种群（100）
		for (j = 0; j < xCity; j++)
		{
			sign = rand() % mod;
			colony[i][j] = array[sign];		//array就是0到442的随机值，第i个种族，第J个城市
			t = array[mod - 1];				//一个city被分配给了一个种族，减一，t为最后一个值
			array[mod - 1] = array[sign];		//array的最后一位，记录为随机的sign值
			array[sign] = t;						//交换了位置
			mod--;
			if (mod == 1) colony[i][++j] = array[0];
		}						//最后为，100个种族分散在442个城市的结果

	for (i = 0; i < xColony; i++)		    /*    init dis_p[]       */
	{
		dis_p[i] = 0;//100个种族适应值,也就是距离
		for (j = 0; j < xCity - 1; j++)
			dis_p[i] = dis_p[i] + city_dis[*(*(colony + i) + j)][*(*(colony + i) + j + 1)];
		dis_p[i] = dis_p[i] + city_dis[**(colony + i)][*(*(colony + i) + xCity - 1)];
	}

	ibest = 0; sumbest = dis_p[0];	    /*  init ibest & sumbest */
	sumTemp = sumbest * 5;

	Ni = 0;               /*   initialize GunNum & Ni    */
	printf("init success!!!\n");
}

void invert(int pos_start, int pos_end)			//两个分割点之间调换位置
{
	int j, k, t;
	if (pos_start < pos_end)
	{
		j = pos_start + 1; k = pos_end;
		for (; j <= k; j++, k--)
		{
			t = temp[j]; temp[j] = temp[k]; temp[k] = t;
			// t为最后一个城市
		 //从j到k的temp交换次序


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

