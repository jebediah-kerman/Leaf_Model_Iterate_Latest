//Source code for spatial variability of the growth in cell area

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (void){
//Restore the data to arrays
	//Col 6:Wild type flower number 6
	type=6;
	//time:Time frames N:Number of mesh points CellNout:The label for the outside of the sepal
	int type,time[10],N[10],CellNout[10];
	time[0]=0;time[1]=12;time[2]=24;time[3]=36;time[4]=48;
	N[0]=224081;N[1]=192757;N[2]=275129;N[3]=320182;N[4]=414384;
//loop of time frames
	int t;
	for (t=1; t<=4; t++) { 
	//Read the data files,Create the output file
		int i,j,k,l;
		FILE *file1,*file11,*file2,*file22,*file111,*file222,*file3;
		char filename1[100],filename11[100],filename2[100],filename22[100],filename111[100],filename222[100],filename3[100];
		sprintf(filename1,"%d_%d_mesh.Int.rot.txt",type,time[t]); 
		sprintf(filename11,"%d_%d_mesh.near.txt",type,time[t]);
		sprintf(filename22,"%d_%d-%d.dat",type,time[t-1],time[t]);
		sprintf(filename111,"%d_%d_area.dat",type,time[t-1]);
		sprintf(filename222,"%d_%d_area.dat",type,time[t]);
		sprintf(filename3,"%d_%d_varea.dat",type,time[t]);
		file1 = fopen(filename1,"r");		//reading the mesh intensity file
		file11 = fopen(filename11,"r");		//reading the mesh neighbor file
		file22 = fopen(filename22,"r");		//reading lineage file
		file111 = fopen(filename111,"r");	//reading area at the current time frame
		file222 = fopen(filename222,"r");	//reading area at the next time frame
		file3 = fopen(filename3,"w");		//writing the spatial variability of areal growh file for each mesh
		//reading the mesh intensity file,reading the mesh neighbor file
		double *x,*y,*z;
		int sute,*MTnum,*Cnum,**MTnear;
		x = (double*)calloc(N[t]+10, sizeof(double));
		y = (double*)calloc(N[t]+10, sizeof(double));
		z = (double*)calloc(N[t]+10, sizeof(double));
		Cnum = (int*)calloc(N[t]+10, sizeof(double));
		MTnear = (int**)calloc(20, sizeof(double));
		for (i=1; i<=16; i++)MTnear[i] = (int*)calloc(N[t]+10, sizeof(double));
		for (i=0; i<N[t]; i++) {
			fscanf(file1,"%d%lf%lf%lf%d\n",&MTnum[i],&x[i],&y[i],&z[i],&Cnum[i]);
			fscanf(file11,"%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d\n",&sute,&sute2,&MTnear[1][i],&MTnear[2][i],&MTnear[3][i],&MTnear[4][i],&MTnear[5][i],&MTnear[6][i],&MTnear[7][i],&MTnear[8][i],&MTnear[9][i],&MTnear[10][i],&MTnear[11][i],&MTnear[12][i],&MTnear[13][i],&MTnear[14][i],&MTnear[15][i],&MTnear[16][i]);
		}
		//reading lineage file
		int *sun,*parent;
		sun = (int*)calloc(1000, sizeof(double));
		parent = (int*)calloc(1000, sizeof(double));
		for (i=0; i<700; i++) {
			fscanf(file22,"%d%d\n",&sun[i],&parent[i]);
		}
		//reading area at the current time frame
		int *Cnumareacur;
		double *areacur;
		Cnumareacur = (int*)calloc(10000+10, sizeof(double));
		areacur = (double*)calloc(10000+10, sizeof(double));
		for (i=0; i<700; i++) {
			fscanf(file111,"%d%lf\n",&Cnumareacur[i],&areacur[i]);
		}
		//reading area at the next time frame
		int *Cnumareanext;
		double *areanext;
		Cnumareanext = (int*)calloc(10000+10, sizeof(double));
		areanext = (double*)calloc(10000+10, sizeof(double));
		for (i=0; i<700; i++) {
			fscanf(file222,"%d%lf\n",&Cnumareanext[i],&areanext[i]);
		}
		//rearrange the order of two area files which always are consistent with next frame
		int *Cnumareacur1,*Cnumareanext1;
		double *areacur1,*areanext1;
		Cnumareacur1 = (int*)calloc(10000+10, sizeof(double));
		Cnumareanext1 = (int*)calloc(10000+10, sizeof(double));
		areacur1 = (double*)calloc(10000+10, sizeof(double));
		areanext1 = (double*)calloc(10000+10, sizeof(double));
		double sumaa;
		for (i=0; i<700; i++) {
			sumaa=0.0;
			for (j=0; j<700; j++) {
				if(sun[j]==Cnumareanext[i] && Cnumareanext[i]>0){
					sumaa=sumaa+areanext[i];
					Cnumareanext1[i]=sun[j];
				}
				for (k=0; k<700; k++) {
					if(parent[j]==Cnumareacur[k] && Cnumareacur[k]>0){
						areacur1[i]=areacur[k];
						Cnumareacur1[i]=parent[j];
					}
				}
			}
			areanext1[i]=sumaa;
		}
		fclose(file1);
		fclose(file11);
		fclose(file2);
		fclose(file22);
		fclose(file111);
		fclose(file222);
//Extract the outline of cells in order to detect the neighboring cells
		int count1,**MTnum1,*count3,**walllabel,**adjcell1,**adjcell2,**adjcell3;
		double **MTx1,**MTy1,**MTz1;
		int **sun1;
		MTx1 = (double**)calloc(10000, sizeof(double));
		MTy1 = (double**)calloc(10000, sizeof(double));
		MTz1 = (double**)calloc(10000, sizeof(double));
		MTnum1 = (int**)calloc(10000, sizeof(double));
		count3 = (int*)calloc(10000, sizeof(double));
		walllabel = (int**)calloc(10000, sizeof(double));
		adjcell1 = (int**)calloc(10000, sizeof(double));
		adjcell2 = (int**)calloc(10000, sizeof(double));
		adjcell3 = (int**)calloc(10000, sizeof(double));
		sun1 = (int**)calloc(10000, sizeof(double));
		for (i=0; i<10000; i++) {
			MTx1[i] = (double*)calloc(5000+10, sizeof(double));
			MTy1[i] = (double*)calloc(5000+10, sizeof(double));
			MTz1[i] = (double*)calloc(5000+10, sizeof(double));
			MTnum1[i] = (int*)calloc(5000+10, sizeof(double));
			walllabel[i] = (int*)calloc(5000, sizeof(double));
			adjcell1[i] = (int*)calloc(5000, sizeof(double));
			adjcell2[i] = (int*)calloc(5000, sizeof(double));
			adjcell3[i] = (int*)calloc(5000, sizeof(double));
			sun1[i] = (int*)calloc(1000, sizeof(double));
		}
		int countsun;
		for (j=1; j<10000; j++) {
			countsun=0;
			for (i=0; i<700; i++) {
				if(parent[i]==j){
					sun1[j][countsun]=sun[i];
					countsun=countsun+1;
				}
				if(sun[i]==0)break;
			}
		}
		int count1,j1;
		//Numbering the sons
		for (k=1; k<10000; k++) {
			count1=0;
			if(k==CellNout[zzz-1])continue;
			for (j1=0; j1<10; j1++) {
				j=sun1[k][j1];
				if(j==0)break;
				for (i=0; i<N[zzz]; i++) {
					if(Cnum[i]==-1){
						for (l=1; l<=16; l++) {
							if(j==Cnum[MTnear[l][i]] && MTnear[l][i]>0){
								MTx1[k][count1]=x1[i];
								MTy1[k][count1]=y1[i];
								MTz1[k][count1]=z1[i];
								MTnum1[k][count1]=MTnum[i];
								count1=count1+1;
							}
							if(j==Cnum[MTnear[l][i]] && MTnear[l][i]>0)break;
						}
					}
				}
			}
			count3[k]=count1;
		}
//Detection of the neighboring cells
		int *nearmesh,countmesh,cell1,cell2,cell3;
		int **nearcell,**nearcell1,countnear;
		nearcell = (int**)calloc(10000, sizeof(double));
		nearcell1 = (int**)calloc(10000, sizeof(double));
		for(i=0;i<10000;i++){
			nearcell[i]= (int*)calloc(50000, sizeof(double));
			nearcell1[i]= (int*)calloc(50000, sizeof(double));
		}
		for (j=1; j<10000; j++) {
			if(count3[j]==0 || j==CellNout[t-1])continue;
			countnear=0;
			for (i=0; i<=count3[j]; i++) {
				nearmesh = (int*)calloc(15, sizeof(double));
				countmesh=0;
				cell1=0;cell2=0;cell3=0;
				j1=MTnum1[j][i];
				for (l=1; l<=16; l++) {
					if(MTnear[l][j1]>0)nearmesh[l-1]=Cnum[MTnear[l][j1]];
				}
				for (k=0; k<16; k++){
					if(nearmesh[k]==-1)countmesh=countmesh+1;
					if(nearmesh[k]!=-1 && nearmesh[k]!=0 && cell1==0){
						cell1=nearmesh[k];
					}
					if(nearmesh[k]!=-1 && nearmesh[k]!=0 && cell1!=0 && cell2==0 && cell1!=nearmesh[k]){
						cell2=nearmesh[k];
					}
					if(nearmesh[k]!=-1 && nearmesh[k]!=0 && cell1!=0 && cell2!=0 && cell3==0 && cell1!=nearmesh[k] && cell2!=nearmesh[k]){
						cell3=nearmesh[k];
					}
				}
				adjcell1[j][i]=cell1;
				adjcell2[j][i]=cell2;
				adjcell3[j][i]=cell3;
				walllabel[j][i]=countmesh;
				if(countmesh==3){
					nearcell[j][countnear]=cell1;
					countnear=countnear+1;
					nearcell[j][countnear]=cell2;
					countnear=countnear+1;
					nearcell[j][countnear]=cell3;
					countnear=countnear+1;
				}
				free(nearmesh);
			}
		}
		int countnear1,*countnear2,i1,npa;
		countnear2 = (int*)calloc(10000, sizeof(double));
		for (j=1; j<10000; j++) {
			if(count3[j]==0 || j==CellNout[t-1])continue;
			countnear1=0;
			for (k=1; k<10000; k++) {
				if(count3[k]==0 || k==CellNout[t-1])continue;
				for (i=0; i<100; i++) {
					for (i1=0; i1<700; i1++) {
						if(nearcell[j][i]==sun[i1])npa=parent[i1];
					}
					if(k==npa && npa!=0){
						nearcell1[j][countnear1]=npa;
						countnear1=countnear1+1;
					}
					if(k==npa && npa!=0)break;
				}
			}
			countnear2[j]=countnear1;
		}
//Calculate spatial varibility in the growth of cell area
		int countarea;
		double sumarea,area0,areanext0,*areaN,*areanextN;
		double *varea,*p_varea,thre_varea,bin1_varea,cast_varea;
		int count_vrate=0,count_vani=0,count_varea=0;
		double bin_varea=0.3/50;
		varea = (double*)calloc(10000+10, sizeof(double));
		p_varea = (double*)calloc(500000, sizeof(double));
		for (j=1; j<10000; j++) {
			if(count3[j]==0 || j==CellNout[t-1])continue;
			areaN = (double*)calloc(30, sizeof(double));
			areanextN = (double*)calloc(30, sizeof(double));
			sumarea=0.0;
			for (i=1; i<700; i++) {
				if(j==Cnumareanext1[i]){
					area0=areanext1[i]/areacur1[i];
				}
			}
			countarea=0;
			for (k=0; k<countnear2[j]; k++) {
				j1=nearcell1[j][k];
				for (i=0; i<700; i++) {
					if(j1==Cnumarea[i]){
						areaN[k]=areanext1[i]/areacur1[i];
					}
				}
				if(j!=j1){
					if(fabs(areaN[k])>0){
						sumarea=sumarea+fabs(area0-areaN[k])/(fabs(area0)+fabs(areaN[k])); //spatial variability of areal growth
						counterarea=countarea+1;
					}
				}
			}
			//calculate probability density of varea
			if(isnan(varea[j])){
			}
			else{
				if(isinf(varea[j])){
				}
				else{
					bin1_varea=varea[j]/bin_varea;
					cast_varea=(int)bin1_varea;
					p_varea[cast_varea]=p_varea[cast_varea]+1.0;
					count_varea=count_varea+1;
				}
			}
			free(areaN);
			free(areanextN);
		}
//Calculate 95% confidence interval
		double Bin_varea,sum_varea=0;
		for (j=0; j<50000; j++) {
			Bin_varea=(j+0.5)*bin_varea;
			sum_varea=sum_varea+p_varea[j]/(count_varea);
			if(sum_varea>0.95)thre_varea=Bin_varea;
			if(sum_varea>0.95)break;
		}
//Output the result of spatial variability 
		for (j=1; j<10000; j++) {
			if(count3[j]==0 || j==CellNout[t-1])continue;
			for (i=0; i<N[t]; i++) {
				if(sun[k]==Cnum[i]){
							if(varea[j]<thre_varea)fprintf(file3,"%lf %lf %lf %lf\n",x[i],y[i],z[i],varea[j]);
						}
					}
				}
			}
		}
		free(varea);
		free(p_varea);
		for(i=0;i<10000;i++){
			free(nearcell[i]);
			free(nearcell1[i]);
		}
		free(nearcell);
		free(nearcell1);
		free(countnear2);
		for(i=0;i<10000;i++){
			matrix0[i];
			matrix00[i];
		}
		free(matrix0);
		free(matrix00);
		for (i=0; i<10000; i++){
			free(MTx1[i]);
			free(MTy1[i]);
			free(MTz1[i]);
			free(MTnum1[i]);
			free(walllabel[i]);
			free(adjcell1[i]);
			free(adjcell2[i]);
			free(adjcell3[i]);
			free(sun1[i]);
		}
		free(MTx1);
		free(MTy1);
		free(MTz1);
		free(MTnum1);
		free(count3);
		free(walllabel);
		free(adjcell1);
		free(adjcell2);
		free(adjcell3);
		free(sun1);
		fclose(file3);
		free(MTnum);
		for (i=1; i<=16; i++)free(MTnear[i]);
		free(MTnear);
		free(x);
		free(y);
		free(z);
		free(Cnum);
		free(Cnumareacur);
		free(Cnumareanext);
		free(areacur);
		free(areacur1);
		free(areanext);
		free(areanext1);
		free(sun);
		free(parent);
	}
	
	return(0);
}
